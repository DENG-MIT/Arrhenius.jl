# Native equilibrium of one ideal gas with a single initially absent,
# fixed-stoichiometry condensed phase (TP and HP). Include after
# Equilibrium.jl and SurfaceKinetics.jl: the gas-side element-potential system
# and the NASA7/constant-cp standard-thermo coefficient helper are reused.

"""
    StoichiometricCondensedPhase(; name, species, elements, molecular_weight,
        molar_volume, thermo_model, thermo_coefficients, temperature_range,
        reference_pressure=one_atm, charge=0.0)
    StoichiometricCondensedPhase(path)

A pure condensed phase of fixed stoichiometry and constant molar volume, in
equilibrium with one ideal gas through
`equilibrate(gas, phase; ...)`. The phase is a single neutral species with
NASA7 or constant-cp standard thermodynamics; it is not a general multiphase,
nonideal, or electrochemical phase model. All inputs are SI.

Keyword arguments:

- `name`, `species`: phase and species names.
- `elements`: dictionary of element names to atom counts, e.g.
  `Dict("C" => 1.0)` for graphite. Counts must be finite and nonnegative with
  at least one positive; negative (charged) counts are rejected.
- `molecular_weight`: kg/kmol, positive and finite.
- `molar_volume`: m^3/kmol, positive and finite. Enters the incompressible
  pressure correction `h(T, P) = h0(T) + (P - reference_pressure) * molar_volume`
  with `h` in J/kmol.
- `thermo_model`: `"NASA7"` or `"constant-cp"`.
- `thermo_coefficients`: NASA7: `[Tmid, a1..a7 high-temperature region,
  a1..a7 low-temperature region]` giving h in J/kmol and s in J/kmol/K;
  constant-cp: `[T0, h0, s0, cp0]` giving `h = h0 + cp0*(T - T0)` and
  `s = s0 + cp0*log(T/T0)`.
- `temperature_range`: `[Tmin, Tmax]` kelvin declaring polynomial validity.
  NASA7 requires finite positive bounds with `Tmid` strictly inside.
  Constant-cp permits zero and `Inf` as open fit bounds; temperatures must
  still be positive. The JSON representation uses `null` for an infinite upper bound.
- `reference_pressure`: Pa, positive, default one atmosphere.
- `charge`: must be zero; charged species are unsupported.

The path constructor loads the numeric JSON written by
`mechanism/export_condensed_phase.py` (format
`arrhenius-condensed-phase-v1`); the same validations are applied.
"""
struct StoichiometricCondensedPhase
    name::String
    species::String
    elements::Dict{String,Float64}
    molecular_weight::Float64
    molar_volume::Float64
    reference_pressure::Float64
    temperature_range::Vector{Float64}
    # Layout expected by `_surface_species_thermo(phase, 1, T)`: model code per
    # species (1 NASA7, 2 constant-cp) and one 15-wide coefficient row.
    thermo_type::Vector{Int}
    thermo_coefficients::Matrix{Float64}
end

function StoichiometricCondensedPhase(; name, species, elements, molecular_weight,
                                      molar_volume, thermo_model, thermo_coefficients,
                                      temperature_range, reference_pressure=one_atm, charge=0.0)
    phase_name, species_name = String(name), String(species)
    !isempty(phase_name) && !isempty(species_name) ||
        throw(ArgumentError("condensed-phase and species names must be nonempty"))
    charge isa Real && charge == 0 ||
        throw(ArgumentError("charged condensed species are unsupported"))
    elements isa AbstractDict && !isempty(elements) ||
        throw(ArgumentError("condensed-phase elements must be a nonempty dictionary of atom counts"))
    composition = Dict{String,Float64}()
    for (key, value) in elements
        key isa AbstractString ||
            throw(ArgumentError("condensed-phase element names must be strings"))
        value isa Real && isfinite(value) ||
            throw(ArgumentError("condensed-phase element counts must be finite numbers"))
        value >= 0 || throw(ArgumentError("charged or negative element counts are unsupported"))
        value > 0 && (composition[String(key)] = Float64(value))
    end
    !isempty(composition) ||
        throw(ArgumentError("the condensed phase must contain at least one element"))
    molecular_weight isa Real && isfinite(molecular_weight) && molecular_weight > 0 ||
        throw(ArgumentError("a positive finite molecular weight (kg/kmol) is required"))
    molar_volume isa Real && isfinite(molar_volume) && molar_volume > 0 ||
        throw(ArgumentError("a positive finite molar volume (m^3/kmol) is required"))
    reference_pressure isa Real && isfinite(reference_pressure) && reference_pressure > 0 ||
        throw(ArgumentError("a positive finite reference pressure (Pa) is required"))
    temperature_range isa AbstractVector && length(temperature_range) == 2 &&
        all(v -> v isa Real, temperature_range) ||
        throw(ArgumentError("temperature_range must be [Tmin, Tmax] in kelvin"))
    bounds = Float64.(temperature_range)
    isfinite(bounds[1]) && 0 <= bounds[1] < bounds[2] ||
        throw(ArgumentError("condensed-phase temperature bounds must be nonnegative and ordered"))
    thermo_coefficients isa AbstractVecOrMat && all(v -> v isa Real, thermo_coefficients) ||
        throw(ArgumentError("thermo_coefficients must be a numeric vector"))
    coefficients = Float64.(vec(collect(thermo_coefficients)))
    all(isfinite, coefficients) ||
        throw(ArgumentError("condensed-phase thermo coefficients must be finite"))
    model = replace(lowercase(String(thermo_model)), "_" => "-")
    matrix = zeros(1, 15)
    if model == "nasa7"
        type_code = 1
        all(isfinite, bounds) && bounds[1] > 0 ||
            throw(ArgumentError("NASA7 requires finite positive temperature bounds"))
        length(coefficients) == 15 || throw(ArgumentError(
            "NASA7 requires 15 coefficients: [Tmid, a1..a7 high, a1..a7 low]"))
        bounds[1] < coefficients[1] < bounds[2] || throw(ArgumentError(
            "the NASA7 mid temperature must lie strictly inside temperature_range"))
        matrix[1, :] .= coefficients
    elseif model == "constant-cp"
        type_code = 2
        length(coefficients) == 4 || throw(ArgumentError(
            "constant-cp requires 4 coefficients: [T0, h0, s0, cp0]"))
        coefficients[1] > 0 || throw(ArgumentError("constant-cp T0 must be positive"))
        coefficients[4] > 0 || throw(ArgumentError("constant-cp cp0 must be positive"))
        matrix[1, 1:4] .= coefficients
    else
        throw(ArgumentError(
            "unsupported condensed-phase thermo model $thermo_model; " *
            "supported models are NASA7 and constant-cp"))
    end
    return StoichiometricCondensedPhase(phase_name, species_name, composition,
        Float64(molecular_weight), Float64(molar_volume), Float64(reference_pressure),
        bounds, [type_code], matrix)
end

const _CONDENSED_PHASE_FORMAT = "arrhenius-condensed-phase-v1"

# SI unit labels required by the JSON loader; the exporter writes the same map.
const _CONDENSED_PHASE_UNITS = Dict(
    "temperature" => "K", "pressure" => "Pa", "molecular_weight" => "kg/kmol",
    "density" => "kg/m^3", "molar_volume" => "m^3/kmol", "enthalpy" => "J/kmol",
    "entropy" => "J/kmol/K", "heat_capacity" => "J/kmol/K")

function _condensed_json_number(data, key)
    value = get(data, key, nothing)
    value isa Real || throw(ArgumentError("condensed-phase field \"$key\" must be a number"))
    return Float64(value)
end

function StoichiometricCondensedPhase(path::AbstractString)
    data = try
        YAML.load_file(path)
    catch exception
        throw(ArgumentError("could not parse condensed-phase file $path: $exception"))
    end
    data isa AbstractDict ||
        throw(ArgumentError("condensed-phase file must contain a JSON object"))
    get(data, "format", nothing) == _CONDENSED_PHASE_FORMAT || throw(ArgumentError(
        "unsupported condensed-phase format; expected $_CONDENSED_PHASE_FORMAT"))
    units = get(data, "units", nothing)
    units isa AbstractDict || throw(ArgumentError("condensed-phase file must declare its units"))
    for (key, expected) in _CONDENSED_PHASE_UNITS
        get(units, key, nothing) == expected || throw(ArgumentError(
            "unsupported condensed-phase units: expected $expected for $key"))
    end
    raw_elements = get(data, "elements", nothing)
    raw_elements isa AbstractDict || throw(ArgumentError("condensed-phase elements must be a dictionary"))
    elements = Dict{String,Any}()
    for (key, value) in raw_elements
        key isa AbstractString && value isa Real || throw(ArgumentError(
            "condensed-phase element counts must be numbers"))
        elements[String(key)] = value
    end
    thermo = get(data, "thermo", nothing)
    thermo isa AbstractDict || throw(ArgumentError("condensed-phase file must contain a thermo object"))
    model = get(thermo, "model", nothing)
    model isa AbstractString || throw(ArgumentError("condensed-phase thermo model must be a string"))
    raw_coefficients = get(thermo, "coefficients", nothing)
    raw_coefficients isa AbstractVector && all(v -> v isa Real, raw_coefficients) ||
        throw(ArgumentError("condensed-phase thermo coefficients must be a numeric vector"))
    raw_range = get(thermo, "temperature_range_K", nothing)
    raw_range isa AbstractVector && length(raw_range) == 2 && raw_range[1] isa Real &&
        (raw_range[2] isa Real || raw_range[2] === nothing) ||
        throw(ArgumentError("condensed-phase temperature_range_K must be [Tmin, Tmax]"))
    molecular_weight = _condensed_json_number(data, "molecular_weight_kg_per_kmol")
    molar_volume = _condensed_json_number(data, "molar_volume_m3_per_kmol")
    if haskey(data, "density_kg_per_m3")
        density = _condensed_json_number(data, "density_kg_per_m3")
        density > 0 && isapprox(molar_volume, molecular_weight / density; rtol=1e-9) ||
            throw(ArgumentError("molar_volume_m3_per_kmol is inconsistent with molecular weight and density"))
    end
    return StoichiometricCondensedPhase(;
        name=_condensed_json_string(data, "phase_name"),
        species=_condensed_json_string(data, "species_name"),
        elements, molecular_weight, molar_volume,
        thermo_model=model,
        thermo_coefficients=Float64.(raw_coefficients),
        temperature_range=[Float64(raw_range[1]), raw_range[2] === nothing ? Inf : Float64(raw_range[2])],
        reference_pressure=_condensed_json_number(thermo, "reference_pressure_Pa"),
        charge=get(data, "charge", 0.0))
end

function _condensed_json_string(data, key)
    value = get(data, key, nothing)
    value isa AbstractString && !isempty(value) ||
        throw(ArgumentError("condensed-phase field \"$key\" must be a nonempty string"))
    return String(value)
end

"Standard enthalpy (J/kmol) and entropy (J/kmol/K), including the
incompressible pressure correction `h += (P - reference_pressure) * molar_volume`."
function _condensed_standard_hs(phase::StoichiometricCondensedPhase, T, P)
    h, s = _surface_species_thermo(phase, 1, T)
    return h + (P - phase.reference_pressure) * phase.molar_volume, s
end

"Dimensionless condensed-phase chemical potential `h/(R*T) - s/R` at `T`, `P`."
function _condensed_gibbs(phase::StoichiometricCondensedPhase, T, P)
    h, s = _condensed_standard_hs(phase, T, P)
    return h / (R * T) - s / R
end

"""
Caller-owned Float64 workspace for the condensed-phase residual/Jacobian.
`A` is the ne×ns gas element matrix, `As` the ne×np condensed stoichiometry
columns, `g`/`gp` the dimensionless standard chemical potentials (gas:
`h°/RT − s°/R`; condensed: `(h° + (P−P°)vm)/RT − s°/R`), `b` the conserved
element totals (kmol), `logpressure` = log(P/P°). The state is
`u = [λ (ne); ν = log N_g; ξ (np)]`. The public driver admits a single
condensed phase (np = 1); the kernel is written for a general active set.
"""
struct CondensedWorkspace
    A::Matrix{Float64}
    As::Matrix{Float64}
    g::Vector{Float64}
    gp::Vector{Float64}
    logb::Vector{Float64}
    logpressure::Float64
    v::Vector{Float64}
    x::Vector{Float64}
    abar::Vector{Float64}
    G::Vector{Float64}
    total::Vector{Float64}
    residual::Vector{Float64}
    jacobian::Matrix{Float64}
    trial::Vector{Float64}
end

function CondensedWorkspace(A, As, b, g, gp, logpressure)
    ne, ns = size(A)
    size(As, 1) == ne ||
        throw(DimensionMismatch("condensed stoichiometry must share the gas element count"))
    np = size(As, 2)
    length(g) == ns || throw(DimensionMismatch("one gas chemical potential per species required"))
    length(gp) == np || throw(DimensionMismatch("one chemical potential per condensed phase required"))
    length(b) == ne || throw(DimensionMismatch("one element total per element required"))
    all(isfinite, A) && all(isfinite, As) && all(isfinite, g) && all(isfinite, gp) ||
        throw(ArgumentError("stoichiometry and chemical potentials must be finite"))
    all(isfinite, b) && all(>(0), b) ||
        throw(ArgumentError("element totals must be positive and finite (log-form residual)"))
    isfinite(logpressure) || throw(ArgumentError("log pressure must be finite"))
    return CondensedWorkspace(Float64.(A), Float64.(As), Float64.(g), Float64.(gp),
        log.(Float64.(b)), Float64(logpressure), zeros(ns), zeros(ns), zeros(ne),
        zeros(ne), zeros(ne), zeros(ne + 1 + np), zeros(ne + 1 + np, ne + 1 + np),
        zeros(ne + 1 + np))
end

# Shared in-place evaluation of v, x = softmax(v), ā, G and the element totals
# `t = N·A·x + As·ξ` with a stable softmax. Returns false without writing the
# residual when an element total leaves the log domain; callers treat that as
# a rejected trial, never a clipped one.
function _condensed_evaluate!(ws::CondensedWorkspace, u)
    ne, ns = size(ws.A)
    np = size(ws.As, 2)
    λ = view(u, 1:ne)
    ν = u[ne+1]
    ξ = view(u, ne+2:ne+1+np)
    v, x, abar, G, total = ws.v, ws.x, ws.abar, ws.G, ws.total
    mul!(v, transpose(ws.A), λ)
    @. v = v - ws.g - ws.logpressure
    vmax = maximum(v)
    s = 0.0
    @inbounds for k in 1:ns
        w = exp(v[k] - vmax)
        x[k] = w
        s += w
    end
    @. x = x / s
    logsum = vmax + log(s)
    mul!(abar, ws.A, x)
    N = exp(ν)
    @. G = N * abar
    copyto!(total, G)
    np > 0 && mul!(total, ws.As, ξ, 1.0, 1.0)
    return logsum, N, all(>(0), total)
end

"In-place residual `F = [log(t) − log(b); logsumexp(v); Asᵀλ − gp]` (Float64 path)."
function _condensed_residual!(ws::CondensedWorkspace, u)
    ne = size(ws.A, 1)
    np = size(ws.As, 2)
    length(u) == ne + 1 + np || throw(DimensionMismatch("state length must be ne+1+np"))
    logsum, _, valid = _condensed_evaluate!(ws, u)
    F = ws.residual
    valid || return fill!(F, Inf)
    @inbounds for i in 1:ne
        F[i] = log(ws.total[i]) - ws.logb[i]
    end
    F[ne+1] = logsum
    if np > 0
        mul!(view(F, ne+2:ne+1+np), transpose(ws.As), view(u, 1:ne))
        @inbounds for p in 1:np
            F[ne+1+p] -= ws.gp[p]
        end
    end
    return F
end

"""
In-place residual and analytic Jacobian. With `x = softmax(v)`, `N = exp(ν)`,
`ā = A·x`, `G = N·ā`, `t = G + As·ξ`:

    ∂F_i/∂λ_j = N (Σ_k A_ik A_jk x_k − ā_i ā_j) / t_i
    ∂F_i/∂ν   = G_i / t_i
    ∂F_i/∂ξ_p = a_ip / t_i
    ∂F_0/∂λ_j = ā_j      ∂F_0/∂ν = 0      ∂F_0/∂ξ_p = 0
    ∂F_p/∂λ_j = a_pj     all other affinity-row entries 0

Views into the workspace; invalidated by the next call. No state clipping.
"""
function _condensed_residual_jacobian!(ws::CondensedWorkspace, u)
    ne, ns = size(ws.A)
    np = size(ws.As, 2)
    length(u) == ne + 1 + np || throw(DimensionMismatch("state length must be ne+1+np"))
    logsum, N, valid = _condensed_evaluate!(ws, u)
    F = ws.residual
    J = ws.jacobian
    if !valid
        fill!(F, Inf)
        fill!(J, Inf)
        return F, J
    end
    total, abar, G, x = ws.total, ws.abar, ws.G, ws.x
    @inbounds for i in 1:ne
        F[i] = log(total[i]) - ws.logb[i]
    end
    F[ne+1] = logsum
    if np > 0
        mul!(view(F, ne+2:ne+1+np), transpose(ws.As), view(u, 1:ne))
        @inbounds for p in 1:np
            F[ne+1+p] -= ws.gp[p]
        end
    end
    fill!(J, 0.0)
    @inbounds for j in 1:ne, i in 1:ne
        moment = 0.0
        for k in 1:ns
            moment += ws.A[i, k] * ws.A[j, k] * x[k]
        end
        J[i, j] = N * (moment - abar[i] * abar[j]) / total[i]
    end
    @inbounds for i in 1:ne
        J[i, ne+1] = G[i] / total[i]
        J[ne+1, i] = abar[i]
        for p in 1:np
            J[i, ne+1+p] = ws.As[i, p] / total[i]
        end
    end
    @inbounds for p in 1:np, j in 1:ne
        J[ne+1+p, j] = ws.As[j, p]
    end
    return F, J
end

"""
Solver-boundary rejection: throws `ArgumentError` on nonfinite state, negative
condensed amounts, or nonpositive element totals. Newton trials failing this
check must be rejected, never accepted or clipped.
"""
function _condensed_valid_state(ws::CondensedWorkspace, u)
    ne = size(ws.A, 1)
    np = size(ws.As, 2)
    length(u) == ne + 1 + np || throw(DimensionMismatch("state length must be ne+1+np"))
    all(isfinite, u) || throw(ArgumentError("state must be finite"))
    ξ = view(u, ne+2:ne+1+np)
    all(>=(0), ξ) || throw(ArgumentError(
        "condensed amounts must be nonnegative; reject the trial instead of clipping"))
    _condensed_evaluate!(ws, u)
    all(>(0), ws.total) ||
        throw(ArgumentError("element totals must stay positive (log-form residual)"))
    return true
end

"Damped Newton with backtracking; trials with negative condensed amounts or
nonpositive element totals are rejected, never clipped into the domain."
function _condensed_newton!(ws::CondensedWorkspace, u; maxiters=100)
    for iteration in 1:maxiters
        F, J = _condensed_residual_jacobian!(ws, u)
        all(isfinite, F) && all(isfinite, J) || return false
        residual = norm(F)
        norm(F, Inf) < 5e-13 && u[end] >= 0 && return true
        step = try
            -(J \ F)
        catch exception
            exception isa SingularException || rethrow()
            return false
        end
        all(isfinite, step) || return false
        alpha = min(1.0, 20 / max(norm(@view(step[1:end-1]), Inf), 1e-30))
        step[end] < 0 && (alpha = min(alpha, 0.995 * u[end] / (-step[end])))
        alpha > 0 || return false
        accepted = false
        trial = ws.trial
        for backtrack in 1:45
            @. trial = u + alpha * step
            if trial[end] >= 0
                ft = _condensed_residual!(ws, trial)
                if all(isfinite, ft) && norm(ft) < residual
                    copyto!(u, trial)
                    accepted = true
                    break
                end
            end
            alpha *= 0.5
        end
        accepted || return false
    end
    return false
end

"""
Map the phase element dictionary onto the selected independent element rows of
the gas equilibrium system. Returns `(a, inactive)`: `a` is the condensed
stoichiometry in the reduced rows, and `inactive` is true when some phase
element has zero initial inventory, so conservation excludes formation.
Elements absent from the gas mechanism, or compositions inconsistent with
element rows removed by QR, are rejected.
"""
function _condensed_element_map(gas::Solution, system, phase::StoichiometricCondensedPhase)
    ne = size(system.A, 1)
    b_all = gas.ele_matrix * system.X
    selected = map(axes(system.A, 1)) do j
        findfirst(i -> gas.ele_matrix[i, system.species] == system.A[j, :],
                  axes(gas.ele_matrix, 1))
    end
    any(isnothing, selected) &&
        error("could not recover the original element rows of the equilibrium system")
    for name in keys(phase.elements)
        isnothing(findfirst(==(name), gas.elements)) && throw(ArgumentError(
            "condensed-phase element $name is absent from the gas mechanism"))
    end
    a = zeros(ne)
    for (name, count) in phase.elements
        e = findfirst(==(name), gas.elements)
        j = findfirst(==(e), selected)
        if isnothing(j)
            b_all[e] > 0 || return a, true
        else
            a[j] = count
        end
    end
    # Rows removed by QR are exact combinations of the selected rows. Any
    # positive condensed amount must satisfy the same combination, or an
    # untracked element balance would break.
    for e in eachindex(gas.elements)
        (e in selected || b_all[e] == 0) && continue
        count = get(phase.elements, gas.elements[e], 0.0)
        row = gas.ele_matrix[e, system.species]
        c = system.At \ row
        norm(system.At * c - row, Inf) <= 1e-10 * max(1.0, norm(row, Inf)) ||
            error("dependent element row is not spanned by the selected rows")
        isapprox(count, dot(c, a); atol=1e-10, rtol=1e-10) || throw(ArgumentError(
            "the condensed-phase composition is inconsistent with the dependent " *
            "element rows of the gas; no positive amount can conserve all elements"))
    end
    return a, false
end

"Relative all-row element conservation error (including rows removed by QR);
throws above the 2e-9 convention shared with the gas-only solver."
function _condensed_element_check(gas::Solution, phase::StoichiometricCondensedPhase,
                                  X, gas_species_moles, condensed_moles)
    before = gas.ele_matrix * X
    after = gas.ele_matrix * gas_species_moles
    for (name, count) in phase.elements
        e = findfirst(==(name), gas.elements)
        isnothing(e) || (after[e] += count * condensed_moles)
    end
    error_elements = maximum(abs.(after - before)) / maximum(abs, before)
    error_elements < 2e-9 ||
        error("condensed-phase equilibrium failed elemental conservation: $error_elements")
    return error_elements
end

"""
TP equilibrium of the gas plus the initially absent condensed phase at `T`,
`P`, per 1 kmol of initial gas. The gas-only Newton runs first; the phase is
admitted only when its affinity `aᵀλ − g_p` is positive (present ⇒ affinity 0,
absent ⇒ affinity ≤ 0, amounts nonnegative without projection). On Newton
failure a homotopy in the condensed chemical potential changes only the
starting path; the final equations keep the supplied phase free energy.
"""
function _condensed_tp!(system, phase::StoichiometricCondensedPhase, a, inactive, T, P)
    gas = system.gas
    gas_state = _equilibrium_at!(system, T, P)
    h, s = _condensed_standard_hs(phase, T, P)
    gp = h / (R * T) - s / R
    ne = size(system.A, 1)
    N = exp(system.state[end])
    ξ = 0.0
    affinity = NaN
    status = inactive ? :inactive_by_elements : :gas_only
    gas_species_per = N .* gas_state.X
    if !inactive
        λ = view(system.state, 1:ne)
        affinity = dot(a, λ) - gp
        if affinity > 1e-11
            ws = CondensedWorkspace(system.A, reshape(a, ne, 1), system.b, system.g,
                                    [gp], log(P / one_atm))
            u = vcat(system.state, 0.0)
            if !_condensed_newton!(ws, u)
                seed = vcat(system.state, 0.0)
                starting_g = dot(a, λ)
                converged = false
                for increment in (2.0, 0.25)
                    copyto!(u, seed)
                    stages = max(2, ceil(Int, abs(starting_g - gp) / increment) + 1)
                    converged = true
                    for fraction in range(0.0, 1.0; length=stages)
                        ws.gp[1] = (1 - fraction) * starting_g + fraction * gp
                        if !_condensed_newton!(ws, u)
                            converged = false
                            break
                        end
                    end
                    converged && break
                end
                converged || error(
                    "condensed-phase equilibrium did not converge at $T K and $P Pa " *
                    "(affinity $affinity)")
            end
            ws.gp[1] = gp
            F, = _condensed_residual_jacobian!(ws, u)
            norm(F, Inf) < 1e-10 && u[end] >= 0 ||
                error("invalid active condensed-phase root at $T K and $P Pa")
            ξ = u[end]
            N = exp(u[ne+1])
            gas_species_per = zeros(gas.n_species)
            gas_species_per[system.species] .= N .* ws.x
            affinity = F[end]
            status = :condensed_present
        end
        affinity <= 1e-10 && abs(ξ * affinity) < 1e-10 ||
            error("condensed-phase equilibrium failed phase complementarity")
    end
    element_error = _condensed_element_check(gas, phase, system.X, gas_species_per, ξ)
    hgas = cal_h_RT(gas, T, P, system.X) .* (R * T)
    H = dot(gas_species_per, hgas) + ξ * h
    return (T=Float64(T), P=Float64(P), N=N, ξ=ξ, gas_species_per=gas_species_per,
            H=H, status=status, affinity=affinity, element_error=element_error)
end

"Scale the per-kmol state by `initial_gas_moles` and build the public result."
function _condensed_assemble(gas::Solution, state, enthalpy_error, initial_gas_moles)
    gas_species_moles = state.gas_species_per .* initial_gas_moles
    gas_moles = state.N * initial_gas_moles
    condensed_moles = state.ξ * initial_gas_moles
    X = state.gas_species_per ./ state.N
    Y = X .* gas.MW ./ dot(X, gas.MW)
    return (T=state.T, P=state.P, X=X, Y=Y, gas_moles=gas_moles,
            condensed_moles=condensed_moles, gas_species_moles=gas_species_moles,
            species_moles=vcat(gas_species_moles, condensed_moles),
            status=state.status, affinity=state.affinity,
            element_error=state.element_error, enthalpy_error=enthalpy_error)
end

"HP outer solve: bracketed safeguarded secant on the conserved enthalpy,
including the condensed-phase enthalpy (with its molar-volume pressure term)."
function _condensed_hp_result(system, phase, a, inactive, T0, P0, lo, hi, rtol,
                              initial_gas_moles)
    gas = system.gas
    target = dot(system.X, cal_h_RT(gas, T0, P0, system.X)) * (R * T0)
    scale = max(abs(target), 1e6 * dot(gas.MW, system.X))
    state(temp) = _condensed_tp!(system, phase, a, inactive, temp, P0)
    assemble(trial) = _condensed_assemble(gas, trial, abs(trial.H - target) / scale,
                                          initial_gas_moles)
    trialT = clamp(max(T0, 3500.0), lo, hi)
    trial = state(trialT)
    residual = trial.H - target
    abs(residual) < rtol * scale && return assemble(trial)
    if residual < 0
        lower, fl = trialT, residual
        upper, fu = trialT, residual
        while fu < 0 && upper < hi
            upper = min(hi, max(upper + 100, 1.4 * upper))
            trial = state(upper)
            fu = trial.H - target
        end
    else
        upper, fu = trialT, residual
        lower, fl = trialT, residual
        while fl > 0 && lower > lo
            lower = max(lo, min(lower - 100, lower / 1.4))
            trial = state(lower)
            fl = trial.H - target
        end
    end
    fl <= 0 <= fu || throw(ArgumentError("HP equilibrium is not bracketed in $(lo)–$(hi) K"))
    for iteration in 1:80
        width = upper - lower
        trialT = clamp((lower * fu - upper * fl) / (fu - fl),
                       lower + 0.05 * width, upper - 0.05 * width)
        trial = state(trialT)
        residual = trial.H - target
        (abs(residual) <= rtol * scale || width < min(1e-7, rtol * trialT)) &&
            return assemble(trial)
        if residual > 0
            upper, fu = trialT, residual
        else
            lower, fl = trialT, residual
        end
    end
    error("condensed-phase HP temperature search did not converge")
end

"""
    equilibrate(gas::Solution, phase::StoichiometricCondensedPhase; T,
        P=one_atm, X, mode=:TP, initial_gas_moles=1.0,
        temperature_bounds=(200.0, 6000.0), property_rtol=1e-10)

Native equilibrium of one ideal gas with a single *initially absent*
fixed-stoichiometry condensed phase. Supported conserved pairs are `:TP` and
`:HP` (symbols or strings); the gas-only `equilibrate` methods and their
tolerances are unchanged. The initial state is `initial_gas_moles` kmol of gas
at `T`, `P`, `X` plus zero kmol of the condensed phase; initially present
solids are not supported, and neither are multiple condensed phases, nonideal
phases, or charged species.

The condensed phase forms only when its element inventory allows it and its
formation affinity is positive. If the phase requires an element with zero
initial inventory, formation is excluded by conservation: the gas-only state
is returned with zero condensed amount, `status == :inactive_by_elements`,
and `affinity == NaN`.

Returns `(; T, P, X, Y, gas_moles, condensed_moles, gas_species_moles,
species_moles, status, affinity, element_error, enthalpy_error)`: gas mole and
mass fractions `X`/`Y`, the total gas amount `gas_moles` (kmol), the condensed
amount `condensed_moles` (kmol), the per-species gas amounts, the combined
amounts (gas species order, then the condensed species), a `status` of
`:gas_only`, `:condensed_present`, or `:inactive_by_elements`, the
dimensionless phase affinity (zero when present, nonpositive when absent), the
relative all-element conservation error, and the relative enthalpy
conservation error (`NaN` for `:TP`).

The system is solved per 1 kmol of initial gas and all amounts are scaled by
`initial_gas_moles`; equilibrium temperatures and compositions do not depend
on it. `temperature_bounds` brackets the `:HP` temperature search; like the
gas-only solver, evaluation outside the species fit ranges extrapolates the
polynomials. `property_rtol` controls the `:HP` conserved-enthalpy solve.
"""
function equilibrate(gas::Solution, phase::StoichiometricCondensedPhase;
                     T, P=one_atm, X, mode=:TP, initial_gas_moles=1.0,
                     temperature_bounds=(200.0, 6000.0), property_rtol=1e-10)
    gas.thermo isa IdealGasThermo || throw(ArgumentError("ideal-gas thermo required"))
    isfinite(T) && T > 0 && isfinite(P) && P > 0 ||
        throw(ArgumentError("positive finite T and P required"))
    mode = Symbol(mode)
    mode in (:TP, :HP) || throw(ArgumentError(
        "unsupported condensed-phase equilibrium mode: $mode (supported: :TP, :HP)"))
    isfinite(property_rtol) && 0 < property_rtol < 1 ||
        throw(ArgumentError("property_rtol must lie between zero and one"))
    length(temperature_bounds) == 2 ||
        throw(ArgumentError("provide lower and upper temperature bounds"))
    lo, hi = Float64.(temperature_bounds)
    isfinite(lo) && isfinite(hi) && 0 < lo < hi ||
        throw(ArgumentError("temperature bounds must be finite, positive, and ordered"))
    initial_gas_moles isa Real && isfinite(initial_gas_moles) && initial_gas_moles > 0 ||
        throw(ArgumentError("initial_gas_moles must be positive and finite"))
    initial = mole_fractions(gas, X)
    system = _equilibrium_system(gas, initial)
    a, inactive = _condensed_element_map(gas, system, phase)
    T0, P0 = Float64(T), Float64(P)
    if mode == :TP
        state = _condensed_tp!(system, phase, a, inactive, T0, P0)
        return _condensed_assemble(gas, state, NaN, Float64(initial_gas_moles))
    end
    return _condensed_hp_result(system, phase, a, inactive, T0, P0, lo, hi,
                                property_rtol, Float64(initial_gas_moles))
end

export StoichiometricCondensedPhase
