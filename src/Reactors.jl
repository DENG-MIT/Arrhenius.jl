"""
    IdealGasReactor(gas; temperature, pressure=one_atm,
                    mole_fractions=nothing, mass_fractions=nothing,
                    constraint=:constant_pressure, energy=:adiabatic,
                    rate_multipliers=nothing)

A homogeneous, closed ideal-gas reactor. Supply exactly one composition as a
species-name dictionary or a vector in mechanism order. Amounts are normalized.
`constraint` is `:constant_pressure` or `:constant_volume`; `energy` is
`:adiabatic` or `:isothermal`. SI units are used throughout (K, Pa, kg, s).

The state is `[Y₁, …, Yₙ, T]`. Constant-pressure reactors conserve specific
enthalpy when adiabatic; constant-volume reactors conserve specific internal
energy. Isothermal reactors impose the initial temperature. Walls, inlet and
outlet flows, surface chemistry, and non-ideal equations of state are excluded.
"""
struct IdealGasReactor{G,T,M}
    gas::G
    temperature::T
    pressure::T
    density::T
    mass_fractions::Vector{T}
    constraint::Symbol
    energy::Symbol
    rate_multipliers::M
end

function _reactor_composition(gas, composition, ::Type{T}) where {T}
    values = zeros(T, gas.n_species)
    if composition isa AbstractDict
        for (name, value) in composition
            index = findfirst(==(String(name)), gas.species_names)
            isnothing(index) && throw(ArgumentError("unknown species: $name"))
            values[index] = value
        end
    elseif composition isa AbstractVector
        length(composition) == gas.n_species ||
            throw(DimensionMismatch("one composition entry per species required"))
        copyto!(values, composition)
    else
        throw(ArgumentError("composition must be a species dictionary or vector"))
    end
    total = sum(values)
    all(isfinite, values) && all(>=(zero(T)), values) && isfinite(total) && total > 0 ||
        throw(ArgumentError("composition must be finite, nonnegative, and nonzero"))
    values ./= total
    return values
end

function IdealGasReactor(gas::Solution; temperature, pressure=one_atm,
                         mole_fractions=nothing, mass_fractions=nothing,
                         constraint=:constant_pressure, energy=:adiabatic,
                         rate_multipliers=nothing)
    gas.thermo isa IdealGasThermo || throw(ArgumentError("ideal-gas thermo required"))
    constraint in (:constant_pressure, :constant_volume) ||
        throw(ArgumentError("constraint must be :constant_pressure or :constant_volume"))
    energy in (:adiabatic, :isothermal) ||
        throw(ArgumentError("energy must be :adiabatic or :isothermal"))
    isfinite(temperature) && temperature > 0 && isfinite(pressure) && pressure > 0 ||
        throw(ArgumentError("temperature and pressure must be positive and finite"))
    xor(isnothing(mole_fractions), isnothing(mass_fractions)) ||
        throw(ArgumentError("supply exactly one of mole_fractions and mass_fractions"))
    T = promote_type(eltype(gas.MW), typeof(float(temperature)), typeof(float(pressure)))
    if isnothing(mass_fractions)
        X = _reactor_composition(gas, mole_fractions, T)
        Y = X .* gas.MW
        Y ./= sum(Y)
    else
        Y = _reactor_composition(gas, mass_fractions, T)
    end
    inverse_mw = sum(Y ./ gas.MW)
    density = T(pressure) / (T(R) * T(temperature) * inverse_mw)
    multipliers = if isnothing(rate_multipliers)
        nothing
    else
        length(rate_multipliers) == gas.n_reactions ||
            throw(DimensionMismatch("one rate multiplier per reaction required"))
        values = T.(rate_multipliers)
        all(isfinite, values) && all(>=(zero(T)), values) ||
            throw(ArgumentError("rate multipliers must be finite and nonnegative"))
        values
    end
    return IdealGasReactor(gas, T(temperature), T(pressure), density, Y,
                           constraint, energy, multipliers)
end

"Return an independent initial state vector `[mass_fractions; temperature]`."
reactor_state(reactor::IdealGasReactor) = vcat(reactor.mass_fractions, reactor.temperature)

struct ReactorWorkspace{T}
    X::Vector{T}
    C::Vector{T}
    cp_R::Vector{T}
    h_mole::Vector{T}
    entropy::Vector{T}
    wdot::Vector{T}
    kinetics::KineticsWorkspace{T}
end

function ReactorWorkspace(gas::Solution, ::Type{T}) where {T}
    return ReactorWorkspace(ntuple(_ -> zeros(T, gas.n_species), 6)...,
                            KineticsWorkspace(gas.reaction, T))
end

"A reusable in-place ODE function. Use a separate instance for each concurrent solve."
struct ReactorRHS{R,W,T}
    reactor::R
    workspace::W
    jac_state::Vector{T}
    jac_base::Vector{T}
    jac_plus::Vector{T}
    jac_minus::Vector{T}
end

"""
    reactor_rhs(reactor)

Prepare a callable `rhs!(du, u, p, t)` using native Arrhenius kinetics and
thermodynamics. `p` is unused. The workspace has the same scalar type as the
initial state; use `reactor_jacobian!` with solvers requiring derivatives.

Negative trial mass fractions are clipped only when evaluating reaction
concentrations. The state and its mass/element balances are never projected.
"""
function reactor_rhs(reactor::IdealGasReactor)
    u0 = reactor_state(reactor)
    return ReactorRHS(reactor, ReactorWorkspace(reactor.gas, eltype(u0)),
                      copy(u0), zero(u0), zero(u0), zero(u0))
end

function _reactor_tpρ(reactor, u)
    ns = reactor.gas.n_species
    length(u) == ns + 1 || throw(DimensionMismatch("state must contain n_species + 1 entries"))
    temperature = reactor.energy === :isothermal ? oftype(u[end], reactor.temperature) : u[end]
    inverse_mw = zero(temperature)
    @inbounds for k in 1:ns
        inverse_mw += u[k] / reactor.gas.MW[k]
    end
    isfinite(temperature) && temperature > 0 && isfinite(inverse_mw) && inverse_mw > 0 ||
        throw(DomainError(temperature, "reactor temperature and inverse molecular weight must be positive and finite"))
    if reactor.constraint === :constant_pressure
        pressure = oftype(temperature, reactor.pressure)
        density = pressure / (oftype(temperature, R) * temperature * inverse_mw)
    else
        density = oftype(temperature, reactor.density)
        pressure = density * oftype(temperature, R) * temperature * inverse_mw
    end
    return temperature, pressure, density, inverse_mw
end

function (rhs::ReactorRHS)(du, u, p, t)
    reactor, workspace = rhs.reactor, rhs.workspace
    gas = reactor.gas
    ns = gas.n_species
    length(du) == ns + 1 || throw(DimensionMismatch("derivative must match reactor state"))
    temperature, pressure, density, inverse_mw = _reactor_tpρ(reactor, u)
    gas_constant = oftype(temperature, R)
    @inbounds for k in 1:ns
        workspace.X[k] = u[k] / (gas.MW[k] * inverse_mw)
        workspace.C[k] = max(u[k], zero(u[k])) * density / gas.MW[k]
    end
    cal_cp_R!(workspace.cp_R, gas, temperature, pressure, workspace.X)
    cal_h_RT!(workspace.h_mole, gas, temperature, pressure, workspace.X)
    cal_s0_R!(workspace.entropy, gas, temperature, pressure, workspace.X)
    @inbounds for k in 1:ns
        workspace.h_mole[k] *= gas_constant * temperature
        workspace.entropy[k] *= gas_constant
    end
    wdot!(workspace.wdot, gas.reaction, temperature, workspace.C,
          workspace.entropy, workspace.h_mole, workspace.kinetics;
          rate_multipliers=reactor.rate_multipliers)
    cv_offset = reactor.constraint === :constant_volume ? one(temperature) : zero(temperature)
    capacity, energy_source = zero(temperature), zero(temperature)
    @inbounds for k in 1:ns
        du[k] = workspace.wdot[k] * gas.MW[k] / density
        capacity += u[k] * (workspace.cp_R[k] - cv_offset) / gas.MW[k]
        energy_source += workspace.wdot[k] *
                         (workspace.h_mole[k] - cv_offset * gas_constant * temperature)
    end
    capacity *= gas_constant
    isfinite(capacity) && capacity > 0 || throw(DomainError(capacity, "positive heat capacity required"))
    du[end] = reactor.energy === :isothermal ? zero(temperature) :
              -energy_source / (density * capacity)
    return nothing
end

"""
    reactor_jacobian!(J, u, rhs, t=0)

Evaluate a dense finite-difference Jacobian of a prepared `ReactorRHS`.
Second-order central differences are used inside the physical domain and
second-order forward differences at zero mass fractions. Negative trial mass
fractions use second-order backward differences, sampling the RHS at
`u[j] - step` and `u[j] - 2 * step` so no positive perturbation crosses the
clipped concentration branch of the RHS. This callback avoids passing
automatic-differentiation scalars through the fixed-type workspace.
The input state is unchanged.
"""
function reactor_jacobian!(J, u, rhs::ReactorRHS, t=0)
    n = length(u)
    size(J) == (n, n) || throw(DimensionMismatch("Jacobian must have state dimensions"))
    copyto!(rhs.jac_state, u)
    rhs(rhs.jac_base, u, nothing, t)
    relative_step = cbrt(eps(eltype(rhs.jac_state)))
    @inbounds for j in 1:n
        if j == n && rhs.reactor.energy === :isothermal
            for i in 1:n
                J[i, j] = zero(eltype(J))
            end
            continue
        end
        scale = j == n ? one(u[j]) : oftype(u[j], 1e-6)
        step = relative_step * max(abs(u[j]), scale)
        rhs.jac_state[j] = u[j] + step
        step = rhs.jac_state[j] - u[j]
        if j < n && u[j] < 0
            rhs.jac_state[j] = u[j] - step
            rhs(rhs.jac_plus, rhs.jac_state, nothing, t)
            rhs.jac_state[j] = u[j] - 2 * step
            rhs(rhs.jac_minus, rhs.jac_state, nothing, t)
            for i in 1:n
                J[i, j] = (3 * rhs.jac_base[i] - 4 * rhs.jac_plus[i] +
                           rhs.jac_minus[i]) / (2 * step)
            end
            rhs.jac_state[j] = u[j]
            continue
        end
        rhs(rhs.jac_plus, rhs.jac_state, nothing, t)
        if u[j] >= step
            rhs.jac_state[j] = u[j] - step
            rhs(rhs.jac_minus, rhs.jac_state, nothing, t)
            for i in 1:n
                J[i, j] = (rhs.jac_plus[i] - rhs.jac_minus[i]) / (2step)
            end
        else
            rhs.jac_state[j] = u[j] + 2step
            rhs(rhs.jac_minus, rhs.jac_state, nothing, t)
            for i in 1:n
                J[i, j] = (-3rhs.jac_base[i] + 4rhs.jac_plus[i] - rhs.jac_minus[i]) / (2step)
            end
        end
        rhs.jac_state[j] = u[j]
    end
    return nothing
end

"""
    reactor_properties(reactor, u=reactor_state(reactor))

Return temperature (K), pressure (Pa), density (kg/m³), composition, specific
heat capacities (J/kg/K), enthalpy and internal energy (J/kg), mass-fraction sum,
and elemental inventories (kmol of element/kg). Values use the unprojected
state, so mass, element, and energy drift can be checked after integration.
"""
function reactor_properties(reactor::IdealGasReactor, u=reactor_state(reactor))
    gas = reactor.gas
    temperature, pressure, density, inverse_mw = _reactor_tpρ(reactor, u)
    Y = copy(@view u[1:gas.n_species])
    mole_per_mass = Y ./ gas.MW
    X = mole_per_mass ./ inverse_mw
    h_mole = cal_h(gas, temperature, pressure, X)
    cp_mole = cal_cp(gas, temperature, pressure, X)
    enthalpy = dot(mole_per_mass, h_mole)
    cp = dot(mole_per_mass, cp_mole)
    return (; temperature, pressure, density, mole_fractions=X, mass_fractions=Y,
             cp, cv=cp - R * inverse_mw, enthalpy,
             internal_energy=enthalpy - R * temperature * inverse_mw,
             mass_fraction_sum=sum(Y), elemental_inventory=gas.ele_matrix * mole_per_mass)
end

"""
    reactor_problem(reactor, tspan)

Prepare a solver-independent named tuple `(f, jac, tgrad, u0, tspan, p)`.
The callbacks use the SciML in-place calling convention. `tgrad` is zero for
this autonomous system. Construct a solver's ODE problem from these fields.
Each call allocates independent workspaces for safe concurrent integrations.
"""
function reactor_problem(reactor::IdealGasReactor, tspan)
    length(tspan) == 2 || throw(ArgumentError("tspan must contain start and end times"))
    all(isfinite, tspan) && tspan[2] > tspan[1] ||
        throw(ArgumentError("tspan must be finite and strictly increasing"))
    rhs = reactor_rhs(reactor)
    jac = (J, u, p, t) -> reactor_jacobian!(J, u, rhs, t)
    tgrad = (du, u, p, t) -> (fill!(du, zero(eltype(du))); nothing)
    return (f=rhs, jac=jac, tgrad=tgrad, u0=reactor_state(reactor),
            tspan=(float(tspan[1]), float(tspan[2])), p=nothing)
end

"""
    solve_reactor(reactor, tspan; integrator, kwargs...)

Integrate with a caller-supplied `integrator(problem; kwargs...)`, where
`problem` is the named tuple from `reactor_problem`. Returns the integrator's
solution unchanged. The core package has no mandatory ODE solver dependency.

For example, after loading `SciMLBase` and `OrdinaryDiffEqBDF`:
```julia
function bdf(problem; kwargs...)
    f = ODEFunction(problem.f; jac=problem.jac, tgrad=problem.tgrad)
    ode = ODEProblem(f, problem.u0, problem.tspan, problem.p)
    solve(ode, QNDF(); kwargs...)
end
solution = solve_reactor(reactor, (0.0, 0.001); integrator=bdf,
                         reltol=1e-8, abstol=1e-14, saveat=1e-5)
```
"""
function solve_reactor(reactor::IdealGasReactor, tspan; integrator, kwargs...)
    return integrator(reactor_problem(reactor, tspan); kwargs...)
end

export IdealGasReactor, ReactorRHS, reactor_state, reactor_rhs
export reactor_jacobian!, reactor_properties, reactor_problem, solve_reactor
