# Native fixed-grid two-term EEDF solver.
#
# The numerical scheme in this file is translated from Cantera's experimental
# EEDFTwoTermApproximation (Cantera 4.0 development source).
#
# Copyright (c) 2001-2009, California Institute of Technology
# All rights reserved.
#
# Copyright (c) 2009 Sandia Corporation. Under the terms of
# Contract AC04-94AL85000 with Sandia Corporation, the U.S. Government
# retains certain rights in this software.
#
# Copyright (c) 2011-2026, Cantera Developers
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
#
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# - Neither the name of the California Institute of Technology, Sandia
#   Corporation nor the names of other contributors may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

const _EEDF_ELECTRON_CHARGE = 1.602176634e-19       # C (J/eV)
const _EEDF_ELECTRON_MASS = 9.1093837015e-31        # kg (pinned Cantera value)
const _EEDF_BOLTZMANN = 1.380649e-23                # J/K
const _EEDF_AVOGADRO_KMOL = 6.02214076e26           # 1/kmol
const _EEDF_GAMMA = sqrt(2 * _EEDF_ELECTRON_CHARGE / _EEDF_ELECTRON_MASS)

"""
    EEDFState(model; T, P, mole_fractions, molecular_weights, reduced_field, frequency=0, number_density=nothing)

Gas conditions for a temporal two-term electron-energy calculation. Temperature
is K, pressure Pa, molecular weights kg/kmol, reduced electric field V·m², and
frequency Hz. Dictionaries are keyed by collision-target names. Mole fractions
are normalized over participating targets; additional species are ignored.
`number_density` is the total gas number density in m⁻³. It defaults to
`P/(kB*T)`; supply the actual density for a two-temperature or non-ideal gas.
"""
struct EEDFState
    T::Float64
    P::Float64
    mole_fractions::Dict{String,Float64}
    molecular_weights::Dict{String,Float64}
    reduced_field::Float64
    frequency::Float64
    number_density::Float64
end

# Preserve positional construction with the original ideal-gas density.
EEDFState(T, P, x, mw, reduced_field, frequency) =
    EEDFState(T, P, x, mw, reduced_field, frequency, Float64(P)/(_EEDF_BOLTZMANN*Float64(T)))

function EEDFState(model::EEDFModel; T, P, mole_fractions::AbstractDict,
                   molecular_weights::AbstractDict, reduced_field, frequency=0,
                   number_density=nothing)
    return _eedf_state(model,T,P,mole_fractions,molecular_weights,reduced_field,
        frequency,number_density,false)
end

function _accepted_eedf_state(model::EEDFModel;T,P,mole_fractions::AbstractDict,
        molecular_weights::AbstractDict,reduced_field,frequency=0,number_density=nothing)
    return _eedf_state(model,T,P,mole_fractions,molecular_weights,reduced_field,
        frequency,number_density,true)
end

function _eedf_state(model,T,P,mole_fractions,molecular_weights,reduced_field,
        frequency,number_density,accepted_signed)
    Tf = Float64(T)
    Pf = Float64(P)
    EN = Float64(reduced_field)
    freq = Float64(frequency)
    isfinite(Tf) && Tf > 0 || throw(ArgumentError("EEDF temperature must be finite and positive"))
    isfinite(Pf) && Pf > 0 || throw(ArgumentError("EEDF pressure must be finite and positive"))
    isfinite(EN) && EN >= 0 || throw(ArgumentError("reduced electric field must be finite and nonnegative"))
    isfinite(freq) && freq >= 0 || throw(ArgumentError("electric field frequency must be finite and nonnegative"))

    density = number_density === nothing ? Pf/(_EEDF_BOLTZMANN*Tf) : Float64(number_density)
    isfinite(density) && density > 0 || throw(ArgumentError("EEDF number density must be finite and positive m⁻³"))

    x = Dict{String,Float64}()
    mw = Dict{String,Float64}()
    total = 0.0
    for target in model.target_names
        haskey(mole_fractions, target) || throw(ArgumentError("missing mole fraction for collision target '$target'"))
        haskey(molecular_weights, target) || throw(ArgumentError("missing molecular weight for collision target '$target'"))
        xi = Float64(mole_fractions[target])
        wi = Float64(molecular_weights[target])
        isfinite(xi) && (accepted_signed || xi >= 0) ||
            throw(ArgumentError("target mole fractions must be finite and nonnegative"))
        isfinite(wi) && wi > 0 || throw(ArgumentError("target molecular weights must be finite and positive kg/kmol"))
        x[target] = xi
        mw[target] = wi
        total += xi
    end
    isfinite(total) && total > 0 || throw(ArgumentError("collision-target mole fractions must have a positive sum"))
    for target in keys(x)
        x[target] /= total
    end
    return EEDFState(Tf, Pf, x, mw, EN, freq, density)
end

"""
    TwoTermOptions(; delta0=1e14, max_iterations=200, factor=4., rtol=1e-5,
                   initial_kTe=2., positivity_floor=1e-300, low_field_threshold=1e-21)

Iteration controls for [`solve_eedf`](@ref). `initial_kTe` is eV and the field
threshold is V·m². The positivity floor regularizes logarithmic slopes only.
"""
Base.@kwdef struct TwoTermOptions
    delta0::Float64 = 1e14
    max_iterations::Int = 200
    factor::Float64 = 4.0
    rtol::Float64 = 1e-5
    initial_kTe::Float64 = 2.0
    positivity_floor::Float64 = 1e-300
    low_field_threshold::Float64 = 1e-21
end

"""
    EEDFResult

Electron-energy distribution on `edges` and `centers` (eV), electron `mobility`
(m²/(V·s)), and convergence diagnostics. The center distribution satisfies
∫√ε f(ε)dε = 1; edge values are interpolated with endpoint holding and are not
renormalized. `errors` and `deltas` record the iteration history.
"""
struct EEDFResult
    edges::Vector{Float64}
    edge_eedf::Vector{Float64}
    centers::Vector{Float64}
    center_eedf::Vector{Float64}
    mobility::Float64
    iterations::Int
    converged::Bool
    errors::Vector{Float64}
    deltas::Vector{Float64}
end

struct _EEDFInterval
    i::Int
    j::Int
    eps_a::Float64
    eps_b::Float64
    sigma_a::Float64
    sigma_b::Float64
end

struct _EEDFCollisionCache
    collision::ElectronCollision
    target_fraction::Float64
    incoming_factor::Float64
    intervals::Vector{_EEDFInterval}
end

mutable struct EEDFWorkspace
    centers::Vector{Float64}
    collision_cache::Vector{_EEDFCollisionCache}
    total_center::Vector{Float64}
    total_edge::Vector{Float64}
    elastic_edge::Vector{Float64}
    g::Vector{Float64}
    pq::Matrix{Float64}
    operator::Matrix{Float64}
    system::Matrix{Float64}
end

"""Cantera-compatible composite Simpson quadrature on a strictly increasing grid."""
function _eedf_simpson(f::AbstractVector, x::AbstractVector)
    length(f) == length(x) || throw(DimensionMismatch("quadrature vectors must have equal lengths"))
    n = length(f)
    n >= 2 || throw(ArgumentError("quadrature needs at least two points"))
    all(isfinite, f) && all(isfinite, x) || throw(ArgumentError("quadrature data must be finite"))
    all(diff(x) .> 0) || throw(ArgumentError("quadrature grid must be strictly increasing"))
    last_simpson = isodd(n) ? n : n - 1
    value = 0.0
    for i in 2:2:last_simpson-1
        h0 = x[i] - x[i-1]
        h1 = x[i+1] - x[i]
        hp = h0 + h1
        value += hp / 6 * ((2 - h1/h0) * f[i-1] +
                          hp^2/(h0*h1) * f[i] + (2 - h0/h1) * f[i+1])
    end
    if iseven(n)
        value += (x[end] - x[end-1]) * (f[end] + f[end-1]) / 2
    end
    return value
end

function _linear_interp_hold(x::Float64, xp::Vector{Float64}, yp::Vector{Float64})
    x <= xp[1] && return yp[1]
    x >= xp[end] && return yp[end]
    hi = searchsortedfirst(xp, x)
    xp[hi] == x && return yp[hi]
    lo = hi - 1
    return yp[lo] + (yp[hi] - yp[lo]) * (x - xp[lo]) / (xp[hi] - xp[lo])
end

function _validate_model(model::EEDFModel)
    edges = model.energy_edges
    length(edges) >= 3 || throw(ArgumentError("EEDF grid must contain at least two cells"))
    all(isfinite, edges) && all(edges .>= 0) || throw(ArgumentError("EEDF grid must be finite and nonnegative"))
    all(diff(edges) .> 0) || throw(ArgumentError("EEDF grid must be strictly increasing"))
    isempty(model.target_names) && throw(ArgumentError("EEDF model has no collision targets"))
    isempty(model.collisions) && throw(ArgumentError("EEDF model has no collisions"))
    length(unique(model.target_names)) == length(model.target_names) || throw(ArgumentError("collision target names must be unique"))
    target_set = Set(model.target_names)
    elastic_seen = Set{String}()
    for (k, collision) in pairs(model.collisions)
        collision.target in target_set || throw(ArgumentError("collision $k has unknown target '$(collision.target)'"))
        length(collision.energy) == length(collision.cross_section) || throw(ArgumentError("collision $k table lengths differ"))
        length(collision.energy) >= 2 || throw(ArgumentError("collision $k table needs at least two points"))
        all(isfinite, collision.energy) && all(diff(collision.energy) .> 0) || throw(ArgumentError("collision $k energies must be finite and strictly increasing"))
        all(isfinite, collision.cross_section) && all(collision.cross_section .>= 0) || throw(ArgumentError("collision $k cross sections must be finite and nonnegative"))
        isfinite(collision.threshold) && collision.threshold >= 0 || throw(ArgumentError("collision $k threshold must be finite and nonnegative"))
        if collision.kind == EffectiveCollision || collision.kind == ElasticCollision
            collision.target in elastic_seen && throw(ArgumentError("duplicate effective/elastic collision for target '$(collision.target)'"))
            push!(elastic_seen, collision.target)
        end
    end
    return nothing
end

function _collision_intervals(edges::Vector{Float64}, collision::ElectronCollision)
    n = length(edges) - 1
    shift = collision.kind == IonizationCollision ? 2.0 : 1.0
    low = edges[1] + 1e-9
    high = edges[end] - 1e-9
    high > low || throw(ArgumentError("EEDF energy span is too small for collision cache clipping"))
    shifted = clamp.(shift .* edges .+ collision.threshold, low, high)
    nodes = copy(shifted)
    append!(nodes, (v for v in edges if shifted[1] <= v <= shifted[end]))
    append!(nodes, (v for v in collision.energy if shifted[1] <= v <= shifted[end]))
    sort!(nodes)
    unique!(nodes)
    intervals = _EEDFInterval[]
    for q in 1:length(nodes)-1
        a = nodes[q]
        b = nodes[q+1]
        b > a || continue
        j = searchsortedfirst(edges, b) - 1
        i = searchsortedfirst(shifted, b) - 1
        (1 <= i <= n && 1 <= j <= n) || throw(ArgumentError("collision cache index outside EEDF grid"))
        push!(intervals, _EEDFInterval(i, j, a, b,
              _linear_interp_hold(a, collision.energy, collision.cross_section),
              _linear_interp_hold(b, collision.energy, collision.cross_section)))
    end
    return intervals
end

function EEDFWorkspace(model::EEDFModel, state::EEDFState)
    _validate_model(model)
    all(t -> haskey(state.mole_fractions, t) && haskey(state.molecular_weights, t), model.target_names) ||
        throw(ArgumentError("EEDF state does not cover every model target"))
    edges = model.energy_edges
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    n = length(centers)
    caches = _EEDFCollisionCache[]
    total_center = zeros(n)
    total_edge = zeros(n + 1)
    elastic_edge = zeros(n + 1)
    for collision in model.collisions
        x = state.mole_fractions[collision.target]
        push!(caches, _EEDFCollisionCache(collision, x,
              collision.kind == IonizationCollision ? 2.0 :
              collision.kind == AttachmentCollision ? 0.0 : 1.0,
              _collision_intervals(edges, collision)))
        for j in eachindex(centers)
            total_center[j] += x * _linear_interp_hold(centers[j], collision.energy, collision.cross_section)
        end
        for j in eachindex(edges)
            sigma = _linear_interp_hold(edges[j], collision.energy, collision.cross_section)
            total_edge[j] += x * sigma
            if collision.kind == EffectiveCollision || collision.kind == ElasticCollision
                molecule_mass = state.molecular_weights[collision.target] / _EEDF_AVOGADRO_KMOL
                elastic_edge[j] += 2 * _EEDF_ELECTRON_MASS / molecule_mass * x * sigma
            end
        end
    end
    all(isfinite, total_center) && all(total_center .> 0) || throw(ArgumentError("total center cross section must be finite and positive"))
    all(isfinite, total_edge) && all(total_edge .> 0) || throw(ArgumentError("total edge cross section must be finite and positive"))
    return EEDFWorkspace(centers, caches, total_center, total_edge, elastic_edge,
                         zeros(n), zeros(n, n), zeros(n, n), zeros(n, n))
end

function _integral_pq(a::Float64, b::Float64, u0::Float64, u1::Float64,
                      g::Float64, x0::Float64)
    if g != 0.0
        ag = a*g
        bg = b*g
        ea = expm1(g*(-a + x0))
        eb = expm1(g*(-b + x0))
        A1 = (ea*(ag + 1) + ag - eb*(bg + 1) - bg) / g^2
        A2 = (ea*(2*(ag + 1) + ag^2) + ag*(ag + 2) -
              eb*(2*(bg + 1) + bg^2) - bg*(bg + 2)) / g^3
    else
        A1 = (b^2 - a^2) / 2
        A2 = (b^3 - a^3) / 3
    end
    c0 = (a*u1 - b*u0) / (a - b)
    c1 = (u0 - u1) / (a - b)
    return c0*A1 + c1*A2
end

function _log_slopes!(g::Vector{Float64}, f::Vector{Float64}, centers::Vector{Float64}, floor::Float64)
    n = length(f)
    g[1] = log(max(f[2], floor) / max(f[1], floor)) / (centers[2] - centers[1])
    for i in 2:n-1
        g[i] = log(max(f[i+1], floor) / max(f[i-1], floor)) / (centers[i+1] - centers[i-1])
    end
    g[n] = log(max(f[n], floor) / max(f[n-1], floor)) / (centers[n] - centers[n-1])
    return g
end

function _add_collision_pq!(M::Matrix{Float64}, cache::_EEDFCollisionCache, g, centers, scale=cache.target_fraction)
    collision = cache.collision
    for interval in cache.intervals
        r = _integral_pq(interval.eps_a, interval.eps_b, interval.sigma_a,
                         interval.sigma_b, g[interval.j], centers[interval.j])
        p = _EEDF_GAMMA * r * scale
        M[interval.j, interval.j] -= p
        M[interval.i, interval.j] += cache.incoming_factor * p
    end
    return M
end

function _assemble_pq!(ws::EEDFWorkspace, f::Vector{Float64}, floor::Float64)
    fill!(ws.pq, 0.0)
    _log_slopes!(ws.g, f, ws.centers, floor)
    for cache in ws.collision_cache
        kind = cache.collision.kind
        if kind != EffectiveCollision && kind != ElasticCollision
            _add_collision_pq!(ws.pq, cache, ws.g, ws.centers)
        end
    end
    return ws.pq
end

function _production_frequency(ws::EEDFWorkspace, f::Vector{Float64}, floor::Float64)
    _log_slopes!(ws.g, f, ws.centers, floor)
    nu = 0.0
    scratch = zeros(length(f), length(f))
    for cache in ws.collision_cache
        kind = cache.collision.kind
        if kind == IonizationCollision || kind == AttachmentCollision
            fill!(scratch, 0.0)
            _add_collision_pq!(scratch, cache, ws.g, ws.centers)
            nu += sum(scratch * f)
        end
    end
    isfinite(nu) || throw(ErrorException("non-finite net electron production frequency"))
    return nu
end

function _sg_coefficients(W::Float64, D::Float64, h::Float64)
    z = W*h/D
    isfinite(z) || throw(ErrorException("non-finite Scharfetter-Gummel Peclet number"))
    if abs(z) < 1e-7
        common = D/h
        z2 = z*z
        return common*(1 + z/2 + z2/12), -common*(1 - z/2 + z2/12)
    end
    a0 = W / (-expm1(-z))
    a1 = W / (-expm1(z))
    return a0, a1
end

function _assemble_operator!(ws::EEDFWorkspace, model::EEDFModel, state::EEDFState,
                             f::Vector{Float64}, floor::Float64)
    n = length(f)
    edges = model.energy_edges
    nu = _production_frequency(ws, f, floor)
    a0 = fill(NaN, n + 1)
    a1 = fill(NaN, n + 1)
    density = state.number_density
    omega = 2*pi*state.frequency
    for e in 2:n
        energy = edges[e]
        sigma_tilde = ws.total_edge[e] + nu / (sqrt(energy) * _EEDF_GAMMA)
        isfinite(sigma_tilde) && sigma_tilde != 0 || throw(ErrorException("invalid diffusion cross-section denominator at edge $e"))
        q = omega / (density * _EEDF_GAMMA * sqrt(energy))
        F = sigma_tilde^2 / (sigma_tilde^2 + q^2)
        W = -_EEDF_GAMMA * energy^2 * ws.elastic_edge[e]
        DA = _EEDF_GAMMA / 3 * state.reduced_field^2 * energy
        DB = _EEDF_GAMMA * state.T * _EEDF_BOLTZMANN / _EEDF_ELECTRON_CHARGE * energy^2 * ws.elastic_edge[e]
        D = DA / sigma_tilde * F + DB
        isfinite(D) && D > 0 || throw(ErrorException("nonpositive or non-finite energy diffusion coefficient at edge $e"))
        a0[e], a1[e] = _sg_coefficients(W, D, ws.centers[e] - ws.centers[e-1])
    end

    A = ws.operator
    fill!(A, 0.0)
    A[1,1] = a0[2]
    for j in 2:n-1
        A[j,j] = a0[j+1] - a1[j]
    end
    for j in 1:n-1
        A[j,j+1] = a1[j+1]
        A[j+1,j] = -a0[j+1]
    end
    A[n,n] = -a1[n]
    for i in 1:n
        A[i,i] += 2/3 * (edges[i+1]^1.5 - edges[i]^1.5) * nu
    end
    _assemble_pq!(ws, f, floor)
    A .-= ws.pq
    all(isfinite, A) || throw(ErrorException("non-finite EEDF operator"))
    return nu
end

function _eedf_norm(f::Vector{Float64}, centers::Vector{Float64})
    return _eedf_simpson(f .* sqrt.(centers), centers)
end

function _normalize_eedf!(f::Vector{Float64}, centers::Vector{Float64}; require_nonnegative=false)
    all(isfinite, f) || throw(ErrorException("EEDF contains non-finite values"))
    if require_nonnegative && any(f .< 0)
        bad = findfirst(f .< 0)
        throw(ErrorException("EEDF nonnegativity invariant failed at center index $bad"))
    end
    value = _eedf_norm(f, centers)
    isfinite(value) && value > 0 || throw(ErrorException("EEDF has a nonpositive or non-finite normalization"))
    f ./= value
    return value
end

function _maxwellian(centers::Vector{Float64}, kTe::Float64)
    isfinite(kTe) && kTe > 0 || throw(ArgumentError("Maxwellian energy must be finite and positive"))
    f = @. 2/sqrt(pi) * kTe^(-1.5) * exp(-centers/kTe)
    _normalize_eedf!(f, centers; require_nonnegative=true)
    return f
end

function _electron_mobility(ws::EEDFWorkspace, model::EEDFModel, state::EEDFState,
                            f::Vector{Float64}, floor::Float64)
    nu = _production_frequency(ws, f, floor)
    edges = model.energy_edges
    y = zeros(length(edges))
    for e in 2:length(edges)-1
        df = (f[e] - f[e-1]) / (ws.centers[e] - ws.centers[e-1])
        denom = ws.total_edge[e] + nu / (_EEDF_GAMMA * sqrt(edges[e]))
        isfinite(denom) && denom != 0 || throw(ErrorException("invalid mobility denominator at edge $e"))
        y[e] = edges[e] * df / denom
    end
    density = state.number_density
    mobility = -_EEDF_GAMMA / 3 * _eedf_simpson(y, edges) / density
    isfinite(mobility) || throw(ErrorException("non-finite electron mobility"))
    return mobility
end

function _edge_distribution(edges::Vector{Float64}, centers::Vector{Float64}, f::Vector{Float64})
    return [_linear_interp_hold(x, centers, f) for x in edges]
end

"""
    solve_eedf(model, state; options=TwoTermOptions(), initial=nothing)

Solve the fixed-grid temporal two-term EEDF problem. Molecular weights in
`state` are kg/kmol and `reduced_field` is in V*m^2. Pass a previous
[`EEDFResult`](@ref) as `initial` to continue from its center distribution on the
same grid. The input is copied. At or below `options.low_field_threshold`,
the calculation always starts from the gas-temperature Maxwellian.
"""
function solve_eedf(model::EEDFModel, state::EEDFState;
                    options::TwoTermOptions=TwoTermOptions(),
                    initial::Union{Nothing,EEDFResult}=nothing)
    return _solve_eedf(model,state,options,initial,false)
end

function _solve_accepted_eedf(model::EEDFModel,state::EEDFState;
        options::TwoTermOptions=TwoTermOptions(),initial::Union{Nothing,EEDFResult}=nothing)
    return _solve_eedf(model,state,options,initial,true)
end

function _solve_eedf(model,state,options,initial,accepted_signed)
    _validate_model(model)
    for target in model.target_names
        haskey(state.mole_fractions, target) || throw(ArgumentError("EEDF state is missing target '$target'"))
        haskey(state.molecular_weights, target) || throw(ArgumentError("EEDF state is missing target molecular weight '$target'"))
        xi = state.mole_fractions[target]
        mw = state.molecular_weights[target]
        isfinite(xi) && (accepted_signed || xi >= 0) ||
            throw(ArgumentError("target mole fractions must be finite and nonnegative"))
        isfinite(mw) && mw > 0 || throw(ArgumentError("target molecular weights must be finite and positive kg/kmol"))
    end
    xsum = sum(state.mole_fractions[target] for target in model.target_names)
    isfinite(xsum) && xsum > 0 || throw(ArgumentError("collision-target mole fractions must have a positive sum"))
    abs(xsum - 1) <= 16eps(Float64) || throw(ArgumentError("EEDF state target mole fractions must remain normalized"))
    isfinite(state.T) && state.T > 0 || throw(ArgumentError("EEDF temperature must be finite and positive"))
    isfinite(state.P) && state.P > 0 || throw(ArgumentError("EEDF pressure must be finite and positive"))
    isfinite(state.reduced_field) && state.reduced_field >= 0 || throw(ArgumentError("reduced electric field must be finite and nonnegative"))
    isfinite(state.frequency) && state.frequency >= 0 || throw(ArgumentError("electric field frequency must be finite and nonnegative"))
    isfinite(state.number_density) && state.number_density > 0 || throw(ArgumentError("EEDF number density must be finite and positive m⁻³"))
    options.max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    isfinite(options.delta0) && options.delta0 > 0 || throw(ArgumentError("delta0 must be finite and positive"))
    isfinite(options.factor) && options.factor > 1 || throw(ArgumentError("factor must be finite and greater than one"))
    isfinite(options.rtol) && options.rtol > 0 || throw(ArgumentError("rtol must be finite and positive"))
    isfinite(options.positivity_floor) && options.positivity_floor > 0 || throw(ArgumentError("positivity_floor must be finite and positive"))
    isfinite(options.initial_kTe) && options.initial_kTe > 0 || throw(ArgumentError("initial_kTe must be finite and positive"))
    isfinite(options.low_field_threshold) && options.low_field_threshold >= 0 || throw(ArgumentError("low_field_threshold must be finite and nonnegative"))
    ws = EEDFWorkspace(model, state)

    if initial !== nothing
        initial.edges == model.energy_edges && initial.centers == ws.centers ||
            throw(ArgumentError("initial EEDF must use the same energy grid"))
        length(initial.center_eedf) == length(ws.centers) ||
            throw(ArgumentError("initial center EEDF length does not match the grid"))
        all(isfinite, initial.center_eedf) && all(initial.center_eedf .>= 0) ||
            throw(ArgumentError("initial center EEDF must be finite and nonnegative"))
        norm = _eedf_norm(initial.center_eedf, ws.centers)
        isfinite(norm) && abs(norm - 1) <= 1e-12 ||
            throw(ArgumentError("initial center EEDF must be normalized"))
    end
    low_field = state.reduced_field <= options.low_field_threshold
    kTe = low_field ? _EEDF_BOLTZMANN*state.T/_EEDF_ELECTRON_CHARGE : options.initial_kTe
    f = !low_field && initial !== nothing ? copy(initial.center_eedf) : _maxwellian(ws.centers, kTe)
    errors = Float64[]
    deltas = Float64[]
    iterations = 0
    converged = low_field

    if !converged
        err0 = 0.0
        err1 = 0.0
        delta = options.delta0
        for iteration in 1:options.max_iterations
            if 0 < err1 < err0
                delta *= log(options.factor) / (log(err0) - log(err1))
            end
            isfinite(delta) && delta != 0 || throw(ErrorException("invalid pseudo-time step at iteration $iteration"))
            push!(deltas, delta)
            old = copy(f)
            _assemble_operator!(ws, model, state, old, options.positivity_floor)
            copyto!(ws.system, ws.operator)
            ws.system .*= delta
            for i in axes(ws.system, 1)
                ws.system[i,i] += 1
            end
            f = ws.system \ old
            _normalize_eedf!(f, ws.centers; require_nonnegative=true)
            err0 = err1
            err1 = _eedf_norm(abs.(old .- f), ws.centers)
            isfinite(err1) || throw(ErrorException("non-finite convergence error at iteration $iteration"))
            push!(errors, err1)
            iterations = iteration
            if err1 < options.rtol
                converged = true
                break
            end
        end
        converged || throw(ErrorException("EEDF convergence failed after $(options.max_iterations) iterations; last error=$(errors[end]), delta=$(deltas[end])"))
    end

    edge_f = _edge_distribution(model.energy_edges, ws.centers, f)
    all(isfinite, edge_f) && all(edge_f .>= 0) || throw(ErrorException("edge EEDF nonnegativity invariant failed"))
    mobility = _electron_mobility(ws, model, state, f, options.positivity_floor)
    return EEDFResult(copy(model.energy_edges), edge_f, copy(ws.centers), f,
                      mobility, iterations, converged, errors, deltas)
end
