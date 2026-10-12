import ForwardDiff

# Native gas thermochemistry composed with the fixed Boltzmann collision model.
struct _PlasmaThermalData
    thermo::IdealGasThermo{Float64}
    reaction::Reaction{Float64}
    chebyshev_indices::Vector{Int}
    chebyshev_coefficients::Vector{Matrix{Float64}}
    chebyshev_temperature_ranges::Vector{Tuple{Float64,Float64}}
    chebyshev_pressure_ranges::Vector{Tuple{Float64,Float64}}
    eedf_model::EEDFModel
    reverse_plan::_ReversibleRatePlan{Reaction{Float64},Vector{Bool}}
end
function _PlasmaThermalData(thermo,reaction,indices,coefficients,tranges,pranges,model)
    _PlasmaThermalData(thermo,reaction,indices,coefficients,tranges,pranges,model,
        _ReversibleRatePlan(reaction))
end

"""
    PlasmaMechanism(path; phase=nothing, data_paths=String[], atomic_weights=Dict())

Read an ideal plasma mechanism directly from Cantera-style YAML.
Imported species files are searched beside the mechanism and in `data_paths`.
Supports isotropic rates and Boltzmann two-term distributions. Boltzmann phases
also support imported gas thermodynamics, reversible Arrhenius reactions, falloff
and Chebyshev rates. Reaction orders are nonnegative integers. Rates use
kmol, m, s, K and J; electron energies are in eV.
"""
struct PlasmaMechanism
    name::String
    n_species::Int
    n_reactions::Int
    species_names::Vector{String}
    element_names::Vector{String}
    MW::Vector{Float64}
    elemental_matrix::Matrix{Float64}
    electron_index::Int
    reactants::Matrix{Int}
    products::Matrix{Int}
    orders::Matrix{Int}
    stoichiometry::Matrix{Float64}
    thirdbody::BitVector
    efficiencies::Matrix{Float64}
    rate_types::Vector{UInt8}
    rate_parameters::Matrix{Float64}
    collision_energy::Vector{Vector{Float64}}
    cross_sections::Vector{Vector{Float64}}
    energy_levels::Vector{Float64}
    shape_factor::Float64
    initial_temperature::Float64
    initial_pressure::Float64
    initial_mean_electron_energy::Float64
    initial_mole_fractions::Vector{Float64}
    thermal::Union{Nothing,_PlasmaThermalData}
end

# Preserve construction of the original isotropic mechanism.
function PlasmaMechanism(name,ns,nr,names,elements,mw,em,electron,react,prod,orders,S,
        third,eff,kinds,params,energies,cross,grid,shape,T,P,E,X)
    PlasmaMechanism(name,ns,nr,names,elements,mw,em,electron,react,prod,orders,S,
        third,eff,kinds,params,energies,cross,grid,shape,T,P,E,X,nothing)
end

_plasma_electron_temperature(E) = (2.0 / 3.0) * E * _EEDF_ELECTRON_CHARGE / _EEDF_BOLTZMANN
_plasma_gamma(x::Float64) = ccall((:tgamma, Base.Math.libm), Float64, (Float64,), x)
function _plasma_positive(value, name)
    x = Float64(value)
    isfinite(x) && x > 0 || throw(ArgumentError("$name must be positive and finite"))
    return x
end

function _plasma_isotropic_eedf!(f, energy, E, shape)
    length(f) == length(energy) || throw(DimensionMismatch("EEDF grid length"))
    E = _plasma_positive(E, "mean electron energy")
    shape = _plasma_positive(shape, "shape factor")
    g1, g2 = _plasma_gamma(1.5 / shape), _plasma_gamma(2.5 / shape)
    c1 = shape * g2^1.5 / g1^2.5
    c2 = (g2 / g1)^shape
    isfinite(c1) && isfinite(c2) && c1 > 0 && c2 > 0 ||
        throw(ArgumentError("isotropic shape factor exceeds floating-point range"))
    inverse_energy = 1.0 / E
    scale = c1 * inverse_energy^1.5
    @inbounds for i in eachindex(energy)
        e = energy[i]
        isfinite(e) && e >= 0 || throw(ArgumentError("nonnegative finite EEDF energies required"))
        f[i] = scale * exp(-c2 * (e * inverse_energy)^shape)
    end
    # This analytic distribution is not renormalized over the truncated grid.
    return f
end

function _plasma_rate_coefficients(m::PlasmaMechanism, T, E)
    T = _plasma_positive(T, "temperature")
    E = _plasma_positive(E, "mean electron energy")
    Te = _plasma_electron_temperature(E)
    f = zeros(length(m.energy_levels))
    _plasma_isotropic_eedf!(f, m.energy_levels, E, m.shape_factor)
    kf = zeros(m.n_reactions)
    integrand = similar(f)
    for j in eachindex(kf)
        A, b, Eg, Ee, bg, inverse_T = m.rate_parameters[:, j]
        kind = m.rate_types[j]
        if kind == 0x01
            kf[j] = A * exp(b * log(T) - Eg / (R * T))
        elseif kind == 0x02
            kf[j] = A * exp(bg * log(T) + b * log(Te) - Eg / (R * T) +
                Ee * (Te - T) / (R * Te * T) - T * inverse_T)
        elseif kind == 0x03
            @inbounds for i in eachindex(f)
                e = m.energy_levels[i]
                sigma = _linear_interp_hold(e, m.collision_energy[j], m.cross_sections[j])
                integrand[i] = e * f[i] * sigma
            end
            kf[j] = _EEDF_GAMMA * _EEDF_AVOGADRO_KMOL * _eedf_simpson(integrand, m.energy_levels)
        else
            throw(ArgumentError("unsupported plasma rate type: $kind"))
        end
    end
    all(isfinite, kf) || throw(DomainError(kf, "nonfinite plasma rate coefficient"))
    return kf, f
end

"""
    PlasmaState(mechanism; temperature, pressure, mole_fractions, mass_fractions,
                mean_electron_energy)

An ideal plasma state with a separate electron temperature. Defaults come from
its YAML phase; specify at most one composition. Input fractions are normalized.
Temperatures are K, pressure Pa, density kg/m³ and mean electron energy eV.
Use [`set_mean_electron_energy!`](@ref) to change electron energy at fixed density.
"""
mutable struct PlasmaState
    mechanism::PlasmaMechanism
    temperature::Float64
    density::Float64
    mean_electron_energy::Float64
    mass_fractions::Vector{Float64}
    electric_field::Float64
    eedf::Union{Nothing,EEDFResult}
end
PlasmaState(m,T,rho,E,Y) = PlasmaState(m,T,rho,E,Y,0.0,nothing)

function _plasma_mass_fractions(m, X, Y)
    if X === nothing
        return _reactor_composition(m, Y, Float64)
    end
    Y === nothing || throw(ArgumentError("supply only one composition"))
    x = _reactor_composition(m, X, Float64)
    return x .* m.MW ./ dot(x, m.MW)
end
function _plasma_density(m, Y, T, Te, P)
    inverse_mw = sum(Y ./ m.MW)
    electron_fraction = Y[m.electron_index] / (m.MW[m.electron_index] * inverse_mw)
    mean_T = T + electron_fraction * (Te - T)
    return P / (R * inverse_mw * mean_T)
end
function PlasmaState(m::PlasmaMechanism; temperature=m.initial_temperature,
                     pressure=m.initial_pressure, mole_fractions=nothing,
                     mass_fractions=nothing,
                     mean_electron_energy=m.initial_mean_electron_energy)
    T = _plasma_positive(temperature, "temperature")
    P = _plasma_positive(pressure, "pressure")
    E = _plasma_positive(mean_electron_energy, "mean electron energy")
    if mole_fractions === nothing && mass_fractions === nothing
        mole_fractions = m.initial_mole_fractions
    end
    Y = _plasma_mass_fractions(m, mole_fractions, mass_fractions)
    rho = _plasma_density(m, Y, T, _plasma_electron_temperature(E), P)
    return PlasmaState(m, T, rho, E, Y)
end

"Return independent composition arrays and the two-temperature ideal-plasma properties."
function plasma_properties(s::PlasmaState)
    m = s.mechanism
    Y = copy(s.mass_fractions)
    inverse_mw = sum(Y ./ m.MW)
    X = Y ./ m.MW ./ inverse_mw
    Te = _plasma_electron_temperature(s.mean_electron_energy)
    P = R * s.density * inverse_mw * (s.temperature + X[m.electron_index] * (Te - s.temperature))
    basic = (T=s.temperature, Te, P, rho=s.density, X, Y,
        mean_electron_energy=s.mean_electron_energy)
    m.thermal === nothing && return basic
    N = s.density * inverse_mw * _EEDF_AVOGADRO_KMOL
    mobility = s.eedf === nothing ? nothing : s.eedf.mobility
    joule = mobility === nothing ? nothing :
        _EEDF_ELECTRON_CHARGE * N * X[m.electron_index] * mobility * s.electric_field^2
    return merge(basic, (electric_field=s.electric_field,
        reduced_electric_field=s.electric_field/N, electron_mobility=mobility,
        joule_heating=joule))
end

"""
    set_plasma_state!(state; temperature, pressure, mole_fractions, mass_fractions)

Set gas temperature, pressure and optionally composition while retaining electron
energy. The density follows the two-temperature equation of state. With no
composition argument the current mass fractions are retained.
"""
function set_plasma_state!(s::PlasmaState; temperature=s.temperature,
                           pressure=plasma_properties(s).P,
                           mole_fractions=nothing, mass_fractions=nothing)
    T = _plasma_positive(temperature, "temperature")
    P = _plasma_positive(pressure, "pressure")
    Y = mole_fractions === nothing && mass_fractions === nothing ?
        copy(s.mass_fractions) : _plasma_mass_fractions(s.mechanism, mole_fractions, mass_fractions)
    rho = _plasma_density(s.mechanism, Y, T, _plasma_electron_temperature(s.mean_electron_energy), P)
    s.temperature = T
    s.density = rho
    copyto!(s.mass_fractions, Y)
    return s
end

"Change electron energy in eV, retaining gas temperature, composition and density."
function set_mean_electron_energy!(s::PlasmaState, energy)
    s.mean_electron_energy = _plasma_positive(energy, "mean electron energy")
    return s
end

# Integer powers preserve the same polynomial at zero and signed solver trial states.
function _plasma_mass_action!(q, C, m, kf)
    @inbounds for j in eachindex(q)
        value = kf[j] * one(eltype(C))
        for i in eachindex(C)
            order = m.orders[i, j]
            if order != 0
                power = one(eltype(C))
                for _ in 1:order
                    power *= C[i]
                end
                value *= power
            end
        end
        if m.thirdbody[j]
            collider = zero(eltype(C))
            for i in eachindex(C)
                collider += m.efficiencies[i, j] * C[i]
            end
            value *= collider
        end
        q[j] = value
    end
    return nothing
end

"Return native rate coefficients, progress/production rates and isotropic EEDF at a state."
function plasma_rates(s::PlasmaState)
    m = s.mechanism
    m.thermal === nothing || return _plasma_thermal_rates(s)
    kf, f = _plasma_rate_coefficients(m, s.temperature, s.mean_electron_energy)
    C = s.mass_fractions .* s.density ./ m.MW
    q = similar(kf)
    _plasma_mass_action!(q, C, m, kf)
    wdot = m.stoichiometry * q
    return (forward_rate_constants=kf, net_rates_of_progress=q,
        net_production_rates=wdot, concentrations=C,
        dYdt=wdot .* m.MW ./ s.density,
        electron_energy_distribution=f, electron_energy_levels=copy(m.energy_levels))
end

"""
    PlasmaReactor(state)

A homogeneous constant-pressure reactor with fixed gas and electron temperatures.
Snapshots the supplied state; later state setters do not modify this reactor.
The ODE state contains species mass fractions only. Chemical energy evolution,
electric-field coupling and evolving EEDF are not part of this isothermal model.
"""
struct PlasmaReactor
    mechanism::PlasmaMechanism
    temperature::Float64
    electron_temperature::Float64
    pressure::Float64
    mean_electron_energy::Float64
    mass_fractions::Vector{Float64}
end
function PlasmaReactor(s::PlasmaState)
    s.mechanism.thermal === nothing || throw(ArgumentError("Boltzmann plasma requires an energy-coupled reactor; PlasmaReactor fixes both temperatures"))
    properties = plasma_properties(s)
    return PlasmaReactor(s.mechanism, s.temperature, properties.Te, properties.P,
        s.mean_electron_energy, copy(s.mass_fractions))
end
reactor_state(r::PlasmaReactor) = copy(r.mass_fractions)

mutable struct PlasmaRHS
    reactor::PlasmaReactor
    rate_constants::Vector{Float64}
    density_weights::Vector{Float64}
    C::Vector{Float64}
    q::Vector{Float64}
    calls::Int
end
function reactor_rhs(r::PlasmaReactor)
    m = r.mechanism
    m.thermal === nothing || throw(ArgumentError("Boltzmann plasma requires an energy-coupled reactor"))
    kf, _ = _plasma_rate_coefficients(m, r.temperature, r.mean_electron_energy)
    a = r.temperature ./ m.MW
    a[m.electron_index] = r.electron_temperature / m.MW[m.electron_index]
    return PlasmaRHS(r, kf, a, zeros(m.n_species), zeros(m.n_reactions), 0)
end
function (rhs::PlasmaRHS)(du, u)
    r, m = rhs.reactor, rhs.reactor.mechanism
    length(u) == length(du) == m.n_species || throw(DimensionMismatch("plasma species state length"))
    rhs.calls += 1
    C = eltype(u) === Float64 ? rhs.C : similar(u)
    q = eltype(u) === Float64 ? rhs.q : similar(u, m.n_reactions)
    denominator = zero(eltype(u))
    @inbounds for i in eachindex(u)
        denominator += rhs.density_weights[i] * u[i]
    end
    isfinite(denominator) && denominator > 0 ||
        throw(DomainError(denominator, "positive finite plasma density denominator required"))
    rho = r.pressure / (R * denominator)
    @inbounds for i in eachindex(u)
        C[i] = rho * u[i] / m.MW[i]
    end
    _plasma_mass_action!(q, C, m, rhs.rate_constants)
    @inbounds for i in eachindex(u)
        total = zero(eltype(u))
        for j in eachindex(q)
            total += m.stoichiometry[i, j] * q[j]
        end
        du[i] = m.MW[i] / rho * total
    end
    return nothing
end
(rhs::PlasmaRHS)(du, u, p, t) = rhs(du, u)

struct _PlasmaJacobian{F,C}
    rhs::F
    config::C
    output::Vector{Float64}
    calls::Base.RefValue{Int}
end
function _PlasmaJacobian(rhs::PlasmaRHS, u)
    output = zeros(length(u))
    config = ForwardDiff.JacobianConfig(rhs, output, u)
    return _PlasmaJacobian(rhs, config, output, Ref(0))
end
function (jac::_PlasmaJacobian)(J, u, p, t)
    jac.calls[] += 1
    ForwardDiff.jacobian!(J, jac.rhs, jac.output, u, jac.config)
    return nothing
end
function reactor_jacobian!(J, u, rhs::PlasmaRHS, t=0)
    _PlasmaJacobian(rhs, u)(J, u, nothing, t)
    return nothing
end
function reactor_problem(r::PlasmaReactor, tspan)
    length(tspan) == 2 && all(isfinite, tspan) && tspan[2] > tspan[1] ||
        throw(ArgumentError("tspan must be finite and strictly increasing"))
    u0 = reactor_state(r)
    rhs = reactor_rhs(r)
    jac = _PlasmaJacobian(rhs, u0)
    tgrad = (du, u, p, t) -> (fill!(du, zero(eltype(du))); nothing)
    return (f=rhs, jac, tgrad, u0, tspan=(float(tspan[1]), float(tspan[2])), p=nothing)
end
function solve_reactor(r::PlasmaReactor, tspan; integrator, kwargs...)
    return integrator(reactor_problem(r, tspan); kwargs...)
end
function reactor_properties(r::PlasmaReactor, u=reactor_state(r))
    m = r.mechanism
    length(u) == m.n_species || throw(DimensionMismatch("plasma species state length"))
    all(isfinite, u) || throw(DomainError(u, "nonfinite plasma state"))
    Y = copy(u)
    inverse_mw = sum(Y ./ m.MW)
    X = Y ./ m.MW ./ inverse_mw
    rho = _plasma_density(m, Y, r.temperature, r.electron_temperature, r.pressure)
    return (T=r.temperature, Te=r.electron_temperature, P=r.pressure, rho, X, Y,
        mass_fraction_sum=sum(Y), elemental_inventory=m.elemental_matrix * (Y ./ m.MW),
        mean_electron_energy=r.mean_electron_energy)
end

export PlasmaMechanism, PlasmaState, PlasmaReactor, PlasmaRHS
export plasma_properties, plasma_rates, set_plasma_state!, set_mean_electron_energy!
