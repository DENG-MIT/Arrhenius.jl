function _check_real_gas_rate_models(reaction)
    if hasproperty(reaction,:blowers_masel) && !isempty(reaction.blowers_masel.reaction_indices)
        throw(ArgumentError("Redlich–Kwong reactors do not yet support Blowers–Masel rates: they require actual partial-molar enthalpies separately from Kc reference enthalpies"))
    end
    return nothing
end

"Reusable native RK property, physical concentration and activity arrays."
struct RealGasKineticsWorkspace{T}
    thermo::RedlichKwongWorkspace{T}
    kinetics::KineticsWorkspace{T}
    X::Vector{T}
    C::Vector{T}
    activity::Vector{T}
    h0::Vector{T}
    s0::Vector{T}
    wdot::Vector{T}
end
function RealGasKineticsWorkspace(gas::Solution,model::RedlichKwongThermo,::Type{T}=Float64) where {T}
    _check_real_gas_rate_models(gas.reaction)
    gas.species_names == model.species_names && gas.MW ≈ model.MW ||
        throw(ArgumentError("RK thermodynamics and kinetic species must match"))
    return RealGasKineticsWorkspace(RedlichKwongWorkspace(model,T),KineticsWorkspace(gas.reaction,T),
        (zeros(T,gas.n_species) for _ in 1:6)...)
end

"""
    redlich_kwong_rates!(output, gas, model, T, rho, X, work; rate_multipliers=nothing)

Evaluate species production in kmol/m³/s. `X` is normalized mole composition.
Third-body and falloff rates use physical concentrations. Reaction mass action
uses fugacity-based activity concentrations, and PLOG uses the EOS pressure.
Equilibrium constants use ideal reference species h/s at one_atm; the phase
standard-pressure and standard-concentration terms cancel to this reference.
Returns the thermodynamic state, with arrays aliasing `work`.
Supported reaction models are elementary Arrhenius, third-body, Lindemann,
Troe and PLOG. Blowers–Masel reactions are explicitly rejected.
"""
function redlich_kwong_rates!(output,gas::Solution,model::RedlichKwongThermo,T,rho,X,
                             work::RealGasKineticsWorkspace;rate_multipliers=nothing)
    _check_real_gas_rate_models(gas.reaction)
    state = redlich_kwong_properties!(work.thermo,model,T,rho,X)
    RT = R*T
    @inbounds for k in eachindex(X)
        work.C[k] = X[k]/state.v
        work.activity[k] = X[k]*state.P/RT*exp(state.lnphi[k])
        work.h0[k] = RT*work.thermo.h0[k]
        work.s0[k] = R*work.thermo.s0[k]
    end
    wdot!(output,gas.reaction,T,work.C,work.s0,work.h0,work.kinetics;
        pressure=state.P,activity_concentrations=work.activity,rate_multipliers)
    return state
end

"""
    RedlichKwongReactor(gas, model; temperature, pressure=one_atm,
                       mole_fractions=nothing, mass_fractions=nothing,
                       rate_multipliers=nothing)

Closed, constant-volume, adiabatic homogeneous reactor with native RK
thermodynamics and native reaction rates. `gas` supplies kinetic sidecar data;
its species order and molecular weights must match `model`. Supply exactly one
composition. State variables are `[Y₁,…,Yₙ,T]`, in SI units. The energy equation
uses the constant-T,V derivative of total energy with respect to each species'
mole amount, rather than its ordinary partial molar internal energy.
Elementary Arrhenius, third-body, Lindemann, Troe and PLOG reactions are
supported. Blowers–Masel rates require additional nonideal enthalpy coupling
and are rejected by this constructor.

Use the existing caller-owned ODE interface `solve_reactor(...; integrator)`.
No ODE package or Cantera runtime is required by the core implementation.
"""
struct RedlichKwongReactor{G,M,T,V}
    gas::G
    model::M
    temperature::T
    pressure::T
    density::T
    mass_fractions::Vector{T}
    energy::Symbol
    rate_multipliers::V
end
function RedlichKwongReactor(gas::Solution,model::RedlichKwongThermo;temperature,pressure=one_atm,
                            mole_fractions=nothing,mass_fractions=nothing,rate_multipliers=nothing)
    _check_real_gas_rate_models(gas.reaction)
    gas.species_names == model.species_names && gas.MW ≈ model.MW ||
        throw(ArgumentError("RK thermodynamics and kinetic species must match"))
    xor(isnothing(mole_fractions),isnothing(mass_fractions)) ||
        throw(ArgumentError("supply exactly one composition"))
    basis = isnothing(mass_fractions) ? :mole : :mass
    composition = isnothing(mass_fractions) ? mole_fractions : mass_fractions
    state = redlich_kwong_state(model;T=temperature,P=pressure,X=composition,basis)
    Y = state.X.*gas.MW/state.MW
    multipliers = isnothing(rate_multipliers) ? nothing : Float64.(rate_multipliers)
    if !isnothing(multipliers)
        length(multipliers) == gas.n_reactions && all(isfinite,multipliers) && all(>=(0),multipliers) ||
            throw(ArgumentError("one finite nonnegative multiplier per reaction required"))
    end
    return RedlichKwongReactor(gas,model,Float64(temperature),Float64(pressure),state.rho,Y,:adiabatic,multipliers)
end

reactor_state(reactor::RedlichKwongReactor) = vcat(reactor.mass_fractions,reactor.temperature)
function reactor_rhs(reactor::RedlichKwongReactor)
    u0 = reactor_state(reactor)
    return ReactorRHS(reactor,RealGasKineticsWorkspace(reactor.gas,reactor.model,eltype(u0)),
        copy(u0),zero(u0),zero(u0),zero(u0))
end

function (rhs::ReactorRHS{<:RedlichKwongReactor})(du,u,p,t)
    reactor,work = rhs.reactor,rhs.workspace
    gas = reactor.gas
    n = gas.n_species
    length(u) == length(du) == n+1 || throw(DimensionMismatch("reactor state must have n_species+1 entries"))
    total = zero(u[end])
    @inbounds for k in 1:n
        work.X[k] = max(u[k],zero(u[k]))/gas.MW[k]
        total += work.X[k]
    end
    total > 0 || throw(DomainError(total,"positive total composition required"))
    work.X ./= total
    state = redlich_kwong_rates!(work.wdot,gas,reactor.model,u[end],reactor.density,work.X,work;
        rate_multipliers=reactor.rate_multipliers)
    source = zero(u[end])
    @inbounds for k in 1:n
        du[k] = work.wdot[k]*gas.MW[k]/reactor.density
        source += work.wdot[k]*state.u_TV[k]
    end
    du[end] = -source/(reactor.density*state.cv_mass)
    return nothing
end

function reactor_properties(reactor::RedlichKwongReactor,u=reactor_state(reactor))
    n = reactor.gas.n_species
    length(u) == n+1 || throw(DimensionMismatch("reactor state must have n_species+1 entries"))
    Y = copy(view(u,1:n))
    # Diagnostics retain the unprojected mass and elemental inventories. The EOS
    # clips roundoff-negative trial amounts consistently with the ODE callback.
    state = redlich_kwong_state(reactor.model;T=u[end],rho=reactor.density,X=max.(Y,0),basis=:mass)
    return (;temperature=state.T,pressure=state.P,density=state.rho,mole_fractions=state.X,mass_fractions=Y,
        cp=state.cp_mass,cv=state.cv_mass,enthalpy=state.h_mass,internal_energy=state.u_mass,
        entropy=state.s_mass,mass_fraction_sum=sum(Y),elemental_inventory=reactor.gas.ele_matrix*(Y./reactor.gas.MW))
end

function reactor_problem(reactor::RedlichKwongReactor,tspan)
    length(tspan) == 2 && all(isfinite,tspan) && tspan[2] > tspan[1] ||
        throw(ArgumentError("two finite increasing times required"))
    rhs = reactor_rhs(reactor)
    jac = (J,u,p,t) -> reactor_jacobian!(J,u,rhs,t)
    tgrad = (du,u,p,t) -> (fill!(du,zero(eltype(du)));nothing)
    return (f=rhs,jac=jac,tgrad=tgrad,u0=reactor_state(reactor),tspan=Tuple(float.(tspan)),p=nothing)
end
function solve_reactor(reactor::RedlichKwongReactor,tspan;integrator,kwargs...)
    return integrator(reactor_problem(reactor,tspan);kwargs...)
end

export RealGasKineticsWorkspace,redlich_kwong_rates!,RedlichKwongReactor
