"""
    PlasmaEnergyReactor(state; volume=1.0)

A closed, constant-pressure Boltzmann-plasma reactor. Its state is
`[mass, total_enthalpy, mass_fractions...]`; electron temperature and the
electric field/EEDF are held fixed between explicit [`update_eedf!`](@ref)
calls on the prepared RHS.
"""
struct PlasmaEnergyReactor
    mechanism::PlasmaMechanism
    pressure::Float64
    electron_temperature::Float64
    mean_electron_energy::Float64
    electric_field::Float64
    eedf::EEDFResult
    temperature::Float64
    mass::Float64
    total_enthalpy::Float64
    mass_fractions::Vector{Float64}
end

function PlasmaEnergyReactor(s::PlasmaState;volume=1.0)
    s.mechanism.thermal===nothing &&
        throw(ArgumentError("PlasmaEnergyReactor requires a Boltzmann plasma mechanism"))
    s.eedf===nothing && throw(ArgumentError("call update_eedf! before constructing the reactor"))
    V=_plasma_positive(volume,"volume")
    p=plasma_properties(s)
    thermo=plasma_thermodynamics(s)
    mass=s.density*V
    return PlasmaEnergyReactor(s.mechanism,p.P,p.Te,s.mean_electron_energy,
        s.electric_field,deepcopy(s.eedf),s.temperature,mass,mass*thermo.h_mass,
        copy(s.mass_fractions))
end

reactor_state(r::PlasmaEnergyReactor)=vcat(r.mass,r.total_enthalpy,r.mass_fractions)

# This private tag owns a Float64-only prepared first derivative. The energy
# reactor does not promise support for differentiation through its Jacobian.
struct _PlasmaEnergyJacobianTag end
mutable struct _PlasmaEnergySpeciesKernel{W}
    mechanism::PlasmaMechanism
    pressure::Float64
    electron_temperature::Float64
    workspace::W
end

function (kernel::_PlasmaEnergySpeciesKernel)(out,z)
    m=kernel.mechanism
    T=z[1]
    Y=@view z[2:end]
    denominator=zero(T)
    @inbounds for k in eachindex(Y)
        amount=Y[k]/m.MW[k]
        denominator+=(k==m.electron_index ? kernel.electron_temperature : T)*amount
    end
    rho=kernel.pressure/(R*denominator)
    _plasma_thermal_sources!(kernel.workspace,m,T,kernel.electron_temperature,
        kernel.pressure,rho,Y)
    @inbounds for k in eachindex(Y)
        out[k]=m.MW[k]*kernel.workspace.wdot[k]/rho
    end
    return nothing
end

struct _PlasmaEnergyJacobian{K,C}
    kernel::K
    config::C
    output::Vector{Float64}
    z::Vector{Float64}
    species_jacobian::Matrix{Float64}
end

function _PlasmaEnergyJacobian(r::PlasmaEnergyReactor)
    return _PlasmaEnergyJacobian(r,Val(min(8,r.mechanism.n_species+1)))
end

function _PlasmaEnergyJacobian(r::PlasmaEnergyReactor,::Val{N}) where {N}
    ns=r.mechanism.n_species
    DualT=ForwardDiff.Dual{_PlasmaEnergyJacobianTag,Float64,N}
    workspace=_PlasmaThermoRateWorkspace(r.mechanism,DualT)
    kernel=_PlasmaEnergySpeciesKernel(r.mechanism,r.pressure,
        r.electron_temperature,workspace)
    output=zeros(ns)
    z=zeros(ns+1)
    config=ForwardDiff.JacobianConfig(kernel,output,z,ForwardDiff.Chunk{N}(),
        _PlasmaEnergyJacobianTag())
    return _PlasmaEnergyJacobian(kernel,config,output,z,zeros(ns,ns+1))
end

mutable struct PlasmaEnergyRHS{W,J}
    reactor::PlasmaEnergyReactor
    workspace::W
    eedf::EEDFResult
    electric_field::Float64
    mobility::Float64
    last_temperature::Float64
    jacobian::J
end

function reactor_rhs(r::PlasmaEnergyReactor)
    w=_PlasmaThermoRateWorkspace(r.mechanism,Float64)
    _plasma_collision_rates!(w.collision_rates,w.collision_integrand,r.mechanism,r.eedf)
    return PlasmaEnergyRHS(r,w,deepcopy(r.eedf),r.electric_field,r.eedf.mobility,
        r.temperature,_PlasmaEnergyJacobian(r))
end

function _plasma_hp_temperature!(w,m,target,Y,Te,guess;rtol=1e-13,maxiter=500)
    isfinite(target) || throw(DomainError(target,"finite specific enthalpy required"))
    T=guess
    e=m.electron_index
    for _ in 1:maxiter
        isfinite(T) && T>0 || throw(DomainError(T,"positive finite gas temperature required"))
        _plasma_species_thermo!(w.h,w.cp,w.s0,m,T,Te)
        value=zero(T)
        scale=one(T)
        slope=zero(T)
        @inbounds for k in eachindex(Y)
            amount=Y[k]/m.MW[k]
            term=amount*w.h[k]
            value+=term
            scale+=abs(term)
            k==e || (slope+=amount*w.cp[k])
        end
        residual=target-value
        abs(residual)<=rtol*max(abs(target),scale) && return T
        isfinite(slope) && slope>0 ||
            throw(DomainError(slope,"positive heavy-species heat capacity required during plasma HP inversion"))
        step=clamp(residual/slope,-100.0,100.0)
        T+=max(step,-0.5T)
    end
    throw(ErrorException("plasma HP inversion did not converge in $maxiter iterations"))
end

function _plasma_energy_state!(rhs::PlasmaEnergyRHS,u)
    r,m,w=rhs.reactor,rhs.reactor.mechanism,rhs.workspace
    length(u)==m.n_species+2 || throw(DimensionMismatch("energy-plasma state must contain mass, total enthalpy, and all species"))
    mass=u[1]
    isfinite(mass) && mass>0 || throw(DomainError(mass,"positive finite reactor mass required"))
    H=u[2]
    isfinite(H) || throw(DomainError(H,"finite total enthalpy required"))
    Y=@view u[3:end]
    all(isfinite,Y) || throw(DomainError(Y,"finite signed mass fractions required"))
    T=_plasma_hp_temperature!(w,m,H/mass,Y,r.electron_temperature,rhs.last_temperature)
    rhs.last_temperature=T
    inverse_mw=zero(T)
    denominator=zero(T)
    @inbounds for k in eachindex(Y)
        amount=Y[k]/m.MW[k]
        inverse_mw+=amount
        denominator+=(k==m.electron_index ? r.electron_temperature : T)*amount
    end
    isfinite(inverse_mw) && inverse_mw>0 || throw(DomainError(inverse_mw,"positive finite inverse molecular weight required"))
    isfinite(denominator) && denominator>0 || throw(DomainError(denominator,"positive finite two-temperature EOS denominator required"))
    rho=r.pressure/(R*denominator)
    @inbounds for k in eachindex(Y)
        w.X[k]=(Y[k]/m.MW[k])/inverse_mw
    end
    return T,rho,inverse_mw,Y
end

function (rhs::PlasmaEnergyRHS)(du,u,p,t)
    m=rhs.reactor.mechanism
    length(du)==m.n_species+2 || throw(DimensionMismatch("derivative must match energy-plasma state"))
    T,rho,_,Y=_plasma_energy_state!(rhs,u)
    _plasma_thermal_sources!(rhs.workspace,m,T,rhs.reactor.electron_temperature,
        rhs.reactor.pressure,rho,Y)
    du[1]=0.0
    ye=Y[m.electron_index]
    if ye>0 && rhs.mobility>0 && rhs.electric_field>0
        factor=_EEDF_ELECTRON_CHARGE*_EEDF_AVOGADRO_KMOL/m.MW[m.electron_index]
        du[2]=u[1]*ye*factor*rhs.mobility*rhs.electric_field^2
    else
        du[2]=0.0
    end
    @inbounds for k in eachindex(Y)
        du[k+2]=m.MW[k]*rhs.workspace.wdot[k]/rho
    end
    return nothing
end
(rhs::PlasmaEnergyRHS)(du,u)=rhs(du,u,nothing,0.0)

function reactor_properties(rhs::PlasmaEnergyRHS,u)
    r,m,w=rhs.reactor,rhs.reactor.mechanism,rhs.workspace
    T,rho,inverse_mw,Y=_plasma_energy_state!(rhs,u)
    h_mass=zero(T)
    @inbounds for k in eachindex(Y)
        h_mass+=(Y[k]/m.MW[k])*w.h[k]
    end
    return (T=T,Te=r.electron_temperature,P=r.pressure,rho=rho,X=copy(w.X),Y=copy(Y),
        mass=u[1],h_mass=h_mass,volume=u[1]/rho,electric_field=rhs.electric_field,
        reduced_electric_field=rhs.electric_field/(rho*inverse_mw*_EEDF_AVOGADRO_KMOL),
        electron_mobility=rhs.mobility,mean_electron_energy=r.mean_electron_energy,
        mass_fraction_sum=sum(Y),elemental_inventory=m.elemental_matrix*(Y./m.MW))
end
reactor_properties(r::PlasmaEnergyReactor,u=reactor_state(r))=reactor_properties(reactor_rhs(r),u)

"""Refresh the accepted-state EEDF and cached collision rate coefficients."""
function update_eedf!(rhs::PlasmaEnergyRHS,u;reduced_field=nothing,options=TwoTermOptions())
    p=reactor_properties(rhs,u)
    m=rhs.reactor.mechanism
    N=p.rho*sum(p.Y./m.MW)*_EEDF_AVOGADRO_KMOL
    isfinite(N) && N>0 || throw(DomainError(N,"positive finite number density required"))
    changing_field=reduced_field!==nothing
    EN=changing_field ? Float64(reduced_field) : rhs.electric_field/N
    isfinite(EN) && EN>=0 || throw(ArgumentError("E/N must be finite and nonnegative"))
    state=_accepted_eedf_state(m.thermal.eedf_model;T=p.T,P=p.P,
        mole_fractions=Dict(zip(m.species_names,p.X)),
        molecular_weights=Dict(zip(m.species_names,m.MW)),reduced_field=EN,number_density=N)
    result=_solve_accepted_eedf(m.thermal.eedf_model,state;options=options,initial=rhs.eedf)
    result.converged || throw(ErrorException("Boltzmann EEDF did not converge"))
    rhs.eedf=result
    rhs.mobility=result.mobility
    changing_field && (rhs.electric_field=EN*N)
    _plasma_collision_rates!(rhs.workspace.collision_rates,
        rhs.workspace.collision_integrand,m,result)
    return rhs
end

function reactor_jacobian!(J,u,rhs::PlasmaEnergyRHS,t=0)
    n=length(u)
    m=rhs.reactor.mechanism
    size(J)==(n,n) || throw(DimensionMismatch("Jacobian must have state dimensions"))
    eltype(u)===Float64 ||
        throw(ArgumentError("prepared energy-plasma Jacobian requires a Float64 state"))

    # Recover the scalar HP root once. This also leaves root-state h and cp in
    # the Float64 workspace for the exact implicit-temperature chain below.
    T,_,_,Y=_plasma_energy_state!(rhs,u)
    hp=rhs.workspace
    capacity=0.0
    @inbounds for k in eachindex(Y)
        k==m.electron_index || (capacity+=Y[k]*hp.cp[k]/m.MW[k])
    end
    isfinite(capacity) && capacity>0 ||
        throw(DomainError(capacity,"positive heavy-species heat capacity required for plasma Jacobian"))

    jac=rhs.jacobian
    jac.z[1]=T
    copyto!(jac.z,2,Y,1,length(Y))
    @inbounds for i in eachindex(hp.collision_rates)
        jac.kernel.workspace.collision_rates[i]=hp.collision_rates[i]
    end
    ForwardDiff.jacobian!(jac.species_jacobian,jac.kernel,jac.output,jac.z,jac.config)

    mass=u[1]
    H=u[2]
    dTdm=-H/(mass*mass*capacity)
    dTdH=1/(mass*capacity)
    fill!(J,0.0)
    @inbounds for i in 1:m.n_species
        row=i+2
        dsource_dT=jac.species_jacobian[i,1]
        J[row,1]=dsource_dT*dTdm
        J[row,2]=dsource_dT*dTdH
        for k in 1:m.n_species
            dTdY=-hp.h[k]/(m.MW[k]*capacity)
            J[row,k+2]=jac.species_jacobian[i,k+1]+dsource_dT*dTdY
        end
    end

    e=m.electron_index
    if Y[e]>0 && rhs.mobility>0 && rhs.electric_field>0
        factor=_EEDF_ELECTRON_CHARGE*_EEDF_AVOGADRO_KMOL/m.MW[e]*
            rhs.mobility*rhs.electric_field^2
        J[2,1]=Y[e]*factor
        J[2,e+2]=mass*factor
    end
    return nothing
end

function reactor_problem(r::PlasmaEnergyReactor,tspan)
    length(tspan)==2 && all(isfinite,tspan) && tspan[2]>tspan[1] ||
        throw(ArgumentError("tspan must be finite and strictly increasing"))
    rhs=reactor_rhs(r)
    jac=(J,u,p,t)->reactor_jacobian!(J,u,rhs,t)
    tgrad=(du,u,p,t)->(fill!(du,zero(eltype(du)));nothing)
    return (f=rhs,jac=jac,tgrad=tgrad,u0=reactor_state(r),
        tspan=(float(tspan[1]),float(tspan[2])),p=nothing)
end

solve_reactor(r::PlasmaEnergyReactor,tspan;integrator,kwargs...)=
    integrator(reactor_problem(r,tspan);kwargs...)

export PlasmaEnergyReactor,PlasmaEnergyRHS
