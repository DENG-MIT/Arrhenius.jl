# Reacting wall boundary equations follow ReactingSurf1D::eval at Cantera
# 726522be4e2a13454d8415b7ef799d621f665cf3. Gas and surface kinetics are native.
"An axisymmetric jet impinging on an isothermal, chemically reacting surface."
mutable struct CatalyticImpingingJet{F} <: AbstractPremixedCounterflow
    flow::F
    surface::SurfaceMechanism
    coverages::Vector{Float64}
    fixed_coverages::Vector{Float64}
    gas_multiplier::Float64
    surface_multiplier::Float64
    coverage_enabled::Bool
    inlet_start::Union{Nothing,Vector{Float64}}
end
@inline Base.getproperty(f::CatalyticImpingingJet,k::Symbol) =
    hasfield(typeof(f),k) ? getfield(f,k) : getproperty(getfield(f,:flow),k)
@inline Base.setproperty!(f::CatalyticImpingingJet,k::Symbol,v) =
    hasfield(typeof(f),k) ? setfield!(f,k,v) : setproperty!(getfield(f,:flow),k,v)
Base.propertynames(f::CatalyticImpingingJet,private::Bool=false) =
    (propertynames(f.flow,private)...,fieldnames(typeof(f))...)

"""
    CatalyticImpingingJet(gas, surface; reactants, mdot, T_inlet=300,
                         T_surface=900, P=one_atm, width=.1, coverages)

Couple an ideal gas jet to stationary surface coverages. The wall has fixed
temperature and zero normal velocity. Surface gas production sets the wall
species flux; surface reaction heat is removed by the prescribed thermostat.
The gas must be the surface mechanism's companion phase, in its original order.
Persistent bulk phases, deposition, and moving walls are not supported.
Initialization uses a linear temperature profile and inlet gas composition.
"""
function CatalyticImpingingJet(gas,surface::SurfaceMechanism;reactants,mdot,
        T_inlet=300.,T_surface=900.,P=one_atm,width=.1,grid=nothing,
        coverages=surface.initial_coverages)
    isempty(surface.bulk_molar_volumes) || throw(ArgumentError("catalytic wall requires no persistent bulk mass exchange"))
    gas.species_names==surface.species_names[surface.n_surface+1:surface.n_surface+surface.n_gas] ||
        throw(ArgumentError("surface companion gas species order required"))
    s=IdealSurface(surface;temperature=T_surface,pressure=P,
        mole_fractions=mole_fractions(gas,reactants),coverages)
    flow=ImpingingJet(gas;reactants,mdot,T_inlet,T_surface,P,width,grid,initial_products=:inlet)
    return CatalyticImpingingJet(flow,surface,copy(s.coverages),copy(s.coverages),1.,1.,true,nothing)
end

"Change reaction multipliers and coverage equation mode, retaining the flow profile."
function set_catalytic_reactions!(f::CatalyticImpingingJet;
        gas_multiplier=f.gas_multiplier,surface_multiplier=f.surface_multiplier,
        coverage_enabled=f.coverage_enabled)
    all(x->isfinite(x)&&x>=0,(gas_multiplier,surface_multiplier)) ||
        throw(ArgumentError("finite nonnegative reaction multipliers required"))
    f.gas_multiplier=Float64(gas_multiplier); f.surface_multiplier=Float64(surface_multiplier)
    f.coverage_enabled=Bool(coverage_enabled)
    !coverage_enabled && (f.fixed_coverages .= f.coverages)
    f.converged=false
    return f
end

"Change the inlet composition while retaining the current catalytic solution."
function set_catalytic_inlet!(f::CatalyticImpingingJet,composition;
        temperature=f.fuel_temperature,mdot=f.fuel_mass_flux)
    isfinite(temperature)&&200<=temperature<=6000 && isfinite(mdot)&&mdot>0 ||
        throw(ArgumentError("inlet temperature must be in 200–6000 K and mass flux positive"))
    X=mole_fractions(f.gas,composition)
    if f.converged && isnothing(f.inlet_start)
        f.inlet_start=mole_fractions(f.gas,f.fuel_Y;basis=:mass)
    end
    f.fuel_Y=X.*f.gas.MW./dot(X,f.gas.MW)
    f.fuel_temperature=Float64(temperature); f.fuel_mass_flux=Float64(mdot)
    f.converged=false
    return f
end

"""
Initialize stationary coverages at the inlet composition and wall temperature.
`integrator` is a caller-owned stiff native Julia surface ODE integrator, as in
`solve_surface`. This initialization never uses a reference flow solution.
"""
function initialize_catalytic_coverages!(f::CatalyticImpingingJet;integrator,
        reltol=1e-10,abstol=1e-16)
    s=IdealSurface(f.surface;temperature=f.oxidizer_temperature,pressure=f.pressure,
        mole_fractions=mole_fractions(f.gas,f.fuel_Y;basis=:mass),coverages=f.coverages)
    solution=solve_surface(s,(0.,1.);integrator,reltol,abstol)
    theta=max.(solution.u[end],0.); theta ./= sum(theta)
    rhs=surface_rhs(s); r=similar(theta); J=zeros(length(theta),length(theta))
    for iteration in 1:30
        rhs(r,theta,nothing,0.)
        norm(r,Inf)<1e-6 && abs(sum(theta)-1)<1e-12 && break
        surface_jacobian!(J,theta,rhs)
        J[1,:].=1.; r[1]=sum(theta)-1
        step=-(J\r)
        alpha=1.
        for k in eachindex(theta)
            step[k]<0 && (alpha=min(alpha,-.99theta[k]/step[k]))
        end
        theta .+= alpha.*step
    end
    rhs(r,theta,nothing,0.)
    norm(r,Inf)<1e-5 || error("initial catalytic coverages did not become stationary")
    f.coverages .= theta; f.fixed_coverages .= theta
    return f
end

struct CatalyticFlameWorkspace{W}
    flow::W
    surface::SurfaceWorkspace
    coverage_residual::Vector{Float64}
end
@inline Base.getproperty(w::CatalyticFlameWorkspace,k::Symbol) =
    hasfield(typeof(w),k) ? getfield(w,k) : getproperty(getfield(w,:flow),k)
CounterflowWorkspace(f::CatalyticImpingingJet)=CatalyticFlameWorkspace(
    CounterflowWorkspace(f.flow),SurfaceWorkspace(f.surface),zeros(f.surface.n_surface))

# Coverage residuals are scaled in seconds, independently of the gas equations.
const _catalytic_coverage_timescale=1e-4
function _catalytic_surface_residual!(r,f,u,w;theta=f.coverages,previous=nothing,dt=Inf)
    n,N=f.gas.n_species,length(f.grid)
    X=[max(u[k+1,N],0.)/f.gas.MW[k] for k in 1:n]; X ./= sum(X)
    surface_rates!(w.surface,f.surface,f.oxidizer_temperature,f.pressure,X,theta;
        gas_temperature=1000*u[1,N])
    if f.coverage_enabled
        r .= _catalytic_coverage_timescale*f.surface_multiplier .* w.surface.coverage_rates
        if previous!==nothing
            r .-= _catalytic_coverage_timescale/dt .* (theta.-previous.coverages)
        end
        r[1]=sum(theta)-1
    else
        r .= theta.-f.fixed_coverages
    end
    return r
end

function counterflow_residual!(r,f::CatalyticImpingingJet,u=f.state,
        w=CounterflowWorkspace(f);previous=nothing,dt=Inf,kwargs...)
    counterflow_residual!(r,f.flow.flow,u,w.flow;
        previous=isnothing(previous) ? nothing : previous.gas,dt,kwargs...)
    p=w.properties; n,N=f.gas.n_species,length(f.grid)
    # Extend the final face's diffusion law smoothly through negative Newton
    # trial traces. At a physical composition this is the unchanged mixture
    # mole-gradient flux. Clipping a wall trace inside this boundary equation
    # would make its derivative vanish and leave an unconstrained wall species.
    left=N-1; dz=f.grid[N]-f.grid[left]
    inverse_left=sum(u[k+1,left]/f.gas.MW[k] for k in 1:n)
    inverse_right=sum(u[k+1,N]/f.gas.MW[k] for k in 1:n)
    signed_flux=[-p.diffusion_prefactor[k,left]*(u[k+1,N]/f.gas.MW[k]/inverse_right-
        u[k+1,left]/f.gas.MW[k]/inverse_left)/dz for k in 1:n]
    fluxsum=sum(signed_flux)
    for k in 1:n
        signed_flux[k]-=u[k+1,left]*fluxsum
    end
    delta_flux=signed_flux.-p.flux[:,left]
    p.flux[:,left].=signed_flux
    cell=.5*(f.grid[N]-f.grid[N-2])
    jm,jp=u[end,left]>0 ? (left-1,left) : (left,N)
    enthalpy_change=0.; cpwall=0.
    for k in 1:n
        if k!=f.dependent_species
            r[k+1,left]-=_flame_timescale/p.rho[left]*delta_flux[k]/cell
            r[k+1,N]+=delta_flux[k]
        end
        enthalpy_change+=.5*delta_flux[k]*(p.h[k,jp]-p.h[k,jm])/((f.grid[jp]-f.grid[jm])*f.gas.MW[k])
        cpwall+=u[k+1,left]*p.cp[k,left]/f.gas.MW[k]
    end
    isempty(f.fixed_temperature) && (r[1,left]-=_flame_timescale/(1000*p.rho[left]*cpwall)*enthalpy_change)
    if f.gas_multiplier!=1
        for j in 2:N-1
            chemical=0.; cp=0.
            for k in 1:n
                k!=f.dependent_species && (r[k+1,j]+=(f.gas_multiplier-1)*
                    _flame_timescale/p.rho[j]*f.gas.MW[k]*p.source[k,j])
                chemical+=p.h[k,j]*p.source[k,j]
                cp+=u[k+1,j]*p.cp[k,j]/f.gas.MW[k]
            end
            isempty(f.fixed_temperature) && (r[1,j]-=(f.gas_multiplier-1)*
                _flame_timescale/(1000*p.rho[j]*cp)*chemical)
        end
    end
    _catalytic_surface_residual!(w.coverage_residual,f,u,w;previous,dt)
    for k in 1:n
        k==f.dependent_species && continue
        r[k+1,N]+=f.surface_multiplier*w.surface.production_rates[f.surface.n_surface+k]*f.gas.MW[k]
    end
    return r
end

function _catalytic_newton!(f,w;previous=nothing,dt=Inf,maxiters=50,tolerance=1e-10,loglevel=0)
    u,theta=f.state,f.coverages
    B,N=size(u); ns=length(theta); n=f.gas.n_species; bandwidth=2B-1
    r=similar(u); trial=similar(u); rt=similar(u); rs=zeros(ns); rst=zeros(ns)
    oldtheta=copy(theta)
    border=zeros(length(u),ns); C=zeros(ns,B); D=zeros(ns,ns)
    gas_weights=ones(size(u)); surface_weights=ones(ns)
    for iteration in 1:maxiters
        counterflow_residual!(r,f,u,w;previous,dt)
        rs .= w.coverage_residual
        residualnorm=max(norm(r,Inf),norm(rs,Inf))
        surface_converged=!f.coverage_enabled || previous!==nothing ||
            f.surface_multiplier*norm(w.surface.coverage_rates,Inf)<1e-6
        residualnorm<tolerance && surface_converged && return true
        oldtheta .= theta
        step=try
            J=_counterflow_jacobian!(f,u,w,r;previous,dt)
            # Gas Jacobian is banded; surface coverages form a small dense border.
            # Re-evaluate the base state after the colored gas perturbations.
            counterflow_residual!(r,f,u,w;previous,dt)
            rs .= w.coverage_residual
            fill!(border,0.)
            for k in 1:ns
                h=1e-7*max(abs(theta[k]),1e-5)
                theta[k]+=h
                _catalytic_surface_residual!(rst,f,u,w;previous,dt)
                D[:,k].=(rst.-rs)./h
                for i in 1:n
                    i==f.dependent_species && continue
                    # Only the final gas node is directly coupled to coverages.
                    row=(N-1)*B+i+1
                    border[row,k]=f.surface_multiplier*w.surface.production_rates[ns+i]*f.gas.MW[i]
                end
                theta[k]=oldtheta[k]
            end
            _catalytic_surface_residual!(rst,f,u,w;previous,dt)
            for k in 1:ns
                h=1e-7*max(abs(theta[k]),1e-5)
                for i in 1:n
                    i==f.dependent_species && continue
                    row=(N-1)*B+i+1
                    border[row,k]=(border[row,k]-f.surface_multiplier*w.surface.production_rates[ns+i]*f.gas.MW[i])/h
                end
            end
            for k in 1:B
                h=1e-7*max(abs(u[k,N]),k==1 ? .1 : k<=n+1 ? 1e-5 : .01)
                old=u[k,N]; u[k,N]+=h
                _catalytic_surface_residual!(rst,f,u,w;previous,dt)
                C[:,k].=(rst.-rs)./h
                u[k,N]=old
            end
            # Row equilibration is used only in the line-search norm. In
            # particular, a stiff adsorption equation must not hide progress
            # in the inlet composition and wall mass-balance equations.
            for column in 1:length(u)
                gas_weights[column]=1/max(abs(J[2bandwidth+1,column]),1e-8)
            end
            for k in 1:ns
                surface_weights[k]=1/max(maximum(abs,@view D[k,:]),1.)
            end
            _,pivots=LinearAlgebra.LAPACK.gbtrf!(bandwidth,bandwidth,length(u),J)
            rhs=hcat(-vec(r),border)
            LinearAlgebra.LAPACK.gbtrs!('N',bandwidth,bandwidth,length(u),J,pivots,rhs)
            x=@view rhs[:,1]; Z=@view rhs[:,2:end]
            rows=(N-1)*B+1:N*B
            dtheta=(D-C*Z[rows,:])\(-rs-C*x[rows])
            (reshape(x-Z*dtheta,size(u)),dtheta)
        catch e
            theta .= oldtheta
            e isa SingularException || e isa LinearAlgebra.LAPACKException || rethrow()
            return false
        end
        du,dtheta=step
        all(isfinite,du)&&all(isfinite,dtheta) || return false
        alpha=1.
        for j in 1:N,k in 1:B
            low=k==1 ? .2 : k<=n+1 ? -1e-7 : -1e5
            high=k==1 ? 6. : k<=n+1 ? 1.00001 : 1e5
            du[k,j]<0 && (alpha=min(alpha,.99*(u[k,j]-low)/(-du[k,j])))
            du[k,j]>0 && (alpha=min(alpha,.99*(high-u[k,j])/du[k,j]))
        end
        for k in 1:ns
            dtheta[k]<0 && (alpha=min(alpha,.99*(theta[k]+1e-14)/(-dtheta[k])))
            dtheta[k]>0 && (alpha=min(alpha,.99*(1.00001-theta[k])/dtheta[k]))
        end
        accepted=false
        merit=norm(r.*gas_weights)^2+norm(rs.*surface_weights)^2
        bounded_alpha=alpha
        for backtrack in 1:28
            @. trial=u+alpha*du
            @. theta=oldtheta+alpha*dtheta
            counterflow_residual!(rt,f,trial,w;previous,dt)
            rst .= w.coverage_residual
            if all(isfinite,rt)&&all(isfinite,rst) && norm(rt.*gas_weights)^2+norm(rst.*surface_weights)^2 < merit*(1-1e-4*alpha)
                u .= trial; accepted=true; break
            end
            alpha*=.5
        end
        loglevel>1 && (println("Catalytic Newton ",iteration," residual=",residualnorm," damping=",alpha," bound=",bounded_alpha);flush(stdout))
        if !accepted || alpha<=1e-12
            theta .= oldtheta
            return false
        end
    end
    return false
end

function _catalytic_steady!(f;max_time_steps=800,timestep=Ref(1e-6),loglevel=0)
    w=CounterflowWorkspace(f)
    _catalytic_newton!(f,w;loglevel) && return true
    dt=timestep[]
    for step in 1:max_time_steps
        previous=(gas=copy(f.state),coverages=copy(f.coverages))
        if _catalytic_newton!(f,w;previous,dt,maxiters=30)
            dt=min(dt*1.5,.1); timestep[]=dt
            if step%10==0
                loglevel>0 && (println("Catalytic transient ",step," dt=",dt);flush(stdout))
                _catalytic_newton!(f,w;loglevel) && return true
            end
        else
            f.state .= previous.gas; f.coverages .= previous.coverages
            dt*=.25; timestep[]=dt
            dt>1e-14 || return false
        end
    end
    return false
end

"Solve coupled native gas and surface equations, with adaptive axial refinement."
function solve!(f::CatalyticImpingingJet;refine_grid=true,ratio=10.,slope=.8,curve=.8,
        prune=0.,grid_min=1e-10,max_points=1200,max_time_steps=800,loglevel=0,
        initial_time_step=1e-6)
    isfinite(ratio)&&ratio>1 && 0<slope<=1 && 0<curve<=1 &&
        isfinite(prune)&&prune<=min(slope,curve) && isfinite(grid_min)&&grid_min>0 &&
        isfinite(initial_time_step)&&initial_time_step>0 || throw(ArgumentError("invalid catalytic solver controls"))
    f.converged=false; timestep=Ref(Float64(initial_time_step))
    f.transport_model==:mixture_averaged && !f.soret_enabled && f.flux_gradient_basis==:mole ||
        throw(ArgumentError("catalytic jet currently requires mixture-averaged mole-gradient transport without Soret"))
    if !isnothing(f.inlet_start)
        # A large inlet change can cross very stiff adsorption regimes. Use
        # composition continuation from the last native solution before mesh
        # refinement, keeping both reaction multipliers at their requested values.
        start=f.inlet_start; target=mole_fractions(f.gas,f.fuel_Y;basis=:mass)
        f.inlet_start=nothing
        for fraction in (.001,.01,.05,.1,.2,.4,.6,.8,1.)
            X=(1-fraction).*start.+fraction.*target
            f.fuel_Y=X.*f.gas.MW./dot(X,f.gas.MW)
            loglevel>0 && (println("Catalytic inlet continuation ",fraction);flush(stdout))
            _catalytic_steady!(f;max_time_steps,timestep,loglevel) || error("catalytic inlet continuation failed")
        end
    end
    for pass in 1:40
        loglevel>0 && (println("Catalytic solve on ",length(f.grid)," points; gas/surface multipliers ",f.gas_multiplier," / ",f.surface_multiplier);flush(stdout))
        _catalytic_steady!(f;max_time_steps,timestep,loglevel) ||
            error("native catalytic impinging jet failed on $(length(f.grid)) points")
        if !refine_grid || !_refine_axisymmetric!(f;ratio,slope,curve,prune,grid_min,max_points)
            f.converged=true
            return f
        end
    end
    error("catalytic refinement did not finish")
end

heat_release_rate(f::CatalyticImpingingJet)=f.gas_multiplier .* heat_release_rate(f.flow)

"Surface rates, wall mass fluxes and elemental balances in SI units."
function catalytic_wall_diagnostics(f::CatalyticImpingingJet)
    w=CounterflowWorkspace(f); r=similar(f.state)
    counterflow_residual!(r,f,f.state,w)
    ns=f.surface.n_surface; ng=f.surface.n_gas
    production=f.surface_multiplier.*w.surface.production_rates
    gas_production=production[ns+1:ns+ng]
    wall_flux=copy(w.properties.flux[:,end])
    return (coverages=copy(f.coverages),coverage_rates=f.surface_multiplier.*w.surface.coverage_rates,
        coverage_sum=sum(f.coverages),surface_production_rates=production,
        gas_production_rates=gas_production,diffusive_mass_flux=wall_flux,
        wall_species_residual=wall_flux.+f.state[end,end].*f.state[2:ng+1,end].+gas_production.*f.gas.MW,
        elemental_production=f.surface.elemental_matrix*production,
        total_mass_production=dot(f.surface.molecular_weights,production),
        surface_heat_release=-dot(w.surface.enthalpy,production),flow_residual=norm(r,Inf))
end

export CatalyticImpingingJet, set_catalytic_reactions!, set_catalytic_inlet!
export initialize_catalytic_coverages!, catalytic_wall_diagnostics
