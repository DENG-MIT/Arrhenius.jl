"""
    SurfaceFlowReactor(surface; temperature, pressure=one_atm, mole_fractions,
                       area, mass_flow_rate, surface_area_per_length)

An isothermal, constant-area, frictionless ideal-gas plug-flow reactor. Surface
coverages satisfy stationary algebraic equations at every axial position.
`surface_area_per_length` is catalyst area per axial metre, in m²/m. Persistent
bulk phases are not supported here: total gas mass flux remains constant.
"""
struct SurfaceFlowReactor{G}
    surface::SurfaceMechanism
    gas::G
    temperature::Float64
    pressure::Float64
    area::Float64
    mass_flow_rate::Float64
    surface_area_per_length::Float64
    mass_fractions::Vector{Float64}
    mass_flux::Float64
    momentum_flux::Float64
end
function SurfaceFlowReactor(surface::SurfaceMechanism; temperature,pressure=one_atm,
        mole_fractions,area,mass_flow_rate,surface_area_per_length,
        gas=CreateSolution(surface.gas_file))
    isempty(surface.bulk_molar_volumes) || throw(ArgumentError("flow model requires no bulk mass exchange"))
    all(x->isfinite(x) && x>0,(temperature,pressure,area,mass_flow_rate,surface_area_per_length)) ||
        throw(ArgumentError("positive finite flow parameters required"))
    names = surface.species_names[surface.n_surface+1:surface.n_surface+surface.n_gas]
    gas.species_names == names || throw(ArgumentError("companion gas species order required"))
    isapprox(gas.MW,surface.molecular_weights[surface.n_surface+1:end];rtol=1e-12) ||
        throw(ArgumentError("companion gas molecular weights required"))
    inlet = IdealGasReactor(gas;temperature,pressure,mole_fractions,energy=:isothermal)
    G = mass_flow_rate/area
    G^2 < inlet.density*pressure || throw(ArgumentError("surface flow inlet must be on the subsonic branch"))
    momentum = pressure + G^2/inlet.density
    return SurfaceFlowReactor(surface,gas,Float64(temperature),Float64(pressure),Float64(area),
        Float64(mass_flow_rate),Float64(surface_area_per_length),
        Float64.(inlet.mass_fractions),Float64(G),Float64(momentum))
end

"Recover the subsonic pressure, density, and speed from conserved momentum and mass flux."
function surface_flow_properties(flow::SurfaceFlowReactor,Y)
    length(Y) == flow.gas.n_species || throw(DimensionMismatch("flow gas composition"))
    MW = 1/sum(Y./flow.gas.MW)
    discriminant = flow.momentum_flux^2-4flow.mass_flux^2*R*flow.temperature/MW
    discriminant > 0 || throw(DomainError(discriminant,"flow reached the sonic limit"))
    P = (flow.momentum_flux+sqrt(discriminant))/2
    rho = P*MW/(R*flow.temperature)
    return (pressure=P,density=rho,speed=flow.mass_flux/rho,
            mass_flux=flow.mass_flux,momentum_flux=flow.momentum_flux,molecular_weight=MW)
end

struct SurfaceFlowRHS{F}
    flow::F
    surface_workspace::SurfaceWorkspace
    gas_workspace::ReactorWorkspace{Float64}
    trial::Vector{Float64}
    plus::Vector{Float64}
    minus::Vector{Float64}
    base::Vector{Float64}
end
function surface_flow_rhs(flow::SurfaceFlowReactor)
    n = flow.gas.n_species+flow.surface.n_surface
    return SurfaceFlowRHS(flow,SurfaceWorkspace(flow.surface),ReactorWorkspace(flow.gas,Float64),
                         zeros(n),zeros(n),zeros(n),zeros(n))
end
function (rhs::SurfaceFlowRHS)(f,u,p,z)
    flow,m,w = rhs.flow,rhs.flow.surface,rhs.gas_workspace
    ng,ns = m.n_gas,m.n_surface
    length(u) == length(f) == ng+ns || throw(DimensionMismatch("surface flow state"))
    Y,theta = @view(u[1:ng]),@view(u[ng+1:end])
    state = surface_flow_properties(flow,Y)
    inverseMW = sum(Y./flow.gas.MW)
    for k in 1:ng
        w.X[k] = Y[k]/(flow.gas.MW[k]*inverseMW)
        w.C[k] = max(Y[k],0)*state.density/flow.gas.MW[k]
    end
    T,P = flow.temperature,state.pressure
    surface_rates!(rhs.surface_workspace,m,T,P,w.X,theta)
    if flow.gas.n_reactions > 0
        cal_h_RT!(w.h_mole,flow.gas,T,P,w.X)
        cal_s0_R!(w.entropy,flow.gas,T,P,w.X)
        w.h_mole .*= R*T
        w.entropy .*= R
        wdot!(w.wdot,flow.gas.reaction,T,w.C,w.entropy,w.h_mole,w.kinetics)
    else
        fill!(w.wdot,0.0)
    end
    ratio = flow.surface_area_per_length/flow.area
    source_mass = 0.0
    for k in 1:ng
        f[k] = flow.gas.MW[k]*(w.wdot[k]+ratio*rhs.surface_workspace.production_rates[ns+k])
        source_mass += f[k]
    end
    # This correction vanishes on the stationary-surface constraint and preserves
    # sum(Y) for the nonlinear solver's intermediate off-constraint states.
    for k in 1:ng
        f[k] = (f[k]-Y[k]*source_mass)/flow.mass_flux
    end
    f[ng+1] = sum(theta)-1
    for k in 2:ns
        f[ng+k] = rhs.surface_workspace.production_rates[k]
    end
    return nothing
end

function surface_flow_jacobian!(J,u,rhs::SurfaceFlowRHS,z=0.0)
    n = length(u)
    size(J) == (n,n) || throw(DimensionMismatch("surface flow Jacobian"))
    copyto!(rhs.trial,u)
    rhs(rhs.base,u,nothing,z)
    for j in 1:n
        h = cbrt(eps(Float64))*max(abs(u[j]),1e-8)
        rhs.trial[j] = u[j]+h
        rhs(rhs.plus,rhs.trial,nothing,z)
        if u[j] > h
            rhs.trial[j] = u[j]-h
            rhs(rhs.minus,rhs.trial,nothing,z)
            for i in 1:n
                J[i,j] = (rhs.plus[i]-rhs.minus[i])/(2h)
            end
        else
            rhs.trial[j] = u[j]+2h
            rhs(rhs.minus,rhs.trial,nothing,z)
            for i in 1:n
                J[i,j] = (-3rhs.base[i]+4rhs.plus[i]-rhs.minus[i])/(2h)
            end
        end
        rhs.trial[j] = u[j]
    end
    return J
end

function _surface_cstr_polish(flow,Yin,volume,area,Y0,theta0;tolerance=1e-9,maxiters=25)
    gas,m = flow.gas,flow.surface
    ng,ns = gas.n_species,m.n_surface
    w,sw = ReactorWorkspace(gas,Float64),SurfaceWorkspace(m)
    u = vcat(Y0,theta0)
    pivot = argmax(Yin)
    function residual!(f,u)
        Y,theta = @view(u[1:ng]),@view(u[ng+1:end])
        inverseMW = sum(Y./gas.MW)
        rho = flow.pressure/(R*flow.temperature*inverseMW)
        for k in 1:ng
            w.X[k] = Y[k]/(gas.MW[k]*inverseMW)
            w.C[k] = max(Y[k],0)*rho/gas.MW[k]
        end
        surface_rates!(sw,m,flow.temperature,flow.pressure,w.X,theta)
        if gas.n_reactions > 0
            cal_h_RT!(w.h_mole,gas,flow.temperature,flow.pressure,w.X)
            cal_s0_R!(w.entropy,gas,flow.temperature,flow.pressure,w.X)
            w.h_mole .*= R*flow.temperature; w.entropy .*= R
            wdot!(w.wdot,gas.reaction,flow.temperature,w.C,w.entropy,w.h_mole,w.kinetics)
        else
            fill!(w.wdot,0)
        end
        for k in 1:ng
            f[k] = Yin[k]-Y[k]+gas.MW[k]*(volume*w.wdot[k]+area*sw.production_rates[ns+k])/flow.mass_flow_rate
        end
        f[pivot] = sum(Y)-1
        for k in 1:ns
            f[ng+k] = sw.coverage_rates[k]
        end
        f[ng+1] = sum(theta)-1
        return nothing
    end
    f,plus,minus,trial,scales = similar(u),similar(u),similar(u),copy(u),similar(u)
    J = zeros(length(u),length(u))
    normscaled = Inf
    for _ in 1:maxiters
        residual!(f,u)
        for j in eachindex(u)
            h = cbrt(eps(Float64))*max(abs(u[j]),1e-8)
            trial .= u; trial[j] += h; residual!(plus,trial)
            if u[j] > h
                trial[j] = u[j]-h; residual!(minus,trial)
                J[:,j] .= (plus.-minus)./(2h)
            else
                trial[j] = u[j]+2h; residual!(minus,trial)
                J[:,j] .= (-3f.+4plus.-minus)./(2h)
            end
        end
        for i in eachindex(u)
            scales[i] = i<=ng || i==ng+1 ? 1.0 :
                max(sum(abs(J[i,j])*max(abs(u[j]),1e-12) for j in eachindex(u)),1.0)
        end
        normscaled = maximum(abs.(f)./scales)
        normscaled <= tolerance && return (converged=true,Y=copy(u[1:ng]),theta=copy(u[ng+1:end]),residual=normscaled)
        step = try
            -(J./scales)\(f./scales)
        catch error
            error isa SingularException || rethrow()
            return (converged=false,Y=Y0,theta=theta0,residual=normscaled)
        end
        all(isfinite,step) || return (converged=false,Y=Y0,theta=theta0,residual=normscaled)
        alpha = 1.0
        for k in eachindex(u)
            step[k] < 0 && (alpha=min(alpha,-0.99max(u[k],0.0)/step[k]))
        end
        accepted = false
        for _ in 1:20
            @. trial = u+alpha*step
            residual!(plus,trial)
            if maximum(abs.(plus)./scales) < normscaled
                u .= trial;accepted=true;break
            end
            alpha *= 0.5
        end
        accepted || break
    end
    return (converged=false,Y=Y0,theta=theta0,residual=normscaled)
end

"""
    surface_flow_problem(flow, length; initial_coverages)

Return a constant-mass-matrix DAE in axial distance. The first gas-species rows
are differential and the remaining coverage rows are algebraic. Supply
stationary inlet coverages, or use `solve_surface_flow` to initialize them.
"""
function surface_flow_problem(flow::SurfaceFlowReactor,length;initial_coverages)
    isfinite(length) && length>0 || throw(ArgumentError("positive finite bed length required"))
    theta = _surface_composition(flow.surface.species_names[1:flow.surface.n_surface],initial_coverages)
    rhs = surface_flow_rhs(flow)
    ng,ns = flow.surface.n_gas,flow.surface.n_surface
    check = surface_flow_diagnostics(flow,vcat(flow.mass_fractions,theta))
    check.stationary_relative_residual <= 1e-7 ||
        throw(ArgumentError("stationary inlet coverages required; use solve_surface_flow for initialization"))
    return (f=rhs,jac=(J,u,p,z)->surface_flow_jacobian!(J,u,rhs,z),
        tgrad=(du,u,p,z)->(fill!(du,0.0);nothing),
        mass_matrix=Diagonal(vcat(ones(ng),zeros(ns))),
        u0=vcat(flow.mass_fractions,theta),tspan=(0.0,Float64(length)),p=nothing,
        isoutofdomain=(u,p,z)->any(x->!isfinite(x)||x < -1e-12,u))
end

"Solve the stationary surface DAE with caller-owned native stiff integrators."
function solve_surface_flow(flow::SurfaceFlowReactor,length;integrator,coverage_integrator,
        initial_coverages=nothing,coverage_options=(;reltol=1e-10,abstol=1e-18),kwargs...)
    if isnothing(initial_coverages)
        surface = IdealSurface(flow.surface;temperature=flow.temperature,pressure=flow.pressure,
            mole_fractions=Y2X(flow.gas,flow.mass_fractions))
        transient = solve_surface(surface,(0.0,1.0);integrator=coverage_integrator,coverage_options...)
        initial_coverages = max.(transient.u[end],0.0)
        initial_coverages ./= sum(initial_coverages)
        rhs = surface_rhs(surface)
        J,f = zeros(flow.surface.n_surface,flow.surface.n_surface),similar(initial_coverages)
        converged = false
        for _ in 1:20
            rhs(f,initial_coverages,nothing,0.0)
            if maximum(abs,f) < 1e-6 && abs(sum(initial_coverages)-1) < 1e-12
                converged = true
                break
            end
            surface_jacobian!(J,initial_coverages,rhs)
            J[1,:] .= 1.0
            f[1] = sum(initial_coverages)-1
            step = -(J\f)
            alpha = 1.0
            for k in eachindex(step)
                step[k] < 0 && (alpha=min(alpha,-0.99initial_coverages[k]/step[k]))
            end
            initial_coverages .+= alpha.*step
        end
        converged || error("surface-flow inlet coverages did not satisfy stationary constraints")
    end
    return integrator(surface_flow_problem(flow,length;initial_coverages);kwargs...)
end

"""
    surface_reactor_chain(flow, length; reactors=201, porosity=0.3,
                         integrator, kwargs...)

March a chain of isothermal stirred reactors to steady state. Each volume is
`area*length/(reactors-1)*porosity`; its surface area is
`surface_area_per_length*length/(reactors-1)`. As in `surf_pfr_chain.py`, the
first reported state is the first reactor outlet at coordinate zero. Returns
all gas mole/mass fractions, coverages, pressures, and dimensionless steady
residuals. A damped Newton solve enforces fixed-pressure species and site
balances; native stiff integration supplies a stable fallback initial guess.
"""
function surface_reactor_chain(flow::SurfaceFlowReactor,length;reactors=201,porosity=0.3,
        integrator,steady_tolerance=1e-9,max_time=100.0,kwargs...)
    reactors >= 2 && isfinite(length) && length>0 && 0<porosity<=1 ||
        throw(ArgumentError("invalid reactor-chain geometry"))
    gas,m = flow.gas,flow.surface
    dz = length/(reactors-1)
    volume,area = flow.area*dz*porosity,flow.surface_area_per_length*dz
    Y,theta = copy(flow.mass_fractions),copy(m.initial_coverages)
    ng,ns = gas.n_species,m.n_surface
    ys,xs,thetas = zeros(ng,reactors),zeros(ng,reactors),zeros(ns,reactors)
    pressure,residuals = zeros(reactors),zeros(reactors)
    initial_pressure = flow.pressure
    for a in 1:reactors
        initial = IdealGasReactor(gas;temperature=flow.temperature,pressure=initial_pressure,
                                 mass_fractions=Y,constraint=:constant_volume,energy=:isothermal)
        vessel = WellStirredReactor(initial;volume)
        exhaust = IdealGasReactor(gas;temperature=flow.temperature,pressure=flow.pressure,
                                  mass_fractions=Y,energy=:isothermal)
        feed = MassFlowController(:inlet,:reactor;mdot=flow.mass_flow_rate)
        outlet = PressureController(:reactor,:outlet;primary=feed,K=1e-6)
        network = ReactorNetwork((inlet=Reservoir(initial),reactor=vessel,outlet=Reservoir(exhaust));
                                 flows=(feed,outlet))
        system = CatalyticNetwork(network;surfaces=(ReactorSurface(:reactor,m;area,coverages=theta),))
        rhs = catalytic_rhs(system)
        u = catalytic_state(system)
        derivative = similar(u)
        time = 0.0
        interval = max(10initial.density*volume/flow.mass_flow_rate,1e-4)
        residual = Inf
        polished = _surface_cstr_polish(flow,Y,volume,area,Y,theta;tolerance=steady_tolerance)
        if polished.converged
            residual = polished.residual
        end
        while !polished.converged && time < max_time
            stop = min(time+interval,max_time)
            solution = solve_catalytic(system,(time,stop);integrator,initial_state=u,kwargs...)
            u = copy(solution.u[end])
            rhs(derivative,u,nothing,stop)
            guessY = max.(u[1:ng],0.0);guessY ./= sum(guessY)
            guessTheta = max.(u[ng+2:end],0.0);guessTheta ./= sum(guessTheta)
            polished = _surface_cstr_polish(flow,Y,volume,area,guessY,guessTheta;tolerance=steady_tolerance)
            residual = polished.residual
            time = stop
            polished.converged && break
            interval *= 2
        end
        residual <= steady_tolerance || error("surface chain reactor $a did not reach steady state")
        Y = max.(polished.Y,0.0); Y ./= sum(Y)
        theta = max.(polished.theta,0.0);theta ./= sum(theta)
        initial_pressure = flow.pressure
        ys[:,a],xs[:,a],thetas[:,a] = Y,Y2X(gas,Y),theta
        pressure[a],residuals[a] = initial_pressure,residual
    end
    return (distance=collect(range(0.0,Float64(length);length=reactors)),mass_fractions=ys,
        mole_fractions=xs,coverages=thetas,pressure=pressure,residuals=residuals,
        mass_flow_rate=flow.mass_flow_rate,species_mass_flux=flow.mass_flux.*ys)
end

"Mass, momentum, elemental fluxes, and stationary-surface residual at one axial state."
function surface_flow_diagnostics(flow::SurfaceFlowReactor,u)
    ng,ns = flow.surface.n_gas,flow.surface.n_surface
    length(u) == ng+ns || throw(DimensionMismatch("surface flow state"))
    Y,theta = u[1:ng],u[ng+1:end]
    state = surface_flow_properties(flow,Y)
    X = Y2X(flow.gas,Y)
    w = surface_rates!(SurfaceWorkspace(flow.surface),flow.surface,flow.temperature,
                      state.pressure,X,theta)
    gross = abs.(flow.surface.stoichiometry[1:ns,:])*(w.forward_rates+w.reverse_rates)
    stationary = maximum(abs,@view(w.production_rates[1:ns]))/max(maximum(gross),1e-300)
    molar_flux = flow.mass_flux.*Y./flow.gas.MW
    return (;state...,species_mass_flux=flow.mass_flux.*Y,species_molar_flux=molar_flux,
        elemental_flux=flow.gas.ele_matrix*molar_flux,mass_fraction_sum=sum(Y),
        coverage_sum=sum(theta),stationary_relative_residual=stationary,
        coverage_rate=copy(w.coverage_rates))
end

export SurfaceFlowReactor, SurfaceFlowRHS, surface_flow_properties, surface_flow_rhs
export surface_flow_jacobian!, surface_flow_problem, solve_surface_flow, surface_reactor_chain
export surface_flow_diagnostics
