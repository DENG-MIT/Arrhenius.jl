# Exact Jacobian for the unchanged source engine topology.
# Only chemical derivatives use AD; flow, caloric, volume and ledger chains
# differentiate the shared MovingWallRHS intermediates directly.
function engine_ad_unseeded_jacobian(problem;chunk=8,rate_evaluator=wdot!,clip_trials=true)
    return engine_ad_unseeded_jacobian(problem,Val(chunk),rate_evaluator,clip_trials)
end
function engine_ad_unseeded_jacobian(problem,::Val{N},rate_evaluator,clip_trials) where N
    rhs=problem.f; model=rhs.model; network=model.network
    keys(network.nodes)==(:cylinder,:inlet,:injector,:outlet,:ambient) || error("source engine nodes required")
    network.offsets==[1,0,0,0,0] || error("one source cylinder required")
    isempty(network.walls) && length(model.walls)==1 || error("source piston required")
    wall=only(model.walls)
    wall isa MovingWall && wall.K==wall.U==0 && wall.heat_flux==0 || error("prescribed adiabatic piston required")
    inlet,fuel,outlet=network.flows
    inlet isa Valve && outlet isa Valve && fuel isa MassFlowController && fuel.mdot isa Real || error("frozen source flow gates required")
    inlet.pressure_function===identity && outlet.pressure_function===identity || error("linear source valves required")
    gas=network.nodes.cylinder.initial.gas; ns=gas.n_species
    network.nodes.cylinder.initial.energy===:adiabatic || error("adiabatic source cylinder required")
    isempty(gas.reaction.blowers_masel.reaction_indices) || error("Blowers-Masel outside engine derivative scope")
    state=rhs.network_rhs.state_vector[1]; floats=rhs.network_rhs.workspaces[1]
    temperature=Ref(problem.u0[ns+1]); concentrations=zeros(ns); source=zeros(ns)
    configuration=ForwardDiff.JacobianConfig(nothing,source,concentrations,ForwardDiff.Chunk{N}())
    D=eltype(typeof(configuration))
    cache=Arrhenius._KineticsTemperatureCache(gas.reaction)
    kinetic=KineticsWorkspace(gas.reaction,D)
    composition_source! = function (out,C)
        rate_evaluator(out,gas.reaction,temperature[],C,floats.entropy,floats.h_mole,kinetic;
            temperature_cache=cache,rate_multipliers=network.nodes.cylinder.initial.rate_multipliers)
        nothing
    end
    concentration_jacobian=zeros(ns,ns)
    DT=EngineTemperatureDual
    temp_kinetic=KineticsWorkspace(gas.reaction,DT)
    td_h=zeros(DT,ns);td_cp=zeros(DT,ns);td_s=zeros(DT,ns);td_source=zeros(DT,ns)
    derivative=zero(problem.u0)
    cv_species=zeros(ns);internal_species=zeros(ns)
    return function (J,u,p,t)
        rhs(derivative,u,p,t)
        fill!(J,0)
        T=state.temperature;P=state.pressure;V=state.volume;m=state.mass
        temperature[]=T
        copyto!(concentrations,floats.C)
        # Mass-action right derivatives at absent species; actual negative
        # mass columns are zeroed by the chain rule below, matching the RHS.
        for k in 1:ns
            clip_trials && concentrations[k]==0 && (concentrations[k]=1e-100)
        end
        rate_evaluator(source,gas.reaction,T,floats.C,floats.entropy,floats.h_mole,floats.kinetics;
            temperature_cache=cache,rate_multipliers=network.nodes.cylinder.initial.rate_multipliers)
        ForwardDiff.jacobian!(concentration_jacobian,composition_source!,source,concentrations,configuration)
        Td=engine_temperature_dual(T)
        Arrhenius.cal_cp_R!(td_cp,gas,Td,P,floats.X)
        Arrhenius.cal_h_RT!(td_h,gas,Td,P,floats.X)
        Arrhenius.cal_s0_R!(td_s,gas,Td,P,floats.X)
        capacity_temperature=0.0
        for k in 1:ns
            td_h[k]*=R*Td;td_s[k]*=R
            cv_species[k]=(floats.cp_R[k]-1)*R/gas.MW[k]
            internal_species[k]=(floats.h_mole[k]-R*T)/gas.MW[k]
            capacity_temperature+=u[k]*R/gas.MW[k]*ForwardDiff.partials(td_cp[k])[1]
        end
        rate_evaluator(td_source,gas.reaction,Td,floats.C,td_s,td_h,temp_kinetic;
            rate_multipliers=network.nodes.cylinder.initial.rate_multipliers)
        for k in 1:ns
            for j in 1:ns
                J[k,j]=(!clip_trials || u[j]>=0) ? gas.MW[k]/gas.MW[j]*concentration_jacobian[k,j] : 0.0
            end
            J[k,ns+1]=V*gas.MW[k]*ForwardDiff.partials(td_source[k])[1]
            chain=0.0
            for j in 1:ns
                chain+=concentration_jacobian[k,j]*floats.C[j]
            end
            J[k,ns+2]=gas.MW[k]*(floats.wdot[k]-chain)
        end
        sinlet,sfuel=rhs.network_rhs.state_vector[2],rhs.network_rhs.state_vector[3]
        inlet_slope=sinlet.pressure>P ? -inlet.K*inlet.time_function(t) : 0.0
        outlet_slope=P>rhs.network_rhs.state_vector[4].pressure ? outlet.K*outlet.time_function(t) : 0.0
        exhaust=rhs.network_rhs.mass_flow_rates[3]
        Vdot=rhs.volume_rates[1];capacity=m*state.cv
        Tdot=derivative[ns+1];ledger=model.ledger_offset
        for j in 1:ns+2
            dp=j<=ns ? R*T/(V*gas.MW[j]) : j==ns+1 ? P/T : -P/V
            dfeed=inlet_slope*dp;dexhaust=outlet_slope*dp
            dh=j<=ns ? (floats.h_mole[j]/gas.MW[j]-state.enthalpy)/m : j==ns+1 ? state.cp : 0.0
            energy_gradient=dfeed*sinlet.enthalpy-dexhaust*state.enthalpy-exhaust*dh
            composition_energy_gradient=0.0
            for k in 1:ns
                dy=j<=ns ? ((k==j)-state.mass_fractions[k])/m : 0.0
                J[k,j]+=dfeed*sinlet.mass_fractions[k]-dexhaust*state.mass_fractions[k]-exhaust*dy
                composition_energy_gradient+=internal_species[k]*J[k,j]
                j==ns+1 && (composition_energy_gradient+=cv_species[k]*derivative[k])
            end
            dc=j<=ns ? cv_species[j] : j==ns+1 ? capacity_temperature : 0.0
            J[ns+1,j]=(energy_gradient-composition_energy_gradient-dp*Vdot-Tdot*dc)/capacity
            J[ledger,j]=dfeed-dexhaust
            J[ledger+1,j]=energy_gradient
            J[ledger+2,j]=dp*Vdot
        end
        nothing
    end
end
