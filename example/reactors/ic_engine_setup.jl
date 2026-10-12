using Arrhenius

const ENGINE_FREQUENCY = 50.0
const ENGINE_DISPLACEMENT = 0.5e-3
const ENGINE_CLEARANCE = ENGINE_DISPLACEMENT/19
const ENGINE_PISTON_AREA = π*0.083^2/4
const ENGINE_INJECTION_RATE = 3.2e-5/((365-350)/360/ENGINE_FREQUENCY)

engine_crank_angle(t) = mod(2π*ENGINE_FREQUENCY*t,4π)
engine_piston_speed(t) = -ENGINE_DISPLACEMENT/ENGINE_PISTON_AREA * π*ENGINE_FREQUENCY * sin(engine_crank_angle(t))
engine_volume(t) = ENGINE_CLEARANCE+ENGINE_DISPLACEMENT/2*(1-cos(engine_crank_angle(t)))
engine_gate(t,opening,closing) = mod(engine_crank_angle(t)-opening*π/180,4π) < mod((closing-opening)*π/180,4π)

"Continuous segments separated by the source's valve and injection switching angles."
function engine_switching_times(stop)
    stops = [0.0,float(stop)]
    for cycle in 0:ceil(Int,stop*ENGINE_FREQUENCY/2), angle in (18,198,350,365,522,702)
        t = (720cycle+angle)/360/ENGINE_FREQUENCY
        0 < t < stop && push!(stops,t)
    end
    return sort!(unique!(stops))
end

"Exact source reservoirs, cylinder geometry, flows, and sinusoidal piston."
function engine_network(gas;segment_time=nothing,time_offset=0.0)
    air = Dict("o2"=>1.0,"n2"=>3.76)
    cylinder = WellStirredReactor(gas;temperature=300.0,pressure=1.3e5,
        mole_fractions=air,volume=ENGINE_CLEARANCE)
    inlet = Reservoir(gas;temperature=300.0,pressure=1.3e5,mole_fractions=air)
    injector = Reservoir(gas;temperature=300.0,pressure=1600e5,mole_fractions=Dict("c12h26"=>1.0))
    outlet = Reservoir(gas;temperature=300.0,pressure=1.2e5,mole_fractions=air)
    ambient = Reservoir(gas;temperature=300.0,pressure=1e5,mole_fractions=air)
    if segment_time === nothing
        inlet_gate = t->engine_gate(t,-18,198)
        outlet_gate = t->engine_gate(t,522,18)
        injection = t->ENGINE_INJECTION_RATE*engine_gate(t,350,365)
    else
        inlet_open = engine_gate(segment_time,-18,198)
        outlet_open = engine_gate(segment_time,522,18)
        fuel_open = engine_gate(segment_time,350,365)
        inlet_gate = t->inlet_open
        outlet_gate = t->outlet_open
        injection = ENGINE_INJECTION_RATE*fuel_open
    end
    flows = (Valve(:inlet,:cylinder;K=1e-6,time_function=inlet_gate),
             MassFlowController(:injector,:cylinder;mdot=injection),
             Valve(:cylinder,:outlet;K=1e-6,time_function=outlet_gate))
    base = ReactorNetwork((cylinder=cylinder,inlet=inlet,injector=injector,outlet=outlet,ambient=ambient);flows)
    return MovingWallNetwork(base;walls=(MovingWall(:ambient,:cylinder;
        area=ENGINE_PISTON_AREA,velocity=t->engine_piston_speed(t+time_offset)),))
end

"""
    solve_ic_engine(mechanism; integrator, times=nothing, ...)

Compute all eight revolutions with native ideal-gas chemistry and moving-wall
energy equations. The prepared YAML/NPZ contains the source nDodecane_IG phase.
Default output uses native accepted states with at most one crank-angle degree
and 20 K between records. Pass observation `times` to sample another grid,
including the source's temperature-limited output times.
Each valve/injection transition starts a new continuous integration segment.
With default output, `integrals` uses every accepted state independently of
profile selection; `coarse_integrals` uses every second state to check quadrature.
"""
function solve_ic_engine(mechanism;integrator,
        times=nothing,reltol=1e-13,species_atol=1e-26,progress=false,on_segment=nothing)
    adaptive_output = times === nothing
    times = adaptive_output ? Float64[0.0] : Float64.(times)
    if !adaptive_output
        length(times)>1 && first(times)==0 && issorted(times) && all(diff(times).>0) ||
            throw(ArgumentError("finite increasing output times starting at zero required"))
        all(isfinite,times) || throw(ArgumentError("finite output times required"))
    end
    gas = CreateSolution(mechanism)
    report_model = engine_network(gas)
    state = moving_wall_state(report_model)
    absolute = fill(species_atol,length(state))
    absolute[gas.n_species+1] = 1e-10
    absolute[report_model.volume_indices[1]] = 1e-20
    ledger = report_model.ledger_offset
    absolute[ledger:ledger+2] .= [1e-15,1e-6,1e-6]
    stops = engine_switching_times(adaptive_output ? 8/ENGINE_FREQUENCY : last(times))
    states = Vector{Vector{Float64}}(undef,length(times))
    states[1] = copy(state)
    quadrature = zeros(4)
    coarse_quadrature = zeros(4)
    integration_points = 1
    for i in 1:length(stops)-1
        start,stop = stops[i],stops[i+1]
        # Local elapsed time retains sub-ULP startup resolution when injection
        # activates nearly absent species late in the global trajectory.
        model = engine_network(gas;segment_time=(start+stop)/2,time_offset=start)
        indices = adaptive_output ? Int[] : findall(t->start<t<=stop,times)
        save_times = adaptive_output ? Float64[] : sort!(unique!(vcat([0.0,stop-start],times[indices].-start)))
        sol = solve_moving_wall(model,(0.0,stop-start);integrator,initial_state=state,
            reltol,abstol=absolute,saveat=save_times,dt=min(1e-10,(stop-start)/100),
            save_everystep=adaptive_output,dense=false,maxiters=1000000)
        state = copy(sol.u[end])
        if adaptive_output
            # Integrals use every accepted state in each continuous segment.
            # Frozen gates give the correct one-sided flow at both endpoints;
            # the profile selection below does not affect these integrals.
            segment = ic_engine_observables((;model,gas,times=sol.t,states=sol.u))
            quadrature .+= engine_quadrature_terms(segment,gas;omit_initial=false,integral=engine_quadratic_integral)
            coarse = unique!(vcat(collect(1:2:length(sol.t)),length(sol.t)))
            coarse_quadrature .+= engine_quadrature_terms(segment,gas;indices=coarse,
                omit_initial=false,integral=engine_quadratic_integral)
            integration_points += length(sol.t)-1
            selected = engine_output_indices(sol.t,[u[gas.n_species+1] for u in sol.u])
            for j in selected[2:end]
                push!(times,start+sol.t[j])
                push!(states,copy(sol.u[j]))
            end
        else
            for index in indices
                local_time = times[index]-start
                j = searchsortedfirst(sol.t,local_time)
                sol.t[j] == local_time || error("integrator omitted a requested observation time")
                states[index] = copy(sol.u[j])
            end
        end
        progress && println("crank angle ",round(stop*ENGINE_FREQUENCY*360;digits=3),
                            " deg, T=",state[gas.n_species+1]," K")
        on_segment === nothing || on_segment(stop,copy(state))
    end
    quadrature_integrals = adaptive_output ? engine_integral_values(quadrature) : nothing
    coarse_integrals = adaptive_output ? engine_integral_values(coarse_quadrature) : nothing
    # The ODE already integrates pressure work at its full state tolerance.
    adaptive_output && (quadrature[2] = state[ledger+2]-1e5*(state[report_model.volume_indices[1]]-ENGINE_CLEARANCE))
    integrals = adaptive_output ? engine_integral_values(quadrature) : nothing
    return (;model=report_model,gas,times,states,integrals,quadrature_integrals,coarse_integrals,integration_points)
end

"Select accepted states without skipping a temperature excursion inside an output interval."
function engine_output_indices(times,temperatures;max_interval=1/(360*ENGINE_FREQUENCY),max_change=20.0)
    selected = [1]
    for j in 2:length(times)
        last_index = selected[end]
        if times[j]-times[last_index] > max_interval*(1+1e-12) ||
           abs(temperatures[j]-temperatures[last_index]) > max_change+1e-8
            j-1 > last_index || error("integrator steps exceed engine output limits")
            push!(selected,j-1)
        end
        if j>2 && (temperatures[j-1]-temperatures[j-2])*(temperatures[j]-temperatures[j-1])<0 &&
           abs(temperatures[j-1]-temperatures[selected[end]])>1e-6
            selected[end] == j-1 || push!(selected,j-1)
        end
    end
    selected[end] == length(times) || push!(selected,length(times))
    return selected
end

"Native source observables and independent cumulative mass/energy balances."
function ic_engine_observables(result)
    rhs = moving_wall_rhs(result.model)
    n,ns = length(result.times),result.gas.n_species
    keys = (:temperature,:pressure,:volume,:mass,:crank_angle,:mean_molecular_weight,:entropy_mass,
            :mdot_in,:mdot_fuel,:mdot_out,:work_rate,:heat_release_rate,:internal_energy,
            :mass_balance,:energy_balance,:integrated_work)
    output = Dict{String,Any}(String(k)=>zeros(n) for k in keys)
    output["time"] = result.times
    output["Y"],output["X"] = zeros(n,ns),zeros(n,ns)
    for (i,(t,u)) in enumerate(zip(result.times,result.states))
        d = moving_wall_diagnostics(rhs,u,t)
        s = d.nodes.cylinder
        # The prepared RHS has just evaluated this state's species enthalpies,
        # standard entropies, and homogeneous production rates.
        work = rhs.network_rhs.workspaces[1]
        inv_mw = sum(s.mass_fractions./result.gas.MW)
        X = s.mass_fractions./result.gas.MW./inv_mw
        entropy = sum(s.mass_fractions[k]/result.gas.MW[k] *
            (work.entropy[k]-Arrhenius.R*log(max(X[k],1e-300)*s.pressure/one_atm)) for k in 1:ns)
        heat_release = -sum(work.h_mole.*work.wdot)*s.volume
        values = (s.temperature,s.pressure,s.volume,s.mass,engine_crank_angle(t),1/inv_mw,entropy,
            d.mass_flow_rates[1],d.mass_flow_rates[2],d.mass_flow_rates[3],
            (s.pressure-1e5)*d.volume_rates[1],heat_release,s.mass*s.internal_energy,
            d.mass_balance,d.energy_balance,d.pressure_work_output-1e5*(s.volume-ENGINE_CLEARANCE))
        for (key,value) in zip(keys,values)
            output[String(key)][i] = value
        end
        output["Y"][i,:] = s.mass_fractions
        output["X"][i,:] = X
    end
    return output
end

engine_trapezoid(y,t) = sum((y[i]+y[i-1])*(t[i]-t[i-1])/2 for i in 2:length(t))
"Integrate consecutive quadratic interpolants on a nonuniform accepted-step grid."
function engine_quadratic_integral(y,t)
    n = length(t)
    n < 3 && return engine_trapezoid(y,t)
    value = zero(eltype(y))
    for i in 1:2:n-2
        a,b = t[i+1]-t[i],t[i+2]-t[i+1]
        span = a+b
        value += span/6*((2-b/a)*y[i]+span^2/(a*b)*y[i+1]+(2-a/b)*y[i+2])
    end
    iseven(n) && (value += (y[end]+y[end-1])*(t[end]-t[end-1])/2)
    return value
end

function engine_quadrature_terms(output,gas;indices=eachindex(output["time"]),omit_initial=true,
        integral=engine_trapezoid)
    # Source estimates omit the initial state from their sampled integrals.
    selected = omit_initial ? indices[2:end] : indices
    t = output["time"][selected]
    heat = integral(output["heat_release_rate"][selected],t)
    work = integral(output["work_rate"][selected],t)
    co = findfirst(==("co"),lowercase.(gas.species_names))
    co === nothing && error("engine mechanism has no CO")
    weights = output["mean_molecular_weight"][selected].*output["mdot_out"][selected]
    return [heat,work,integral(weights.*output["X"][selected,co],t),integral(weights,t)]
end
engine_integral_values(q) = (;heat_J=q[1],work_J=q[2],efficiency=q[2]/q[1],CO_ppm=1e6q[3]/q[4])
"Source-style trapezoidal estimates on a supplied profile grid; full-step values are in `result.integrals`."
ic_engine_integrals(output,gas) = engine_integral_values(engine_quadrature_terms(output,gas))

"Accepted-state quadrature and independent display/conservation diagnostics."
function ic_engine_summary(result,output)
    result.integrals === nothing && throw(ArgumentError("summary requires native adaptive output"))
    energy_scale = max(maximum(abs,output["integrated_work"]),maximum(abs,output["internal_energy"]),1.0)
    relative_quadrature_change = Dict(String(k)=>abs(getproperty(result.coarse_integrals,k)/v-1)
        for (k,v) in pairs(result.quadrature_integrals))
    return Dict("integration_points"=>result.integration_points,"display_points"=>length(result.times),
        "end_time_s"=>last(result.times),"maximum_output_interval_s"=>maximum(diff(result.times)),
        "maximum_output_temperature_change_K"=>maximum(abs,diff(output["temperature"])),
        "volume_identity_error_m3"=>maximum(abs,output["volume"]-engine_volume.(result.times)),
        "mass_balance_relative_drift"=>maximum(abs,output["mass_balance"].-output["mass_balance"][1])/output["mass"][1],
        "energy_balance_relative_drift"=>maximum(abs,output["energy_balance"].-output["energy_balance"][1])/energy_scale,
        "minimum_mass_fraction"=>minimum(output["Y"]),
        "integrals"=>Dict(String(k)=>v for (k,v) in pairs(result.integrals)),
        "accepted_state_quadrature"=>Dict(String(k)=>v for (k,v) in pairs(result.quadrature_integrals)),
        "every_second_state_integrals"=>Dict(String(k)=>v for (k,v) in pairs(result.coarse_integrals)),
        "relative_quadrature_change"=>relative_quadrature_change,
        "integrated_work_ledger_J"=>last(output["integrated_work"]),
        "work_quadrature_ledger_relative_error"=>abs(result.quadrature_integrals.work_J/last(output["integrated_work"])-1))
end

function write_ic_engine_csv(path,output,gas)
    scalar_keys = ("time","crank_angle","temperature","pressure","volume","mass","entropy_mass",
        "mdot_in","mdot_fuel","mdot_out","work_rate","heat_release_rate","integrated_work")
    species = ("o2","co2","co","c12h26")
    indices = [findfirst(==(name),lowercase.(gas.species_names)) for name in species]
    open(path,"w") do io
        println(io,join(vcat(collect(scalar_keys),["X_"*name for name in species]),','))
        for i in eachindex(output["time"])
            row = [output[key][i] for key in scalar_keys]
            append!(row,output["X"][i,indices])
            println(io,join(row,','))
        end
    end
    return path
end
