mutable struct EngineQNDFAdapter{J,R,A,C,D}
    records::Vector{Dict{String,Any}}
    saved::Dict{String,Any}
    signed_jacobian_factory::J
    rhs_observer::R
    accepted_observer::A
    coordinate_factory::C
    diagnostics::D
end
EngineQNDFAdapter(;signed_jacobian_factory=nothing,rhs_observer=nothing,accepted_observer=nothing,coordinate_factory=nothing,diagnostics=EngineDiagnosticLog())=
    EngineQNDFAdapter(Dict{String,Any}[],Dict{String,Any}(),
        signed_jacobian_factory,rhs_observer,accepted_observer,coordinate_factory,diagnostics)

function (adapter::EngineQNDFAdapter)(problem;reltol,abstol,kwargs...)
    reltol==1e-12 || error("source relative tolerance required")
    index=length(adapter.records)+1
    stops=engine_switching_times(.16)
    start,stop=stops[index],stops[index+1]
    problem.tspan==(0.,stop-start) || error("shared source interval mismatch")
    signed=signed_engine_problem(problem)
    if adapter.signed_jacobian_factory!==nothing
        signed=merge(signed,(;jac=adapter.signed_jacobian_factory(signed)))
    end
    scaled,scale=engine_scaled_problem(signed,signed.jac)
    ode_rhs=adapter.rhs_observer===nothing ? scaled.f :
        adapter.rhs_observer(scaled.f,scale,index,start,stop)
    model=problem.f.model;gas=model.network.nodes.cylinder.initial.gas;ns=gas.n_species
    source_atol=copy(abstol)
    source_atol[ns+1]=1e-16
    source_atol[model.volume_indices[1]]=1e-16
    physical_mass=sum(problem.u0[1:ns])
    dynamic=model.network.flows[3].time_function(0.)!=0
    weights=EngineMassWeights(ns,scale,physical_mass;dynamic)
    scaled_atol=source_atol./scale
    engine_set_mass_atol!(weights,scaled_atol)
    coordinates=adapter.coordinate_factory===nothing ? nothing : adapter.coordinate_factory(problem)
    coordinates!==nothing && adapter.accepted_observer!==nothing && error("full-state accepted observers cannot use restricted coordinates")
    parts=coordinates===nothing ? (;f=ode_rhs,jac=scaled.jac,tgrad=scaled.tgrad,
        u0=scaled.u0,atol=scaled_atol,options=(;)) :
        (;f=EngineRestrictedRHS(ode_rhs,coordinates),jac=engine_restricted_jacobian(coordinates,scaled.jac),
            tgrad=engine_restricted_tgrad(coordinates,scaled.tgrad),u0=scaled.u0[coordinates.indices],
            atol=scaled_atol[coordinates.indices],options=(;internalnorm=FullEngineRMS(length(scale))))
    initial_reference=weights.reference
    previous_temperature=Ref(scaled.u0[ns+1])
    first_accepted=Ref(UInt64(0))
    function outside(u,p,t)
        full=engine_expand!(coordinates,u)
        scaled.isoutofdomain(full,p,t) ||
            abs(full[ns+1]-previous_temperature[])>20/scale[ns+1] || engine_mass_guard!(weights,full)
    end
    function accepted!(integrator)
        stamp=time_ns()
        first_accepted[]==0 && (first_accepted[]=stamp)
        (stamp-first_accepted[])/1e9<60 || error("one engine interval exceeded 60 s after its first accepted state")
        adapter.accepted_observer===nothing ||
            adapter.accepted_observer(integrator,scale,index,start)
        if coordinates===nothing
            engine_accept_mass_weights!(weights,integrator)
        else
            engine_accept_restricted_weights!(weights,integrator,coordinates,scaled_atol)
        end
        previous_temperature[]=engine_expand!(coordinates,integrator.u)[ns+1]
        SciMLBase.u_modified!(integrator,false)
        nothing
    end
    limiter=DiscreteCallback((u,t,integrator)->true,accepted!;save_positions=(false,false))
    ode=ODEProblem(ODEFunction(parts.f;jac=parts.jac,tgrad=parts.tgrad),
        parts.u0,scaled.tspan,scaled.p)
    trace=Dict{String,Any}("initial_state"=>copy(problem.u0),"interval"=>[start,stop],
        "scale"=>scale,"integrated_initial_state"=>copy(parts.u0),"initial_scaled_atol"=>copy(parts.atol),
        "retained_indices"=>(coordinates===nothing ? collect(eachindex(scale)) : coordinates.indices),
        "rms_norm_denominator"=>[length(scale)])
    tracename="raw-interval"*string(index)
    engine_emit!(adapter.diagnostics,tracename,trace)
    reconstruct=u->engine_expand!(coordinates,u).*scale
    timed=@timed engine_probe_solve(ode,QNDF(),trace,adapter.diagnostics,tracename,reconstruct;time_offset=start,
        reltol,abstol=parts.atol,isoutofdomain=outside,callback=limiter,
        dtmax=1/(360*50),parts.options...,kwargs...)
    sol=timed.value
    states=[engine_expand!(coordinates,u).*scale for u in sol.u]
    trace["saved_times"]=start.+sol.t;trace["saved_states"]=reduce(hcat,states)
    trace["retcode_utf8"]=collect(codeunits(string(sol.retcode)))
    trace["accepted_mass_reference"]=copy(weights.references_used)
    trace["accepted_mass"]=copy(weights.accepted_masses)
    trace["solver_counts"]=[sol.stats.naccept,sol.stats.nreject,sol.stats.nf,sol.stats.njacs,sol.stats.nnonlinconvfail]
    engine_emit!(adapter.diagnostics,tracename,trace)
    if !SciMLBase.successful_retcode(sol)
        engine_emit!(adapter.diagnostics,"failed-state",Dict("state"=>engine_expand!(coordinates,last(sol.u)).*scale,
            "time"=>[start+last(sol.t)],"interval"=>[start,stop]))
        error("full engine interval $index failed: $(sol.retcode)")
    end
    weights.accepted_times==sol.t[2:end] || error("accepted-state weight audit is incomplete")
    refs=vcat(initial_reference,weights.references_used)
    masses=[sum(u[1:ns]) for u in states]
    all(refs.<=masses) || error("nonconservative species mass error weights")
    maximum(abs,diff([u[ns+1] for u in states]))<=20+1e-8 || error("20 K limit failed")
    maximum(diff(sol.t))<=1/(360*50)*(1+1e-12) || error("one-degree limit failed")
    last(sol.t)==stop-start || error("source switch endpoint omitted")
    # Compute independent history checks using shared physical observables.
    output=ic_engine_observables((;model,gas,times=sol.t,states))
    minimum(output["Y"])>=-1e-13 || error("accepted species positivity gate failed")
    maximum(abs,output["volume"].-engine_volume.(start.+sol.t))<1e-10 || error("accepted piston volume gate failed")
    fuel_y=copy(signed.f.network_rhs.state_vector[3].mass_fractions)
    inlet_y=copy(signed.f.network_rhs.state_vector[2].mass_fractions)
    atom_rates=zeros(length(sol.t),length(gas.elements))
    atoms=zeros(size(atom_rates))
    for i in eachindex(sol.t)
        flow=output["mdot_in"][i].*inlet_y.+output["mdot_fuel"][i].*fuel_y.-
            output["mdot_out"][i].*view(output["Y"],i,:)
        atom_rates[i,:]=gas.ele_matrix*(flow./gas.MW)
        atoms[i,:]=gas.ele_matrix*(states[i][1:ns]./gas.MW)
    end
    atom_integrals=[engine_quadratic_integral(view(atom_rates,:,e),sol.t) for e in axes(atoms,2)]
    coarse=unique!(vcat(collect(1:2:length(sol.t)),length(sol.t)))
    atom_coarse=[engine_quadratic_integral(atom_rates[coarse,e],sol.t[coarse]) for e in axes(atoms,2)]
    atom_reference=max.(maximum(abs,atoms;dims=1)[:],abs.(atom_integrals))
    # Use the smallest actual component mass atol applied in this interval;
    # this is no looser than any accepted state's species absolute budget.
    component_mass_atols=fill(minimum(refs)*1e-16,ns)
    element_check=engine_element_budget(gas.ele_matrix,gas.MW,component_mass_atols,
        atom_reference,atoms[end,:].-atoms[1,:].-atom_integrals,atom_coarse.-atom_integrals)
    atom_scale=element_check.total_budget./element_check.relative_limit
    element_error=maximum(element_check.balance_ratio)*element_check.relative_limit
    element_quadrature_error=maximum(element_check.quadrature_ratio)*element_check.relative_limit
    # Preserve evidence before applying a post-solve validation gate.
    diagnostic=Dict{String,Any}("time"=>sol.t,"global_time"=>start.+sol.t,
        "states"=>reduce(hcat,states),"mass_reference"=>refs,"mass"=>masses,
        "atom_rates"=>atom_rates,"atoms"=>atoms,"atom_integrals"=>atom_integrals,
        "atom_coarse_integrals"=>atom_coarse,"atom_scale"=>atom_scale,
        "atom_reference_scale"=>atom_reference,"atom_absolute_budget"=>element_check.absolute_budget,
        "atom_total_budget"=>element_check.total_budget,"component_mass_atols"=>component_mass_atols,
        "element_error"=>[element_error],"element_quadrature_error"=>[element_quadrature_error],
        "inlet_Y"=>inlet_y,"fuel_Y"=>fuel_y,"molecular_weights"=>gas.MW,
        "element_matrix"=>gas.ele_matrix,"interval"=>[start,stop])
    for (name,value) in output;diagnostic["output_"*name]=value;end
    engine_emit!(adapter.diagnostics,"interval"*string(index),diagnostic)
    element_check.balance_pass && element_check.quadrature_pass || error("element flux integration did not converge")
    key="segment_"*lpad(index,2,'0')
    adapter.saved[key*"_time"]=start.+sol.t
    adapter.saved[key*"_state"]=reduce(hcat,states)
    adapter.saved[key*"_mass_reference"]=refs
    adapter.saved[key*"_mass"]=masses
    adapter.saved[key*"_mass_balance"]=output["mass_balance"]
    adapter.saved[key*"_energy_balance"]=output["energy_balance"]
    adapter.saved[key*"_pressure_work"]=output["integrated_work"]
    adapter.saved[key*"_internal_energy"]=output["internal_energy"]
    stats=sol.stats
    record=Dict{String,Any}("index"=>index,"start"=>start,"stop"=>stop,
        "integrated_coordinate_count"=>length(parts.u0),"full_state_count"=>length(scale),
        "rms_norm_denominator"=>length(scale),
        "absent_elements"=>(coordinates===nothing ? String[] : coordinates.absent_elements),
        "removed_species"=>(coordinates===nothing ? String[] : gas.species_names[coordinates.removed]),
        "dynamic_mass_reference"=>dynamic,"maximum_reference_to_mass_ratio"=>maximum(refs./masses),
        "minimum_reference_mass_kg"=>minimum(refs),"maximum_reference_mass_kg"=>maximum(refs),
        "mass_guard_rejections"=>weights.mass_guard_rejections,
        "element_flux_relative_error"=>element_error,"element_flux_quadrature_change"=>element_quadrature_error,
        "element_reference_inventory_kmol"=>atom_reference,
        "element_absolute_budget_kmol"=>element_check.absolute_budget,
        "element_total_budget_kmol"=>element_check.total_budget,
        "minimum_component_mass_atol_kg"=>first(component_mass_atols),
        "accepted_steps"=>stats.naccept,"rejected_steps"=>stats.nreject,
        "rhs_evaluations"=>stats.nf,"jacobian_evaluations"=>stats.njacs,
        "nonlinear_iterations"=>stats.nnonliniter,"nonlinear_failures"=>stats.nnonlinconvfail,
        "seconds_including_first_specialization"=>timed.time)
    push!(adapter.records,record)
    for (name,value) in output
        adapter.saved[key*"_output_"*name]=value
    end
    engine_finish_interval!(adapter.coordinate_factory.audit,index,states,record)
    return (;t=sol.t,u=states,stats=sol.stats,retcode=sol.retcode)
end
