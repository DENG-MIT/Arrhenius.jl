#!/usr/bin/env julia
# """Compare default and Cantera-compatible Newton policy at snapshot and adaptive levels.
#
# Input is TOML with repeated `[[case]]` entries: `id`, `spec`, and `capture_dir`.
# Capture directories use the format produced by `test/cantera_parity/snapshots/capture.py`.
# Run under `run_capped`; snapshot replay runs one fixed-grid transition per policy.
# """
using Arrhenius, YAML, NPZ, TOML, LinearAlgebra
BLAS.set_num_threads(1)
include(joinpath(@__DIR__,"..","test","cantera_parity","snapshots","replay.jl"))

function json_value(x)
    x === nothing && return "null"
    x isa Bool && return x ? "true" : "false"
    x isa Number && return isfinite(x) ? string(x) : "null"
    x isa AbstractString && return "\""*replace(x,"\\"=>"\\\\","\""=>"\\\"","\n"=>"\\n")*"\""
    x isa AbstractDict && return "{"*join((json_value(string(k))*":"*json_value(v) for (k,v) in x),",")*"}"
    if x isa Tuple || x isa AbstractVector
        return "["*join(json_value.(collect(x)),",")*"]"
    end
    error("unsupported JSON value $(typeof(x))")
end

function case_flame(spec, root)
    mechanism=normpath(get(spec,"mechanism_path",joinpath(root,spec["mechanism"]*".yaml")))
    gas=CreateSolution(mechanism)
    P=spec["P_atm"]*one_atm
    kind=spec["kind"]
    if kind=="free"
        f=FreeFlame(gas;T=spec["T"],P,X=spec["X"],width=spec["width"],discretization=:finite_difference)
    elseif kind=="burner"
        f=BurnerFlame(gas;mdot=spec["mdot"],T=spec["T"],P,X=spec["X"],width=spec["width"])
    elseif kind=="diffusion"
        f=CounterflowDiffusionFlame(gas;fuel=spec["fuel"],oxidizer=spec["oxidizer"],
            mdot_fuel=spec["mdot_fuel"],mdot_oxidizer=spec["mdot_oxidizer"],T_fuel=spec["T"],
            T_oxidizer=get(spec,"T_oxidizer",spec["T"]),P,width=spec["width"])
    else
        error("unsupported policy_evidence kind: $kind")
    end
    (f isa FreeFlame || f isa BurnerFlame) &&
        set_transport!(f,Symbol(replace(get(spec,"transport","mixture-averaged"),"-"=>"_")))
    return f
end

function element_profile_drift(gas,Y)
    base=zeros(length(gas.elements))
    maxerr=0.0
    for e in eachindex(gas.elements)
        for k in eachindex(gas.species_names)
            base[e]+=Y[k,1]/gas.MW[k]*gas.ele_matrix[e,k]
        end
        for j in axes(Y,2)
            value=sum(Y[k,j]/gas.MW[k]*gas.ele_matrix[e,k] for k in eachindex(gas.species_names))
            maxerr=max(maxerr,abs(value-base[e])/max(abs(base[e]),1e-300))
        end
    end
    return maxerr
end
element_profile_drift(f) = f isa CounterflowDiffusionFlame ? nothing : element_profile_drift(f.gas,mass_fractions(f))

function element_source_balance(f)
    r=similar(f.state)
    if f isa CounterflowDiffusionFlame
        w=Arrhenius.CounterflowWorkspace(f)
        Arrhenius.counterflow_residual!(r,f,f.state,w)
        source=w.properties.source
    else
        w=Arrhenius.FlameWorkspace(f)
        flame_residual!(r,f,f.state,w)
        source=w.source
    end
    max_relative=0.0
    max_absolute=0.0
    for e in eachindex(f.gas.elements), j in axes(source,2)
        imbalance=0.0; magnitude=0.0
        for k in eachindex(f.gas.species_names)
            contribution=f.gas.ele_matrix[e,k]*source[k,j]
            imbalance+=contribution
            magnitude+=abs(contribution)
        end
        max_absolute=max(max_absolute,abs(imbalance))
        max_relative=max(max_relative,abs(imbalance)/max(magnitude,1e-300))
    end
    return (relative=max_relative,absolute=max_absolute)
end

function interpolate_profile(x,y,xnew)
    out=similar(xnew,Float64)
    for (j,xj) in enumerate(xnew)
        i=clamp(searchsortedlast(x,xj),1,length(x)-1)
        q=clamp((xj-x[i])/(x[i+1]-x[i]),0.0,1.0)
        out[j]=(1-q)*y[i]+q*y[i+1]
    end
    return out
end

function reference_accuracy(f,capture_dir)
    manifest=YAML.load_file(joinpath(capture_dir,"manifest.json"))
    stage=manifest["stages"][end]
    ref=NPZ.npzread(joinpath(capture_dir,stage["file"]))
    zref=getarr(ref,"flame/grid"); z=collect(f.grid)
    Tref=interpolate_profile(zref,getarr(ref,"flame/T"),z)
    terr=maximum(abs.(temperature(f).-Tref)./max.(abs.(Tref),1.0))
    Yref=reduce(hcat,[interpolate_profile(zref,getarr(ref,"flame/"*s),z) for s in f.gas.species_names])
    yerr=maximum(abs.(mass_fractions(f).-Yref')./max.(abs.(Yref'),1e-300))
    uref=interpolate_profile(zref,getarr(ref,"flame/velocity"),z)
    uerr=maximum(abs.(velocity(f).-uref)./max.(abs.(uref),1e-300))
    return Dict("max_relative_temperature_error"=>terr,
        "max_relative_species_error"=>yerr,"max_relative_velocity_error"=>uerr)
end

function adaptive_result(specfile, capture_dir, policy, label, root)
    spec=YAML.load_file(specfile); f=case_flame(spec,root); t0=time(); err=""
    try
        solve!(f;ratio=get(spec,"ratio",3.0),slope=get(spec,"slope",.06),
            curve=get(spec,"curve",.12),prune=get(spec,"prune",0.0),
            max_points=1200,max_time_steps=800,newton_options=policy)
    catch e
        err=sprint(showerror,e)
    end
    stats=solver_stats(f); Y=mass_fractions(f)
    element_balance=element_source_balance(f)
    return Dict("level"=>"adaptive","case"=>label,"policy"=>policy.cantera_compatible ? "cantera" : "default",
        "converged"=>isempty(err) && f.converged,"error"=>err,"seconds"=>time()-t0,
        "grid_size"=>length(f.grid),"max_sumY_error"=>maximum(abs.(sum(Y;dims=1).-1)),
        "max_element_profile_drift"=>element_profile_drift(f),
        "max_element_source_balance_relative_error"=>element_balance.relative,
        "max_element_source_balance_absolute"=>element_balance.absolute,
        "solver_stats"=>Dict("newton_calls"=>stats.newton_calls,
            "newton_iterations"=>stats.newton_iterations,"jacobians"=>stats.jacobians,
            "time_steps"=>stats.time_steps,"time_step_failures"=>stats.time_step_failures,
            "steady_failures"=>stats.steady_failures),
        "reference_accuracy"=>isempty(err) ? reference_accuracy(f,capture_dir) : nothing)
end

function snapshot_result(specfile,capture_dir,policy,label,transition)
    out=mktempdir(); csv=joinpath(out,"replay.csv")
    rows=main(specfile,capture_dir,csv,1,transition,policy)
    values=rows[1]
    native=NPZ.npzread(values["native_state"])
    Y=permutedims(native["Y"])
    spec=YAML.load_file(specfile)
    mechanism=normpath(get(spec,"mechanism_path",joinpath(@__DIR__,"..",spec["mechanism"]*".yaml")))
    gas=CreateSolution(mechanism)
    return Dict("level"=>"snapshot_replay","case"=>label,"transition"=>transition,
        "policy"=>policy.cantera_compatible ? "cantera" : "default",
        "converged"=>values["error"]=="",
        "within_snapshot_tolerance"=>values["within_tolerance"],
        "grid_size"=>length(native["grid"]),"seconds"=>values["native_seconds"],
        "max_sumY_error"=>maximum(abs.(sum(Y;dims=1).-1)),
        "max_element_profile_drift"=>spec["kind"]=="diffusion" ? nothing : element_profile_drift(gas,Y),
        "snapshot_relative_state_error"=>values["relative_state_error"],
        "solver_stats"=>Dict("newton_iterations"=>values["native_iterations"],
            "jacobians"=>values["native_jacobians"],
            "cantera_jacobians"=>values["cantera_jacobians"],
            "native_steps"=>values["native_steps"],
            "cantera_steps"=>values["cantera_steps"]),"error"=>values["error"])
end

function main(casefile,outfile)
    cases=TOML.parsefile(casefile)["case"]
    isempty(cases) && error("case list is empty")
    root=normpath(joinpath(@__DIR__,".."))
    policies=(NewtonPolicyOptions(),Arrhenius._newton_policy(nothing,true))
    results=Any[]
    for c in cases
        spec=normpath(joinpath(dirname(abspath(casefile)),c["spec"]))
        capture=normpath(joinpath(dirname(abspath(casefile)),c["capture_dir"]))
        isfile(spec) || error("missing spec: $spec")
        isfile(joinpath(capture,"manifest.json")) || error("missing snapshot manifest: $capture")
        for policy in policies
            transition=get(c,"transition",1)
            push!(results,snapshot_result(spec,capture,policy,c["id"],transition))
            push!(results,adaptive_result(spec,capture,policy,c["id"],root))
        end
    end
    open(outfile,"w") do io; write(io,json_value(Dict("cases"=>results))); write(io,"\n"); end
    println("policy evidence: ",length(cases)," cases, ",length(results)," records -> ",outfile)
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("usage: julia --project validation/policy_evidence.jl CASES.toml OUTPUT.json")
    main(ARGS[1],ARGS[2])
end
