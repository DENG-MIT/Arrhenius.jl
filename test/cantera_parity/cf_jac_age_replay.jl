#!/usr/bin/env julia
# Replay Cantera's rejected and accepted counterflow attempts from a captured
# initial snapshot, retaining the native Jacobian workspace between attempts.
using Arrhenius, NPZ, YAML, LinearAlgebra

getvec(npz,key)=vec(npz[key])
function build_diffusion(spec,snapshot,gas)
    z=getvec(snapshot,"flame/grid")
    names=gas.species_names
    Ys=reduce(hcat,[getvec(snapshot,"flame/"*name) for name in names])
    f=CounterflowDiffusionFlame(gas;fuel=spec["fuel"],oxidizer=spec["oxidizer"],
        mdot_fuel=spec["mdot_fuel"],mdot_oxidizer=spec["mdot_oxidizer"],
        T_fuel=spec["T"],T_oxidizer=get(spec,"T_oxidizer",spec["T"]),
        P=spec["P_atm"]*one_atm,width=z[end]-z[1],grid=z)
    f.state[1,:].=getvec(snapshot,"flame/T")./1000
    f.state[2:end-3,:].=permutedims(Ys)
    f.state[end-2,:].=getvec(snapshot,"flame/spreadRate")./Arrhenius._counterflow_Vscale
    f.state[end-1,:].=getvec(snapshot,"flame/Lambda")./Arrhenius._counterflow_Lscale
    f.state[end,:].=getvec(snapshot,"flame/density").*getvec(snapshot,"flame/velocity")
    return f
end
function first_damping(log)
    m=match(r"Counterflow Newton [0-9]+ .*?damping=([0-9.eE+-]+)",log)
    return m===nothing ? nothing : parse(Float64,m.captures[1])
end
function cantera_attempt_trace(path)
    rows=Dict{String,Any}[]
    current=nothing; last_jac=0; pending_refreshes=0
    for line in eachline(path)
        if occursin("Timestep (0) dt=",line)
            dt=first(split(strip(split(line,"dt=")[2])))
            current=Dict{String,Any}("dt"=>parse(Float64,dt),"jacobians"=>0,
                "jacobian_age_start"=>nothing)
            push!(rows,current)
            continue
        end
        if current!==nothing && (occursin("Re-evaluating Jacobian",line) ||
                occursin("Maximum Jacobian age reached",line))
            current["jacobians"]+=1
            pending_refreshes+=1
            # A max-age refresh before the first Newton row makes entry age zero.
            current["jacobian_age_start"]===nothing &&
                occursin("Maximum Jacobian age reached",line) && (current["jacobian_age_start"]=0)
        end
        cols=split(strip(line))
        if length(cols)==8 && tryparse(Int,cols[1])!==nothing &&
                tryparse(Float64,cols[2])!==nothing && occursin("/",cols[end])
            njac=tryparse(Int,cols[7]); age=tryparse(Int,first(split(cols[8],'/')))
            if njac!==nothing && age!==nothing
                delta=max(0,njac-last_jac)
                if current!==nothing
                    newly_observed=max(0,delta-pending_refreshes)
                    current["jacobians"]+=newly_observed
                    if current["jacobian_age_start"]===nothing
                        current["jacobian_age_start"]=newly_observed>0 ? 0 : age-1
                    end
                end
                pending_refreshes=max(0,pending_refreshes-delta)
                last_jac=njac
            end
        end
    end
    return rows
end
function native_entry_age(log)
    m=match(r"Counterflow Jacobian decision: (?:fresh|reuse) age=([0-9]+)",log)
    return m===nothing ? nothing : parse(Int,m.captures[1])
end
function first_difference(ct, native)
    ct_damp=isempty(ct["damping"]) ? nothing : Float64(first(ct["damping"]))
    nd=native["damping"]
    damping_ok = ct_damp===nothing || nd===nothing ? ct_damp===nd :
        abs(nd-ct_damp)<=0.5*10.0^(floor(log10(abs(ct_damp)))-3)
    age_ok=ct["jacobian_age_start"]==native["jacobian_age_start"]
    return ct["accepted"]==native["accepted"] && ct["jacobians"]==native["jacobians"] && age_ok && damping_ok,
        damping_ok,age_ok
end
function main(specpath,capture,outpath)
    spec=YAML.load_file(specpath)
    manifest=YAML.load_file(joinpath(capture,"manifest.json"))
    trace=YAML.load_file(joinpath(capture,"ct_decisions.json"))
    attempts=trace["attempts"]
    ct_age_trace=cantera_attempt_trace(joinpath(capture,"ct_decisions.log"))
    length(ct_age_trace)>=length(attempts) || error("Cantera log has fewer timestep entries than JSON attempts")
    gas=CreateSolution(spec["mechanism"]*".yaml")
    snap=npzread(joinpath(capture,manifest["stages"][1]["file"]))
    f=build_diffusion(spec,snap,gas); initial=copy(f.state)
    w=CounterflowWorkspace(f); rows=Any[]
    steady_log=open(joinpath(capture,"native_steady.log"),"w")
    steady_ok=redirect_stdout(steady_log) do
        Arrhenius._counterflow_newton!(f,w;loglevel=3)
    end
    close(steady_log); copyto!(f.state,initial)
    for (i,ct) in enumerate(attempts)
        previous=copy(f.state); logpath=joinpath(capture,"native_attempt_$(lpad(i,2,'0')).log")
        log=open(logpath,"w"); iterations=Ref(0); before=w.jac_evals
        accepted=redirect_stdout(log) do
            Arrhenius._counterflow_newton!(f,w;previous,dt=Float64(ct["dt"]),
                maxiters=200,loglevel=3,iterations)
        end
        close(log)
        native_log=read(logpath,String)
        native=Dict("accepted"=>accepted,"jacobians"=>w.jac_evals-before,
            "jacobian_age_start"=>native_entry_age(native_log),
            "damping"=>first_damping(native_log),"iterations"=>iterations[])
        # ct_decisions.json stores the per-solve Jacobian delta. Log markers
        # are cumulative across retries and may be emitted before the next
        # attempt consumes a refreshed Jacobian.
        logrow=ct_age_trace[i]
        ctrow=Dict("dt"=>ct["dt"],"accepted"=>ct["accepted"],
            "jacobians"=>logrow["jacobians"],
            "jacobian_age_start"=>logrow["jacobian_age_start"],"damping"=>ct["damping"])
        matched,damping_ok,age_ok=first_difference(ctrow,native)
        push!(rows,Dict("attempt"=>i,"cantera"=>ctrow,"native"=>native,
            "jacobian_age_matches"=>age_ok,"damping_matches_print_precision"=>damping_ok,
            "decisions_match"=>matched))
        if accepted
            # An accepted attempt advances the state and clears the retry count.
        else
            copyto!(f.state,previous)
            get(ct,"reset_after",false) && Arrhenius._counterflow_reset_species!(f.state,gas.n_species)
        end
        matched || break
    end
    report=Dict("case"=>get(manifest,"case","unknown"),"initial_steady_failed"=>!steady_ok,
        "replayed_attempts"=>length(rows),"matched_attempts"=>count(r->r["decisions_match"],rows),
        "first_difference"=>(any(!r["decisions_match"] for r in rows) ? length(rows) : nothing),
        "attempts"=>rows)
    open(outpath,"w") do io
        YAML.write(io,report)
    end
    println(YAML.write(report))
end
length(ARGS)==3 || error("usage: cf_jac_age_replay.jl spec.json capture_dir report.yaml")
main(ARGS...)
