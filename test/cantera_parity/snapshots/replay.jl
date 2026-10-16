#!/usr/bin/env julia
## Replay each captured Cantera transition on the same native grid/state.

## Usage: julia --project=WORKTREE replay.jl spec.json capture_dir report.json
## Runs native solve! from each captured input stage with grid refinement disabled,
## then compares the result with the corresponding Cantera output snapshot.

using Arrhenius, YAML, NPZ, TOML, LinearAlgebra
BLAS.set_num_threads(1)
include("jacobian_trace.jl")

function getarr(npz, key)
    haskey(npz, key) || error("snapshot missing $key")
    return vec(npz[key])
end
function build_state(spec, snap, gas)
    kind=spec["kind"]; names=gas.species_names; z=getarr(snap,"flame/grid")
    T=getarr(snap,"flame/T"); Ys=reduce(hcat,[getarr(snap,"flame/"*s) for s in names])
    if kind in ("free","burner")
        f = kind == "free" ? FreeFlame(gas;T=spec["T"],P=spec["P_atm"]*one_atm,
            X=spec["X"],width=z[end]-z[1],grid=z) : BurnerFlame(gas;
            mdot=spec["mdot"],T=spec["T"],P=spec["P_atm"]*one_atm,
            X=spec["X"],width=z[end]-z[1],grid=z)
        f.state[1,:] .= T ./ 1000; f.state[2:end-1,:] .= permutedims(Ys)
        if kind == "free"
            f.state[end,:] .= getarr(snap,"flame/density")
            f.anchor = clamp(argmin(abs.(f.state[1,:] .- f.fixed_temperature/1000)),2,length(z)-1)
            f.fixed_temperature = 1000*f.state[1,f.anchor]

        else
            f.state[end,:] .= spec["mdot"]
        end
        return f
    elseif kind == "diffusion"
        f=CounterflowDiffusionFlame(gas;fuel=spec["fuel"],oxidizer=spec["oxidizer"],
            mdot_fuel=spec["mdot_fuel"],mdot_oxidizer=spec["mdot_oxidizer"],
            T_fuel=spec["T"],T_oxidizer=get(spec,"T_oxidizer",spec["T"]),
            P=spec["P_atm"]*one_atm,width=z[end]-z[1],grid=z)
        u=getarr(snap,"flame/velocity"); V=getarr(snap,"flame/spreadRate")
        Lam=getarr(snap,"flame/Lambda"); rho=getarr(snap,"flame/density")
        f.state[1,:] .= T ./ 1000; f.state[2:end-3,:] .= permutedims(Ys)
        f.state[end-2,:] .= V ./ Arrhenius._counterflow_Vscale
        f.state[end-1,:] .= Lam ./ Arrhenius._counterflow_Lscale
        f.state[end,:] .= rho .* u; f.grid=z
        return f
    end
    error("unsupported flame kind $kind")
end
function state_matrix(f)
    return f.state
end
function total_stat(stage, key)
    stats=get(stage,"solver_stats",Dict())
    return sum(Int.(get(stats,key,Int[])); init=0)
end
function main(specfile, dir, reportfile, max_transitions=typemax(Int))
    spec=YAML.load_file(specfile); manifest=YAML.load_file(joinpath(dir,"manifest.json"))
    stages=manifest["stages"]; gas=CreateSolution(spec["mechanism"]*".yaml")
    rows=Any[]; firstdiv=nothing
    for i in 1:min(length(stages)-1,max_transitions)
        a=stages[i]; b=stages[i+1]
        input=npzread(joinpath(dir,a["file"])); target=npzread(joinpath(dir,b["file"]))
        f=build_state(spec,input,gas); tic=time()
        err=nothing; dt=get(b,"dt",nothing); refine_changed=false; refine_grid_equal=true
        attempted_dt=dt; decision_note=""
        ct_first_damping=NaN; native_first_damping=NaN
        # For the first counterflow transient, compare the very first Cantera
        # attempt (including a rejected attempt), not only its eventual
        # accepted timestep. Both solvers start from the identical snapshot.
        ct_attempt=nothing
        tracepath=joinpath(dir,"ct_decisions.log")
        if b["event"]=="time_step" && a["event"]=="initial_guess" && isfile(tracepath)
            lines=readlines(tracepath)
            ct_rows=cantera_attempt_trace(tracepath)
            isempty(ct_rows) && error("Cantera Newton log has no transient attempts")
            for (j,line) in enumerate(lines)
                m=match(r"Timestep \(0\) dt=\s*([0-9.eE+-]+)",line)
                if m !== nothing
                    accepted=false; ct_damping=Float64[]
                    for k in j+1:length(lines)
                        cols=split(strip(lines[k]))
                        if length(cols)>=8 && tryparse(Int,cols[1])!==nothing
                            try push!(ct_damping,parse(Float64,cols[2])) catch end
                        end
                        if occursin("Timestep (0) failed",lines[k]) || occursin("Timestep (0) succeeded",lines[k])
                            accepted=occursin("succeeded",lines[k]); break
                        end
                    end
                    ct_attempt=(dt=parse(Float64,m.captures[1]),accepted=accepted,damping=ct_damping,
                                jacobians=ct_rows[1]["jacobians"])
                    attempted_dt=ct_attempt.dt; decision_note="first Cantera dt attempt"
                    break
                end
            end
        end
        iterations=Ref(0); steady_iterations=Ref(0); steady_ok=nothing
        native_jac_before=0
        native_attempt_jacobians=0
        native_log_path=joinpath(dir,"native_stage_$(i).log")
        native_log=open(native_log_path,"w"); old_stdout=stdout; redirect_stdout(native_log)
        try
            ztarget=getarr(target,"flame/grid")
            if length(f.grid)!=length(ztarget) || maximum(abs.(f.grid.-ztarget))>1e-12
                if f isa CounterflowDiffusionFlame
                    refine_changed=Arrhenius._refine_flame!(f;ratio=spec["ratio"],
                        slope=spec["slope"],curve=spec["curve"],prune=get(spec,"prune",0.0),max_points=1200)
                else
                    refine_changed=Arrhenius._refine_flame!(f;ratio=spec["ratio"],
                        slope=spec["slope"],curve=spec["curve"],prune=get(spec,"prune",0.0),max_points=1000)
                end
                refine_grid_equal=length(f.grid)==length(ztarget) && maximum(abs.(f.grid.-ztarget))<=1e-12
                refine_grid_equal || error("native refinement grid differs from Cantera")
            end
            if f isa CounterflowDiffusionFlame
                w=CounterflowWorkspace(f)
                previous=copy(f.state)
                if b["event"]=="time_step" && a["event"]=="initial_guess"
                    steady_input=copy(f.state)
                    steady_ok=Arrhenius._counterflow_newton!(f,w;loglevel=3,iterations=steady_iterations)
                    println("STAGE steady success=",steady_ok," iterations=",steady_iterations[])
                    steady_ok || copyto!(f.state,steady_input)
                end
                native_jac_before=solver_stats(f).jacobians
                println("STAGE transient attempt dt=",attempted_dt)
                ok=attempted_dt===nothing ? Arrhenius._counterflow_newton!(f,w;loglevel=0,iterations) :
                    Arrhenius._counterflow_newton!(f,w;previous,dt=Float64(attempted_dt),maxiters=20,loglevel=3,iterations)
            else
                w=Arrhenius.FlameWorkspace(f)
                previous=copy(f.state)
                if b["event"]=="time_step" && a["event"]=="initial_guess"
                    steady_input=copy(f.state)
                    steady_ok=Arrhenius._flame_newton!(f,w;loglevel=3,iterations=steady_iterations)
                    println("STAGE steady success=",steady_ok," iterations=",steady_iterations[])
                    steady_ok || copyto!(f.state,steady_input)
                end
                native_jac_before=solver_stats(f).jacobians
                println("STAGE transient attempt dt=",attempted_dt)
                ok=attempted_dt===nothing ? Arrhenius._flame_newton!(f,w;loglevel=0,iterations) :
                    Arrhenius._flame_newton!(f,w;previous,dt=Float64(attempted_dt),maxiters=20,loglevel=3,iterations)
            end
            ok || error("native stage Newton solve failed")
        catch e
            err=sprint(showerror,e)
        end
        redirect_stdout(old_stdout); close(native_log)
        elapsed=time()-tic
        native_path=joinpath(dir,"native_stage_$(i).npz")
        native=Dict{String,Any}("grid"=>copy(f.grid),"T"=>temperature(f),
            "velocity"=>velocity(f),"Y"=>permutedims(mass_fractions(f)))
        if f isa CounterflowDiffusionFlame
            native["spreadRate"]=spread_rate(f)
            native["Lambda"]=pressure_curvature(f)
        end
        npzwrite(native_path,native)
        ytarget=reduce(hcat,[getarr(target,"flame/"*s) for s in gas.species_names])
        ttarget=getarr(target,"flame/T")
        dy=length(f.grid)==length(ttarget) ? maximum(abs.(mass_fractions(f).-ytarget')./max.(abs.(ytarget'),1e-300)) : Inf
        dT=length(f.grid)==length(ttarget) ? maximum(abs.(temperature(f).-ttarget)./max.(abs.(ttarget),1.0)) : Inf
        dU=Inf; dV=0.0; dL=0.0
        if length(f.grid)==length(ttarget)
            uct=getarr(target,"flame/velocity")
            dU=maximum(abs.(velocity(f).-uct)./max.(abs.(uct),1e-300))
            if f isa CounterflowDiffusionFlame
                vct=getarr(target,"flame/spreadRate")
                lct=getarr(target,"flame/Lambda")
                dV=maximum(abs.(Arrhenius._counterflow_Vscale .* f.state[end-2,:] .- vct)./max.(abs.(vct),1e-300))
                dL=maximum(abs.(Arrhenius._counterflow_Lscale .* f.state[end-1,:] .- lct)./max.(abs.(lct),1e-300))
            end
        end
        input_grid=getarr(input,"flame/grid"); target_grid=getarr(target,"flame/grid")
        gridok=length(f.grid)==length(target_grid) && maximum(abs.(f.grid.-target_grid)) <= 1e-12
        relative=max(dy,dT,dU,dV,dL)
        same_grid=refine_grid_equal
        ct_steps=max(0,total_stat(b,"steps")-total_stat(a,"steps"))
        ct_jacobians=max(0,total_stat(b,"jacobian_evals")-total_stat(a,"jacobian_evals"))
        native_jacobians=solver_stats(f).jacobians
        native_steps=attempted_dt===nothing || err!==nothing ? 0 : 1
        decisions_equal=ct_steps==native_steps && ct_jacobians==native_jacobians
        decision_first=""
        if ct_attempt!==nothing
            cantera_accept=ct_attempt.accepted
            native_accept=err===nothing
            decisions_equal &= cantera_accept==native_accept
            # Compare applied Newton damping from the attempt-specific log.
            loglines=readlines(native_log_path)
            marker=findfirst(x->occursin("STAGE transient attempt",x),loglines)
            ndamp=Float64[]
            if marker!==nothing
                for line in loglines[marker+1:end]
                    occursin("Newton",line) || continue
                    m=match(r"damping=([0-9.eE+-]+)",line)
                    m===nothing || push!(ndamp,parse(Float64,m.captures[1]))
                end
            end
            cdamp=ct_attempt.damping
            !isempty(cdamp) && (ct_first_damping=first(cdamp))
            !isempty(ndamp) && (native_first_damping=first(ndamp))
            native_attempt_jacobians=max(0,solver_stats(f).jacobians-native_jac_before)
            if native_attempt_jacobians != ct_attempt.jacobians && isempty(decision_first)
                decision_first="Jacobian refresh count: Cantera=$(ct_attempt.jacobians), native=$(native_attempt_jacobians), dt=$(attempted_dt)"
                decisions_equal=false
            end
            if !isempty(cdamp) && !isempty(ndamp) && abs(first(cdamp)-first(ndamp))>1e-6
                isempty(decision_first) && (decision_first="damping: Cantera=$(first(cdamp)), native=$(first(ndamp)), dt=$(attempted_dt)")
                decisions_equal=false
            end
            if cantera_accept != native_accept
                isempty(decision_first) && (decision_first="acceptance: Cantera=$(cantera_accept), native=$(native_accept), dt=$(attempted_dt)")
            end
        end
        ok=err===nothing && gridok && same_grid && relative<=1e-8 && decisions_equal
        if !ok && firstdiv===nothing
            firstdiv=i+1
        end
        push!(rows, Dict("input"=>i,"output"=>i+1,"input_event"=>a["event"],"output_event"=>b["event"],
            "native_seconds"=>elapsed,"relative_state_error"=>relative,"species_error"=>dy,
            "temperature_error"=>dT,"velocity_error"=>dU,"spread_rate_error"=>dV,
            "lambda_error"=>dL,"grid_equal"=>gridok,"input_grid_equal"=>same_grid,
            "native_state"=>native_path,"refined"=>refine_changed,
            "native_log"=>native_log_path,"native_steady_success"=>something(steady_ok,false),
            "native_steady_iterations"=>steady_iterations[],"native_iterations"=>iterations[],
            "native_jacobians"=>native_jacobians,"cantera_jacobians"=>ct_jacobians,
            "native_steps"=>native_steps,"cantera_steps"=>ct_steps,
            "decisions_equal"=>decisions_equal,
            "attempt_dt"=>something(attempted_dt,0.0),"decision_first_difference"=>decision_first,
            "cantera_attempt_accepted"=>(ct_attempt===nothing ? "" : string(ct_attempt.accepted)),
            "cantera_attempt_jacobians"=>(ct_attempt===nothing ? 0 : ct_attempt.jacobians),
            "native_attempt_jacobians"=>native_attempt_jacobians,
            "cantera_first_damping"=>ct_first_damping,
            "native_first_damping"=>native_first_damping,
            "decision_note"=>decision_note,
            "within_tolerance"=>ok,"error"=>something(err,"")))
        !ok && break
    end
    open(reportfile,"w") do io
        fields=("input","output","input_event","output_event","native_seconds",
            "relative_state_error","species_error","temperature_error","velocity_error",
            "spread_rate_error","lambda_error","grid_equal","input_grid_equal",
            "native_state","refined","native_iterations","native_jacobians",
            "cantera_jacobians","native_steps","cantera_steps","decisions_equal",
            "attempt_dt","decision_first_difference","cantera_attempt_accepted",
            "cantera_attempt_jacobians","native_attempt_jacobians","decision_note",
            "cantera_first_damping","native_first_damping",
            "within_tolerance","error")
        println(io,join(fields,","))
        for r in rows
            println(io,join(("\""*replace(string(r[k]),"\""=>"\"\"")*"\"" for k in fields),","))
        end
        println(io,"first_divergence,",something(firstdiv,0),",tolerance,1e-8")
    end






    println("replayed transitions=$(length(rows)) first_divergence=$firstdiv report=$reportfile")
end
main(ARGS[1],ARGS[2],ARGS[3],length(ARGS)>=4 ? parse(Int,ARGS[4]) : typemax(Int))
