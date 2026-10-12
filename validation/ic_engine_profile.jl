# One cold invocation and one warm profiled invocation of the shared engine.
# This is diagnostic, never a qualifying timing result. See Julia's manual:
# https://docs.julialang.org/en/v1/manual/performance-tips/
# julia --project=SOLVER_ENV validation/ic_engine_profile.jl MECHANISM.yaml OUTPUT_PREFIX
include("ic_engine_timing.jl")
using Profile, InteractiveUtils

struct EngineProfileIntegrator{F}
    integrator::F
    reports::Vector{Dict{String,Float64}}
end
function (adapter::EngineProfileIntegrator)(problem;kwargs...)
    measurement = @timed adapter.integrator(problem;kwargs...)
    stats = measurement.value.stats
    push!(adapter.reports,Dict("seconds"=>measurement.time,"allocated_bytes"=>measurement.bytes,
        "accepted_steps"=>stats.naccept,"rejected_steps"=>stats.nreject,
        "solver_rhs_evaluations"=>stats.nf,"jacobian_evaluations"=>stats.njacs,
        "linear_solves"=>stats.nsolve,"matrix_updates"=>stats.nw,
        "nonlinear_iterations"=>stats.nnonliniter,"nonlinear_convergence_failures"=>stats.nnonlinconvfail))
    return measurement.value
end

function profiled_engine_calculation(mechanism)
    adapter = EngineProfileIntegrator(native_engine_sdirk,Dict{String,Float64}[])
    result = solve_ic_engine(mechanism;integrator=adapter)
    output = ic_engine_observables(result)
    return (;result,output,summary=ic_engine_summary(result,output),segments=adapter.reports)
end

function hot_engine_callbacks(calculation,path)
    gas = calculation.result.gas
    # Initial mixing, fuel injection, and hot expansion exercise different
    # pressure-flow/chemical branches of the same native callbacks.
    measurements = Dict{String,Any}()
    open(path,"w") do io
        for target in (0.0005,0.0199,0.025)
            index = searchsortedfirst(calculation.result.times,target)
            t = calculation.result.times[index]
            state = copy(calculation.result.states[index])
            model = engine_network(gas;segment_time=t)
            problem = moving_wall_problem(model,(0.,0.01);initial_state=state)
            du,J = zero(state),zeros(length(state),length(state))
            problem.f(du,state,nothing,t)
            problem.jac(J,state,nothing,t)
            rhs_measurement = @timed for _ in 1:1000
                problem.f(du,state,nothing,t)
            end
            jac_measurement = @timed for _ in 1:10
                problem.jac(J,state,nothing,t)
            end
            label = string(t)
            measurements[label] = Dict("rhs_seconds_per_call"=>rhs_measurement.time/1000,
                "rhs_bytes_per_call"=>rhs_measurement.bytes/1000,
                "jacobian_seconds_per_call"=>jac_measurement.time/10,
                "jacobian_bytes_per_call"=>jac_measurement.bytes/10)
            println(io,"Native moving-wall RHS at t = ",t)
            code_warntype(io,problem.f,Tuple{typeof(du),typeof(state),Nothing,Float64};debuginfo=:source)
            println(io,"Native moving-wall Jacobian at t = ",t)
            code_warntype(io,problem.jac,Tuple{typeof(J),typeof(state),Nothing,Float64};debuginfo=:source)
        end
    end
    return measurements
end

function engine_profile_main(args)
    length(args)==2 || error("supply MECHANISM.yaml OUTPUT_PREFIX")
    mechanism,prefix = args
    realpath(pathof(Arrhenius))==realpath(joinpath(@__DIR__,"..","src","Arrhenius.jl")) ||
        error("the profiler and Arrhenius package must use the same checkout")
    thread_checks = Dict("before_first"=>benchmark_julia_thread_settings())
    cold = @timed profiled_engine_calculation(mechanism)
    engine_require_checks(cold.value,"cold profiling invocation")
    thread_checks["before_warm"] = benchmark_julia_thread_settings()
    println("Cold engine including JIT: ",cold.time," s; checks ",engine_timing_checks(cold.value)); flush(stdout)
    GC.gc()
    Profile.init(n=10^7,delay=0.005)
    Profile.clear()
    warm = @timed Profile.@profile profiled_engine_calculation(mechanism)
    thread_checks["after_warm"] = benchmark_julia_thread_settings(;enforce=false)
    println("One warm profiled engine: ",warm.time," s; ",warm.bytes," allocated bytes"); flush(stdout)
    repeated = engine_require_replay(warm.value,cold.value,"warm profiling invocation")
    valid = engine_require_checks(warm.value,"warm profiling invocation") && repeated
    open(prefix*".profile.txt","w") do io
        Profile.print(io;format=:flat,sortedby=:count,mincount=20)
    end
    callbacks = hot_engine_callbacks(warm.value,prefix*".inference.txt")
    totals = Dict(key=>sum(segment[key] for segment in warm.value.segments)
                  for key in keys(warm.value.segments[1]))
    report = Dict("qualification"=>"diagnostic_only","correctness_pass"=>valid,
        "warm_matches_cold"=>repeated,"cold_seconds_including_JIT"=>cold.time,
        "warm_seconds_with_sampling_profiler"=>warm.time,"warm_allocated_bytes"=>warm.bytes,
        "warm_gc_seconds"=>warm.gctime,"solver_totals"=>totals,
        "segments"=>warm.value.segments,"hot_callbacks"=>callbacks,"native_summary"=>warm.value.summary,
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "numerical_thread_checks"=>thread_checks,
        "thread_helper_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"numerical_threads.jl")))),
        "julia_version"=>string(VERSION),"harness_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "source_hashes"=>engine_source_hashes(normpath(joinpath(@__DIR__,".."))))
    open(prefix*".toml","w") do io
        TOML.print(io,report;sorted=true)
    end
    valid || error("warm profiling calculation failed replay or conservation checks")
end

if abspath(PROGRAM_FILE)==@__FILE__
    engine_profile_main(ARGS)
end
