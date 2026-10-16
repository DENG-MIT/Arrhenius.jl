# Full calculation from Cantera's non_ideal_shock_tube.py, without plots or I/O.
# Usage: julia --project=SOLVER_ENV real_gas_timing.jl MECHANISM_DIR OUTPUT.npz
#        [REPETITIONS=7] [informational|controlled|validate-only] [--smoke] [--finite-difference] [--qndf]
# SOLVER_ENV needs SciMLBase, OrdinaryDiffEqSDIRK and ForwardDiff; prepare MECHANISM_DIR with
# real_gas_reactor_cases.py. --qndf additionally needs OrdinaryDiffEqBDF.
# Run real_gas_timing.py afterward for comparison.
for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    ENV[name] = "1"
end
const SHOCKTUBE_EXAMPLE_PATH = joinpath(@__DIR__,"..","example","reactors","non_ideal_shock_tube.jl")
const SHOCKTUBE_AD_PATH = joinpath(@__DIR__,"..","example","reactors","real_gas_ad_jacobian.jl")
import_seconds = @elapsed begin
    @eval using Arrhenius, NPZ, SHA, Dates, TOML, LinearAlgebra, Libdl
    include(SHOCKTUBE_EXAMPLE_PATH)
    "--qndf" in ARGS && include(joinpath(dirname(SHOCKTUBE_EXAMPLE_PATH),"real_gas_qndf_solver.jl"))
end

function shocktube_thread_settings(;enforce=true)
    Threads.nthreads()==1 || error("launch Julia with JULIA_NUM_THREADS=1")
    enforce && BLAS.set_num_threads(1)
    result = Dict("julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads())
    result["blas_threads"]==1 || error("Julia BLAS must use one thread")
    if Sys.isapple()
        # Apple thread_api.h: SINGLE_THREADED=1. This is local to the calling
        # thread, which runs the sequential ODE calculations below.
        if enforce
            status=ccall((:BLASSetThreading,"/System/Library/Frameworks/Accelerate.framework/Accelerate"),Cint,(Cuint,),1)
            status==0 || error("Accelerate single-thread setting failed")
        end
        result["accelerate_threading_mode"]=Int(ccall((:BLASGetThreading,"/System/Library/Frameworks/Accelerate.framework/Accelerate"),Cuint,()))
        result["accelerate_threading_mode"]==1 || error("Accelerate must use one thread")
    end
    # LinearSolve may select MKL independently of Julia's default BLAS.
    for path in unique(filter(path -> occursin("mkl_rt",lowercase(basename(path))),Libdl.dllist()))
        handle=Libdl.dlopen(path)
        try
            setter=Libdl.dlsym(handle,:MKL_Set_Num_Threads)
            getter=Libdl.dlsym(handle,:MKL_Get_Max_Threads)
            enforce && ccall(setter,Cvoid,(Cint,),1)
            count=Int(ccall(getter,Cint,()))
            result["mkl_threads"]=count
            count==1 || error("LinearSolve MKL must use one thread")
        finally
            Libdl.dlclose(handle)
        end
    end
    return result
end

function shocktube_require_replay(actual,expected,repetition)
    isequal(actual,expected) || error("native warm repetition $repetition differs from the checked first trajectory")
    return true
end

function shocktube_source_hashes(root)
    paths=[joinpath(folder,file) for (folder,_,files) in walkdir(joinpath(root,"src")) for file in files]
    append!(paths,[joinpath(root,"example","reactors",file) for file in
        ("non_ideal_shock_tube.jl","real_gas_ode_solver.jl","real_gas_ad_jacobian.jl",
         "real_gas_trial_states.jl","real_gas_qndf_solver.jl")])
    push!(paths,joinpath(root,"validation","real_gas_timing.jl"))
    Dict(replace(relpath(path,root),'\\'=>'/')=>bytes2hex(sha256(read(path))) for path in paths)
end

shocktube_mechanism_hashes(directory)=Dict(file=>bytes2hex(sha256(read(joinpath(directory,file))))
    for file in ("dodecane_RK.yaml","dodecane_IG.yaml","dodecane_IG.yaml.npz"))

function shocktube_timing_main(args)
    length(args)>=2 || error("supply MECHANISM_DIR OUTPUT.npz [REPETITIONS] [QUALIFICATION] [--smoke] [--finite-difference] [--qndf]")
    directory,output = args[1:2]
    positional = filter(x -> x ∉ ("--smoke","--ad","--finite-difference","--qndf"),args[3:end])
    repetitions = length(positional)>=1 ? parse(Int,positional[1]) : 7
    qualification = length(positional)>=2 ? positional[2] : "informational"
    qualification in ("informational","controlled","validate-only") || error("invalid qualification")
    repetitions>=7 || error("at least seven warm repetitions required")
    smoke = "--smoke" in args
    ad = !("--finite-difference" in args)
    qndf = "--qndf" in args
    qndf && !ad && error("--qndf requires the AD Jacobian")
    qndf && !isdefined(@__MODULE__,:shocktube_qndf_integrator) &&
        error("include real_gas_qndf_solver.jl before calling this driver with --qndf")
    solver = qndf ? :qndf : :sdirk
    root = normpath(joinpath(@__DIR__,".."))
    realpath(pathof(Arrhenius))==realpath(joinpath(root,"src","Arrhenius.jl")) ||
        error("the benchmark and Arrhenius package must use the same checkout")
    hashes_before=shocktube_source_hashes(root)
    mechanisms_before=shocktube_mechanism_hashes(directory)
    gas = CreateSolution(joinpath(directory,"dodecane_IG.yaml"))
    model = RedlichKwongThermo(joinpath(directory,"dodecane_RK.yaml"))
    trial_policy = qndf ? SignedIntegerShockTubeTrials(gas,joinpath(directory,"dodecane_IG.yaml")) : ClippedShockTubeTrials()
    thread_checks = Dict("before_first"=>shocktube_thread_settings())
    started = time_ns()
    first_result = shocktube_calculations(gas,model;smoke,jacobian=ad ? :ad : :finite_difference,solver,trial_policy)
    first_seconds = (time_ns()-started)/1e9
    # Check any backend loaded lazily during the first/JIT invocation.
    thread_checks["before_warm"] = shocktube_thread_settings()
    println("First calculation complete: ",length(first_result)," trajectories; ",first_seconds," s including JIT"); flush(stdout)
    warm_seconds = Float64[]
    warm_delays = Vector{Float64}[]
    warm_matches_first = Bool[]
    if qualification != "validate-only"
        GC.gc()
        for repetition in 1:repetitions
            started = time_ns()
            result = shocktube_calculations(gas,model;smoke,jacobian=ad ? :ad : :finite_difference,solver,trial_policy)
            push!(warm_seconds,(time_ns()-started)/1e9)
            thread_checks["after_warm_"*string(repetition)] = shocktube_thread_settings(;enforce=false)
            push!(warm_delays,[r.delay for r in result])
            # Compare every saved state, time, endpoint and solver count outside
            # the timing interval. Exact replay transfers the first-run checks
            # to every measured repetition without another ODE calculation.
            push!(warm_matches_first,shocktube_require_replay(result,first_result,repetition))
            println("Warm repetition ",repetition,": ",last(warm_seconds)," s"); flush(stdout)
        end
    end
    utf8(value) = collect(codeunits(string(value)))
    cpu = if Sys.islinux()
        strip(split(first(filter(l -> startswith(l,"model name"),readlines("/proc/cpuinfo"))),':';limit=2)[2])
    elseif Sys.isapple()
        readchomp(`sysctl -n machdep.cpu.brand_string`)
    else
        Sys.CPU_NAME
    end
    hashes_after=shocktube_source_hashes(root)
    mechanisms_after=shocktube_mechanism_hashes(directory)
    hashes_before==hashes_after || error("native source inventory or bytes changed during calculations")
    mechanisms_before==mechanisms_after || error("native mechanism bytes changed during calculations")
    toml_string(value) = sprint(io -> TOML.print(io,value))
    manifest = joinpath(dirname(Base.active_project()),"Manifest.toml")
    data = Dict{String,Any}(
        "cpu_utf8"=>utf8(cpu),"julia_version_utf8"=>utf8(VERSION),
        "system_utf8"=>utf8(Sys.isapple() ? "Darwin" : Sys.islinux() ? "Linux" : string(Sys.KERNEL)),
        "kernel_release_utf8"=>utf8(Sys.isunix() ? readchomp(`uname -r`) : "unknown"),
        "platform_utf8"=>utf8(Sys.MACHINE),"qualification_utf8"=>utf8(qualification),
        "timestamp_utc_utf8"=>utf8(Dates.now(Dates.UTC)),"scope_utf8"=>utf8(smoke ? "1000_K_pair" : "full_34_trajectories"),
        "source_hashes_toml_utf8"=>utf8(toml_string(hashes_before)),
        "source_hashes_after_toml_utf8"=>utf8(toml_string(hashes_after)),
        "mechanism_hashes_toml_utf8"=>utf8(toml_string(mechanisms_before)),
        "mechanism_hashes_after_toml_utf8"=>utf8(toml_string(mechanisms_after)),
        "manifest_toml_utf8"=>utf8(read(manifest,String)),
        "thread_environment_toml_utf8"=>utf8(toml_string(Dict(name=>get(ENV,name,"")
            for name in ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")))),
        "thread_checks_toml_utf8"=>utf8(toml_string(thread_checks)),
        "blas_config_utf8"=>utf8(BLAS.get_config()),
        "loaded_libraries_utf8"=>utf8(join(Libdl.dllist(),"\n")),
        "harness_sha256_utf8"=>utf8(bytes2hex(sha256(read(@__FILE__)))),
        "jacobian_helper_sha256_utf8"=>utf8(ad ? bytes2hex(sha256(read(SHOCKTUBE_AD_PATH))) : "unused"),
        "julia_threads"=>[Threads.nthreads()],"blas_threads"=>[BLAS.get_num_threads()],
        "logical_cpus"=>[Sys.CPU_THREADS],"import_seconds"=>[import_seconds],
        "first_seconds"=>[first_seconds],"warm_seconds"=>warm_seconds,
        "warm_matches_first"=>UInt8.(warm_matches_first),
        "warm_delays"=>isempty(warm_delays) ? zeros(length(first_result),0) : reduce(hcat,warm_delays),
        "native_rtol"=>[qndf ? 1e-9 : 1e-13],"native_atol"=>[qndf ? 1e-19 : 1e-26],
        "native_temperature_atol_K"=>[qndf ? 1e-6 : 1e-26],
        "absolute_tolerance_formula_utf8"=>utf8(qndf ?
            "C atol = 1e-19 kmol/m^3; T atol = 1e-6 K; z = [C; T/1000 K], so z_T atol = 1e-9" :
            "scalar 1e-26 tolerance for every physical [Y; T] state component"),
        "trial_policy_utf8"=>utf8(qndf ? "mechanism-bound signed C1/C2/C3" : "clipped concentrations"),
        "solver_utf8"=>utf8(qndf ? "OrdinaryDiffEqBDF.QNDF with scaled concentrations, cached composition AD and nonlinear coefficient 0.01" :
            "OrdinaryDiffEqSDIRK.KenCarp4 with "*(ad ? "cached composition AD" : "finite-difference Jacobian")),
        "case_order_utf8"=>utf8(join([r.phase*"_"*string(r.temperature) for r in first_result],",")),
    )
    for r in first_result
        key = r.phase*"_"*string(r.temperature)
        data[key*"_time"],data[key*"_state"] = r.times,reduce(hcat,r.states)
        data[key*"_delay"],data[key*"_steps"] = [r.delay],[r.steps]
        data[key*"_final_time"],data[key*"_final_state"] = [r.final_time],r.final_state
        data[key*"_density"] = [r.density]
        data[key*"_rhs_evaluations"],data[key*"_jacobian_evaluations"] = [r.rhs_evaluations],[r.jacobian_evaluations]
        for field in (:accepted_steps,:rejected_steps,:linear_solves,:matrix_updates,
            :nonlinear_iterations,:nonlinear_convergence_failures)
            data[key*"_"*string(field)] = [getproperty(r,field)]
        end
    end
    npzwrite(output,data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    shocktube_timing_main(ARGS)
end
