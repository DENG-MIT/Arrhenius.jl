# Export the complete accepted history from the public QNDF calculation.
# julia --project=CALLER_ENV validation/ic_engine_qndf_case.jl dodecane_IG.yaml native.npz
# Caller dependencies: Arrhenius, ForwardDiff, OrdinaryDiffEqBDF, SciMLBase.
using Arrhenius, TOML, SHA, LinearAlgebra, Libdl
import Arrhenius.NPZ
include(joinpath(@__DIR__,"..","example","reactors","ic_engine_qndf_solver.jl"))
include(joinpath(@__DIR__,"numerical_threads.jl"))

engine_case_hashfile(path) = bytes2hex(sha256(read(path)))
function engine_case_inputs(mechanism)
    root = dirname(@__DIR__)
    paths = [joinpath(root,"example","reactors","ic_engine_setup.jl"),
             joinpath(root,"example","reactors","ic_engine_qndf_solver.jl")]
    append!(paths,readdir(joinpath(root,"example","reactors","engine_qndf");join=true))
    sources = Dict(replace(relpath(p,root),'\\'=>'/')=>engine_case_hashfile(p)
        for p in paths if endswith(p,".jl"))
    core = pkgdir(Arrhenius)
    core_hashes = Dict(replace(relpath(joinpath(folder,name),core),'\\'=>'/')=>engine_case_hashfile(joinpath(folder,name))
        for (folder,_,names) in walkdir(joinpath(core,"src")) for name in names if endswith(name,".jl"))
    project = Base.active_project()
    project === nothing && error("run the exporter in the caller's solver project")
    manifest = joinpath(dirname(project),"Manifest.toml")
    Dict("mechanism_sha256"=>engine_case_hashfile(mechanism),
        "sidecar_sha256"=>engine_case_hashfile(mechanism*".npz"),
        "source_sha256"=>core_hashes,"helper_sha256"=>sources,
        "exporter_sha256"=>engine_case_hashfile(@__FILE__),
        "thread_guard_sha256"=>engine_case_hashfile(joinpath(@__DIR__,"numerical_threads.jl")),
        "caller_project_sha256"=>engine_case_hashfile(project),
        "caller_manifest_sha256"=>engine_case_hashfile(manifest))
end

function engine_case_runtime(;enforce=false)
    settings = benchmark_julia_thread_settings(;enforce)
    numerical = filter(path->occursin(r"blas|mkl|accelerate|iomp|libgomp"i,basename(path)),Libdl.dllist())
    Dict("checks_pass"=>true,"settings"=>settings,"blas_config"=>string(BLAS.get_config()),
        "environment"=>Dict(name=>get(ENV,name,"") for name in
            ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")),
        "loaded_numerical_library_sha256"=>Dict(basename(path)=>engine_case_hashfile(path)
            for path in numerical if isfile(path)))
end

function engine_case_archive_data(calculation)
    data = copy(calculation.accepted)
    data["species_names_utf8"] = collect(codeunits(join(calculation.result.gas.species_names,"\n")))
    data["exact_source_initial_state"] = calculation.initial_state
    data["states"] = reduce(hcat,calculation.result.states)
    data["quadrature_terms"] = calculation.checks["quadrature_terms"]
    data["coarse_quadrature_terms"] = calculation.checks["coarse_quadrature_terms"]
    data["ledger_work"] = [calculation.checks["ledger_work"]]
    data["checks_pass"] = [UInt8(calculation.checks["checks_pass"])]
    for (key,value) in calculation.output
        data["output_"*key] = value
    end
    data
end

function write_engine_case(calculation,destination,before,after,runtime_before,runtime_after)
    # Serialize before rejecting a changed input/backend or failed calculation.
    NPZ.npzwrite(destination,engine_case_archive_data(calculation))
    unchanged = before==after
    valid = calculation.checks["checks_pass"] && unchanged && runtime_before["checks_pass"] && runtime_after["checks_pass"]
    record = Dict("checks_pass"=>valid,"checks"=>calculation.checks,
        "summary"=>calculation.summary,"segments"=>calculation.records,
        "mechanism_sha256"=>before["mechanism_sha256"],"sidecar_sha256"=>before["sidecar_sha256"],
        "native_sha256"=>engine_case_hashfile(destination),"julia_version"=>string(VERSION),
        "inputs_unchanged"=>unchanged,"inputs_before"=>before,"inputs_after"=>after,
        "runtime_before"=>runtime_before,"runtime_after"=>runtime_after,
        "stored_trajectory_input_to_calculation"=>false)
    open(replace(destination,r"\.npz$"=>".toml"),"w") do io
        TOML.print(io,record;sorted=true)
    end
    valid || error("engine calculation, input-identity or numerical-backend checks failed; outputs were saved")
    record
end

function export_ic_engine_qndf(mechanism,destination)
    endswith(destination,".npz") || throw(ArgumentError("destination must end in .npz"))
    mkpath(dirname(abspath(destination)))
    before = engine_case_inputs(mechanism)
    runtime_before = engine_case_runtime(;enforce=true)
    open(destination*".invocation.toml","w") do io
        TOML.print(io,Dict("inputs_before"=>before,"runtime_before"=>runtime_before);sorted=true)
    end
    sink = NativeEngineQNDF.EngineFileSink(destination*".diagnostics")
    calculation = NativeEngineQNDF.solve_ic_engine_qndf(mechanism;diagnostic_sink=sink)
    # A failed post-run observation must not discard a completed calculation.
    after = try
        engine_case_inputs(mechanism)
    catch error
        Dict("capture_error"=>sprint(showerror,error))
    end
    runtime_after = try
        engine_case_runtime(;enforce=false)
    catch error
        Dict("checks_pass"=>false,"capture_error"=>sprint(showerror,error))
    end
    write_engine_case(calculation,destination,before,after,runtime_before,runtime_after)
    println("Saved complete engine history: ",destination)
    nothing
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("usage: ic_engine_qndf_case.jl MECHANISM.yaml OUTPUT.npz")
    export_ic_engine_qndf(ARGS[1],ARGS[2])
end
