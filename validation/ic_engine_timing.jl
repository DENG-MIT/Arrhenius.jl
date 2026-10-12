# Public complete QNDF engine timing. Caller dependencies remain external.
# julia --project=CALLER_ENV validation/ic_engine_timing.jl MECHANISM.yaml OUT.npz REPEATS MODE PRIOR_NATIVE.toml PRIOR_SHA256
# smoke is exactly first+one warm; controlled requires exactly nine warm calls.
const ENGINE_TIMING_IMPORT_SECONDS=@elapsed begin
    using Arrhenius, SHA, TOML, Dates, LinearAlgebra, Libdl, Pkg, Serialization
    import Arrhenius.NPZ
    include(joinpath(@__DIR__,"ic_engine_qndf_case.jl"))
end

engine_hash(path)=bytes2hex(sha256(read(path)))
function engine_tree_hashes(directory)
    isdir(directory) || error("missing required source directory: $directory")
    Dict(replace(relpath(joinpath(folder,name),directory),'\\'=>'/')=>engine_hash(joinpath(folder,name))
        for (folder,_,names) in walkdir(directory) for name in names)
end
function engine_timing_source_files(root)
    files=Dict("src/"*p=>h for (p,h) in engine_tree_hashes(joinpath(root,"src")))
    merge!(files,Dict("example/reactors/engine_qndf/"*p=>h for (p,h) in engine_tree_hashes(joinpath(root,"example/reactors/engine_qndf"))))
    for name in ("Project.toml","example/reactors/ic_engine.jl","example/reactors/ic_engine_setup.jl",
            "example/reactors/ic_engine_qndf_solver.jl","validation/ic_engine_timing.jl",
            "validation/ic_engine_timing.py","validation/ic_engine_source.py","validation/ic_engine_qndf_case.jl",
            "validation/ic_engine_qndf_case.py","validation/numerical_threads.jl","validation/benchmark_environment.py")
        files[name]=engine_hash(joinpath(root,name))
    end
    files
end
function engine_dependency_files()
    result=Dict{String,Any}()
    for (uuid,package) in Pkg.dependencies()
        path=package.source
        isnothing(path) && error("dependency source unavailable: $(package.name)")
        files=Dict{String,String}()
        for folder in ("src","ext","deps","lib")
            isdir(joinpath(path,folder)) || continue
            merge!(files,Dict(folder*"/"*p=>h for (p,h) in engine_tree_hashes(joinpath(path,folder))))
        end
        for name in ("Project.toml","Artifacts.toml")
            isfile(joinpath(path,name)) && (files[name]=engine_hash(joinpath(path,name)))
        end
        result[string(uuid)]=Dict("name"=>package.name,"version"=>string(package.version),
            "tree_hash"=>string(package.tree_hash),"path"=>realpath(path),"files"=>files)
    end
    result
end
function engine_timing_inputs(mechanism,root,receipt)
    basic=engine_case_inputs(mechanism)
    basic["source_inventory"]=engine_timing_source_files(root)
    basic["dependencies"]=engine_dependency_files()
    basic["julia_executable_sha256"]=engine_hash(joinpath(Sys.BINDIR,Sys.iswindows() ? "julia.exe" : "julia"))
    basic["file_paths"]=Dict("mechanism"=>abspath(mechanism),"sidecar"=>abspath(mechanism*".npz"),
        "caller_project"=>Base.active_project(),"caller_manifest"=>joinpath(dirname(Base.active_project()),"Manifest.toml"),
        "julia_executable"=>joinpath(Sys.BINDIR,Sys.iswindows() ? "julia.exe" : "julia"),"prior_library_receipt"=>abspath(receipt))
    basic["prior_library_receipt_sha256"]=engine_hash(receipt)
    basic
end

engine_mapped_libraries()=sort!(unique(realpath(p) for p in Libdl.dllist() if isfile(p)))
engine_numeric_library(p)=occursin(r"blas|mkl|accelerate|libgomp|libiomp|libomp"i,basename(p))
function engine_library_pins(receipt,expected_sha; mapped=engine_mapped_libraries())
    engine_hash(receipt)==expected_sha || error("prior independent library receipt hash differs")
    proof=TOML.parsefile(receipt)
    proof["checks_pass"] || error("prior library receipt must be a successful independent public run")
    expected=proof["runtime_after"]["loaded_numerical_library_sha256"]
    isempty(expected) && error("prior numerical-library inventory is empty")
    # Resolve only numerical basenames already pinned by the independent public
    # engine proof. In particular, MKL dispatch additions are not an open-ended
    # allowance for every sibling in an artifact directory.
    directories=unique(dirname(p) for p in mapped if engine_numeric_library(p))
    permitted=Dict{String,String}()
    for (name,hash) in expected
        basename(name)==name || error("prior library receipt must use basenames")
        matches=filter(isfile,[joinpath(directory,name) for directory in directories])
        isempty(matches) && error("cannot resolve pinned numerical library: $name")
        for path in matches
            engine_hash(path)==hash || error("independent numerical library bytes differ: $path")
            permitted[realpath(path)]=hash
        end
    end
    all(p->!engine_numeric_library(p) || get(permitted,p,"")==engine_hash(p),mapped) ||
        error("currently mapped numerical library absent from the prior independent proof")
    merge(Dict(p=>engine_hash(p) for p in mapped),permitted)
end
function engine_runtime_snapshot()
    settings=benchmark_julia_thread_settings(;enforce=false)
    # Also query mapped OpenMP runtimes and OpenBLAS instances directly. A
    # requested environment value alone is not an actual-thread observation.
    for path in engine_mapped_libraries()
        name=lowercase(basename(path))
        symbols=occursin("openblas",name) ? (:openblas_get_num_threads64_,:openblas_get_num_threads) :
            occursin(r"libgomp|libiomp|libomp",name) ? (:omp_get_max_threads,) : ()
        isempty(symbols) && continue
        handle=Libdl.dlopen(path)
        try
            query=C_NULL
            for symbol in symbols
                query=Libdl.dlsym(handle,symbol;throw_error=false)
                query==C_NULL || break
            end
            query==C_NULL && error("cannot query actual numerical threads: $path")
            settings["mapped_threads:"*basename(path)]=Int(ccall(query,Cint,()))
        finally
            Libdl.dlclose(handle)
        end
    end
    all(==(1),values(settings)) || error("mapped numerical backend must use one thread")
    Dict("checks_pass"=>true,"settings"=>settings,"mapped"=>engine_mapped_libraries(),
        "mapped_sha256"=>Dict(p=>engine_hash(p) for p in engine_mapped_libraries()),
        "environment"=>Dict(k=>get(ENV,k,"") for k in ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")))
end
function engine_require_runtime(before,after,pins;allow_first_load=false)
    before["environment"]==after["environment"] || error("thread environment changed")
    all(==(1),values(after["settings"])) || error("actual numerical thread settings changed")
    initial,final=Set(before["mapped"]),Set(after["mapped"])
    (allow_first_load ? issubset(initial,final) : initial==final) || error("loaded-library membership changed")
    all(haskey(pins,p) && pins[p]==h for (p,h) in after["mapped_sha256"]) || error("unrecognized or changed mapped library")
    all(isfile(p) && engine_hash(p)==h for (p,h) in pins) || error("pre-pinned library bytes changed")
    true
end

# Only this explicitly identified performance field is nondeterministic.
const ENGINE_WALL_CLOCK_FIELDS=("seconds_including_first_specialization",)
function engine_deterministic_records(records)
    [Dict(k=>v for (k,v) in row if !(k in ENGINE_WALL_CLOCK_FIELDS)) for row in records]
end
function engine_replay_payload(calculation)
    Dict("archive"=>engine_case_archive_data(calculation),"summary"=>calculation.summary,
        "checks"=>calculation.checks,"records"=>engine_deterministic_records(calculation.records),
        "diagnostics"=>calculation.diagnostics,
        "integrals"=>calculation.result.integrals,"quadrature_integrals"=>calculation.result.quadrature_integrals,
        "coarse_integrals"=>calculation.result.coarse_integrals,
        "display_times"=>calculation.result.times,
        "integration_points"=>calculation.result.integration_points)
end
function engine_exact_equal(a,b)
    typeof(a)===typeof(b) || return false
    if a isa AbstractDict
        Set(keys(a))==Set(keys(b)) && all(engine_exact_equal(a[k],b[k]) for k in keys(a))
    elseif a isa AbstractArray
        axes(a)==axes(b) || return false
        if isbitstype(eltype(a)) && a isa Array && b isa Array
            reinterpret(UInt8,vec(a))==reinterpret(UInt8,vec(b))
        else
            all(engine_exact_equal(x,y) for (x,y) in zip(a,b))
        end
    elseif a isa NamedTuple || a isa Tuple
        all(engine_exact_equal(x,y) for (x,y) in zip(a,b))
    else
        isequal(a,b)
    end
end
function engine_require_replay(actual,expected,label)
    engine_exact_equal(actual,expected) || error("$label differs from first complete engine output/counters")
    true
end
function engine_capture_completed(calculation,path)
    NPZ.npzwrite(path,engine_case_archive_data(calculation))
    payload=engine_replay_payload(calculation)
    serialize(path*".replay.jls",payload)
    payload
end
function engine_require_complete(calculation)
    records=calculation.records
    length(records)==25 || error("25 complete intervals required")
    calculation.summary["end_time_s"]==.16 || error("eight complete revolutions required")
    [r["index"] for r in records]==collect(1:25) || error("interval indices changed")
    [r["integrated_coordinate_count"] for r in records]==vcat(fill(8,3),fill(105,22)) || error("coordinate activation changed")
    calculation.checks["checks_pass"] || error("public QNDF physical checks failed")
    NativeEngineQNDF.engine_qndf_summary_checks(calculation) || error("engine summary gates failed")
    true
end
@noinline engine_timing_calculation(mechanism)=NativeEngineQNDF.solve_ic_engine_qndf(mechanism)

function engine_timing_main(args;calculate=engine_timing_calculation)
    length(args)==6 || error("MECHANISM.yaml OUT.npz REPEATS MODE PRIOR_NATIVE.toml PRIOR_SHA256")
    mechanism,destination,repeat_arg,mode,receipt,receipt_sha=args
    repeats=parse(Int,repeat_arg)
    mode in ("smoke","controlled") || error("unknown engine timing mode")
    repeats==(mode=="smoke" ? 1 : 9) || error("smoke requires one warm call; controlled requires nine")
    endswith(destination,".npz") || error("destination must end in .npz")
    root=realpath(dirname(@__DIR__))
    realpath(pkgdir(Arrhenius))==root || error("public helpers and package must come from one frozen checkout")
    paths=[destination,replace(destination,r"\.npz$"=>".toml"),destination*".timing.toml"]
    any(ispath,paths) && error("preserve existing attempts")
    mkpath(dirname(abspath(destination)))
    benchmark_julia_thread_settings(;enforce=true)
    keys=Sys.isapple() ? ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS") :
        ("JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")
    all(k->get(ENV,k,"")=="1",keys) || error("all numerical thread environment variables must be 1")
    before=engine_timing_inputs(mechanism,root,receipt);pins=engine_library_pins(receipt,receipt_sha);runtime=engine_runtime_snapshot()
    report=Dict{String,Any}("complete"=>false,"performance_qualified"=>false,"mode"=>mode,"warm_repetitions"=>repeats,
        "scope"=>"fresh public QNDF construction, all eight revolutions/25 intervals, outputs and required integrals",
        "public_callable"=>"NativeEngineQNDF.solve_ic_engine_qndf","input_before"=>before,"library_pins"=>pins,
        "runtime_before"=>runtime,"import_seconds"=>ENGINE_TIMING_IMPORT_SECONDS,"wall_clock_exclusions"=>collect(ENGINE_WALL_CLOCK_FIELDS),
        "first_call_scope"=>"invokelatest of the complete public wrapper; includes wrapper specialization, excludes earlier driver imports and guard compilation",
        "first_seconds"=>0.,"warm_seconds"=>Float64[],"runs"=>Dict{String,Any}[],"julia_version"=>string(VERSION),
        "cpu"=>Sys.cpu_info()[1].model,"system"=>string(Sys.KERNEL),"kernel_release"=>Sys.isunix() ? readchomp(`uname -r`) : "unknown",
        "native_policy"=>"public QNDF rtol1e-12, source-derived conservative dynamic mass weights, existing per-component scaling and 105-coordinate RMS")
    record_path=destination*".timing.toml"
    checkpoint()=open(io->TOML.print(io,report;sorted=true),record_path,"w")
    checkpoint();first_payload=nothing
    try
        for index in 0:repeats
            label=index==0 ? "first" : "warm-$index"
            path=index==0 ? destination : replace(destination,r"\.npz$"=>".warm-$index.npz")
            ispath(path) && error("preserve existing $label artifact")
            inputs=engine_timing_inputs(mechanism,root,receipt)
            inputs==before || error("source/input/dependency membership or bytes changed before $label")
            prior=engine_runtime_snapshot()
            engine_require_runtime(runtime,prior,pins)
            started=time_ns()
            calculation=Base.invokelatest(calculate,mechanism)
            seconds=(time_ns()-started)*1e-9
            # Save completed arrays and complete deterministic payload BEFORE
            # any replay, input, runtime or physical postcondition can throw.
            payload=engine_capture_completed(calculation,path)
            run=Dict{String,Any}("label"=>label,"seconds"=>seconds,"archive"=>abspath(path),"archive_sha256"=>engine_hash(path),
                "payload_sha256"=>engine_hash(path*".replay.jls"),"runtime_before"=>prior,
                "checks_pass"=>false,"replay_pass"=>false)
            push!(report["runs"],run);checkpoint()
            after=try
                engine_timing_inputs(mechanism,root,receipt)
            catch error
                Dict("capture_error"=>sprint(showerror,error))
            end
            observed=try
                engine_runtime_snapshot()
            catch error
                Dict("checks_pass"=>false,"capture_error"=>sprint(showerror,error))
            end
            run["inputs"]=after;run["runtime_after"]=observed;checkpoint()
            # The existing exporter record remains consumable by the current
            # independent accepted-history comparator.
            write_engine_case(calculation,path,inputs,after,prior,observed)
            # The existing exporter rewrites NPZ outside the timer. Pin its
            # final bytes (ZIP metadata need not be identical on the rewrite).
            run["archive_sha256"]=engine_hash(path)
            run["metadata_sha256"]=engine_hash(replace(path,r"\.npz$"=>".toml"));checkpoint()
            isfinite(seconds) && seconds>0 || error("invalid native call duration")
            engine_require_complete(calculation)
            inputs==after==before || error("source/input/dependency bytes changed in $label")
            engine_require_runtime(prior,observed,pins;allow_first_load=index==0)
            run["checks_pass"]=true
            if index==0
                first_payload=payload;report["first_seconds"]=seconds
            else
                engine_require_replay(payload,first_payload,label)
                push!(report["warm_seconds"],seconds)
            end
            run["replay_pass"]=true;runtime=observed
            report["runtime_after"]=observed;checkpoint()
            println("Engine ",label," complete: ",seconds," s; exact output/counter and physical checks pass");flush(stdout)
        end
        report["input_after"]=engine_timing_inputs(mechanism,root,receipt)
        report["input_after"]==before || error("final input inventory changed")
        report["complete"]=true
    catch error
        report["error"]=sprint(showerror,error)
        if error isa NativeEngineQNDF.EngineCalculationFailure
            serialize(destination*".failure.jls",error.diagnostics)
        end
        rethrow()
    finally
        checkpoint()
    end
    nothing
end
if abspath(PROGRAM_FILE)==@__FILE__
    engine_timing_main(ARGS)
end
