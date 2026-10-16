# Private diagnostic timing driver for the pinned transport/multiprocessing_viscosity workload.
# julia --threads=4 --project=<SOURCE_ROOT> diagnostic_timing.jl SOURCE_ROOT INPUT_DIR OUTPUT_NPZ
# INPUT_DIR holds the staged gri30.yaml, gri30.yaml.npz and gri30.yaml.multicomponent.npz.
# Diagnostic only: raw timings are recorded and are never performance-qualified.
using Arrhenius, Arrhenius.NPZ, LinearAlgebra, SHA, Libdl
length(ARGS)==3 || error("supply SOURCE_ROOT INPUT_DIR OUTPUT_NPZ")
include(joinpath(ARGS[1],"example","transport","multiprocessing_viscosity.jl"))

# Compiled once via Base.invokelatest after the cold first call: one untimed
# complete warmup, then nine timed complete repetitions of the public function.
# Every measured call creates its own phases and workspaces inside
# parallel_transport_calculation; nothing is hoisted out of the timer.
function _complete_warm_loop(mechanism,transport;repetitions::Integer=9)
    warmup = parallel_transport_calculation(mechanism,transport)
    seconds = Vector{Float64}(undef,repetitions)
    results = Vector{Any}(undef,repetitions)
    for rep in 1:repetitions
        start = time_ns()
        results[rep] = parallel_transport_calculation(mechanism,transport)
        seconds[rep] = (time_ns()-start)/1e9
    end
    return warmup,seconds,results
end

function source_inventory(root)
    files = sort!([replace(relpath(joinpath(dir,name),root),'\\'=>'/')
        for (dir,_,names) in walkdir(root) for name in names])
    return files,[bytes2hex(sha256(read(joinpath(root,file)))) for file in files]
end

function main(args)
    length(args)==3 || error("supply SOURCE_ROOT INPUT_DIR OUTPUT_NPZ")
    source_root,input_dir,destination = args
    Threads.nthreads()==4 || error("launch Julia with --threads=4")
    BLAS.set_num_threads(1)
    accelerate_mode = Sys.isapple() ? begin
        status = ccall((:BLASSetThreading,"/System/Library/Frameworks/Accelerate.framework/Accelerate"),Cint,(Cuint,),1)
        status==0 || error("Accelerate single-thread setting failed")
        Int(ccall((:BLASGetThreading,"/System/Library/Frameworks/Accelerate.framework/Accelerate"),Cuint,()))
    end : -1
    paths = (mechanism=joinpath(input_dir,"gri30.yaml"),sidecar=joinpath(input_dir,"gri30.yaml.npz"),
        transport=joinpath(input_dir,"gri30.yaml.multicomponent.npz"),
        example=joinpath(source_root,"example","transport","multiprocessing_viscosity.jl"),driver=@__FILE__)
    filehash(path) = bytes2hex(sha256(read(path)))
    hashes_before = map(filehash,paths)
    core_root = joinpath(pkgdir(Arrhenius),"src")
    core_before = source_inventory(core_root)

    first_start = time_ns()
    initial = Base.invokelatest(parallel_transport_calculation,paths.mechanism,paths.transport)
    first_seconds = (time_ns()-first_start)/1e9
    fields = (:T,:conductivity_parallel,:conductivity_serial,:viscosity_parallel,:viscosity_serial)
    data = Dict{String,Any}(String(key)=>getproperty(initial,key) for key in fields)
    data["first_call_seconds"] = [first_seconds]
    data["native_checks_pass"] = [false]
    data["complete_timing_workload"] = [true]
    data["performance_qualified"] = [false]
    npzwrite(destination,data)
    warmup,warm_seconds,warm_results = Base.invokelatest(_complete_warm_loop,paths.mechanism,paths.transport)
    data["warm_seconds"] = warm_seconds
    npzwrite(destination,data)

    # Preserve a mismatching output before rejecting any replay, outside timers.
    for (index,candidate) in enumerate((warmup,warm_results...))
        if !all(field->isequal(getproperty(candidate,field),getproperty(initial,field)),fields)
            for field in fields
                data["failed_"*String(field)] = getproperty(candidate,field)
            end
            data["failed_replay_index"] = [index]
            npzwrite(destination,data)
            error("repetition output differs from the first complete call")
        end
    end

    @assert isequal(initial.conductivity_parallel,initial.conductivity_serial)
    @assert isequal(initial.viscosity_parallel,initial.viscosity_serial)
    for property in (:conductivity_parallel,:conductivity_serial,:viscosity_parallel,:viscosity_serial)
        values = getproperty(initial,property)
        @assert length(values)==5000 && all(isfinite,values) && all(>(0),values)
    end
    @assert map(filehash,paths)==hashes_before && source_inventory(core_root)==core_before
    @assert Threads.nthreads()==4 && BLAS.get_num_threads()==1
    @assert !Sys.isapple() || accelerate_mode==1
    @assert !Sys.isapple() || Int(ccall((:BLASGetThreading,"/System/Library/Frameworks/Accelerate.framework/Accelerate"),Cuint,()))==accelerate_mode
    @assert !any(path->occursin(r"(?i)cantera|coolprop|python",path),Libdl.dllist())

    utf8(value) = collect(codeunits(string(value)))
    for key in keys(paths)
        data[String(key)*"_sha256_utf8"] = utf8(getproperty(hashes_before,key))
    end
    cpu = Sys.isapple() ? readchomp(`sysctl -n machdep.cpu.brand_string`) :
        Sys.islinux() ? strip(split(first(filter(l->startswith(l,"model name"),readlines("/proc/cpuinfo"))),':';limit=2)[2]) :
        Sys.cpu_info()[1].model
    data["native_checks_pass"] = [true]
    data["julia_threads"],data["blas_threads"] = [Threads.nthreads()],[BLAS.get_num_threads()]
    data["accelerate_threading_mode"] = [accelerate_mode]
    data["pressure_Pa"] = [101325.0]
    data["core_source_paths_utf8"],data["core_source_sha256_utf8"] = utf8(join(core_before[1],"\n")),utf8(join(core_before[2],"\n"))
    data["julia_version_utf8"],data["cpu_utf8"] = utf8(VERSION),utf8(cpu)
    data["system_utf8"] = utf8(Sys.isapple() ? "Darwin" : Sys.islinux() ? "Linux" : "Windows")
    data["kernel_release_utf8"] = utf8(Sys.iswindows() ? string(Sys.windows_version()) : readchomp(`uname -r`))
    data["composition_utf8"] = utf8("CH4:1.0, O2:1.0, N2:3.76")
    data["transport_model_utf8"] = utf8("multicomponent")
    npzwrite(destination,data)
    println("PASS: complete native timing workload; first call plus nine warm repetitions recorded.")
end
main(ARGS)
