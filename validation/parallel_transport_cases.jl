# julia --threads=4 --project=. validation/parallel_transport_cases.jl INPUT_DIR OUTPUT.npz
using Arrhenius, Arrhenius.NPZ, LinearAlgebra, SHA, Libdl
include(joinpath(@__DIR__,"..","example","transport","multiprocessing_viscosity.jl"))

function source_inventory(root)
    files = sort!([replace(relpath(joinpath(dir,name),root),'\\'=>'/')
        for (dir,_,names) in walkdir(root) for name in names])
    return files,[bytes2hex(sha256(read(joinpath(root,file)))) for file in files]
end

function main(args)
    length(args)==2 || error("supply INPUT_DIR and OUTPUT.npz")
    Threads.nthreads()==4 || error("launch Julia with --threads=4")
    BLAS.set_num_threads(1)
    directory,destination = args
    paths = (mechanism=joinpath(directory,"gri30.yaml"),sidecar=joinpath(directory,"gri30.yaml.npz"),
        transport=joinpath(directory,"gri30.yaml.multicomponent.npz"),
        example=joinpath(@__DIR__,"..","example","transport","multiprocessing_viscosity.jl"),harness=@__FILE__)
    filehash(path) = bytes2hex(sha256(read(path)))
    hashes = map(filehash,paths)
    core_root = joinpath(pkgdir(Arrhenius),"src")
    core = source_inventory(core_root)
    result = parallel_transport_calculation(paths.mechanism,paths.transport)
    data = Dict{String,Any}(String(key)=>getproperty(result,key) for key in keys(result))
    data["native_checks_pass"] = [false]
    npzwrite(destination,data)
    @assert isequal(result.conductivity_parallel,result.conductivity_serial)
    @assert isequal(result.viscosity_parallel,result.viscosity_serial)
    for property in (:conductivity_parallel,:conductivity_serial,:viscosity_parallel,:viscosity_serial)
        values = getproperty(result,property)
        @assert length(values)==5000 && all(isfinite,values) && all(>(0),values)
    end
    @assert map(filehash,paths)==hashes && source_inventory(core_root)==core
    @assert Threads.nthreads()==4 && BLAS.get_num_threads()==1
    @assert !any(path->occursin(r"(?i)cantera|coolprop|python",path),Libdl.dllist())
    utf8(value) = collect(codeunits(string(value)))
    for key in keys(paths)
        data[String(key)*"_sha256_utf8"] = utf8(getproperty(hashes,key))
    end
    cpu = Sys.isapple() ? readchomp(`sysctl -n machdep.cpu.brand_string`) :
        Sys.islinux() ? strip(split(first(filter(l->startswith(l,"model name"),readlines("/proc/cpuinfo"))),':';limit=2)[2]) :
        Sys.cpu_info()[1].model
    data["native_checks_pass"] = [true]
    data["julia_threads"],data["blas_threads"] = [Threads.nthreads()],[BLAS.get_num_threads()]
    data["core_source_paths_utf8"],data["core_source_sha256_utf8"] = utf8(join(core[1],"\n")),utf8(join(core[2],"\n"))
    data["julia_version_utf8"],data["cpu_utf8"] = utf8(VERSION),utf8(cpu)
    data["system_utf8"] = utf8(Sys.isapple() ? "Darwin" : Sys.islinux() ? "Linux" : "Windows")
    data["kernel_release_utf8"] = utf8(Sys.iswindows() ? string(Sys.windows_version()) : readchomp(`uname -r`))
    npzwrite(destination,data)
    println("PASS: all four 5,000-temperature sweeps; serial/parallel outputs are bitwise identical.")
end
main(ARGS)
