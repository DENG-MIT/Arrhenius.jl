# Complete warm mixing-network calculation on prepared mechanisms.
using Arrhenius, LinearAlgebra, SHA, TOML, Libdl, Statistics
import Arrhenius.NPZ
length(ARGS)==4 || error("SOURCE_ROOT INPUT_DIR VALIDATED_NATIVE_DIR OUTPUT_DIR")
include(joinpath(ARGS[1],"example","reactors","mixing_solver.jl"))
include(joinpath(ARGS[1],"validation","numerical_threads.jl"))

function mixing_timing(source,input,reference,output)
    ispath(output) && error("refusing existing timing output")
    mkpath(output)
    filehash(path)=bytes2hex(sha256(read(path)))
    paths=[joinpath(input,name) for name in ("gri30.yaml","gri30.yaml.npz","air.yaml","air.yaml.npz")]
    append!(paths,[joinpath(source,"example","reactors",name) for name in ("mixing_solver.jl","network_cases.jl")])
    append!(paths,[@__FILE__,joinpath(source,"validation","numerical_threads.jl"),
        joinpath(reference,"native.npz"),joinpath(reference,"native.toml")])
    for (folder,_,files) in walkdir(joinpath(pkgdir(Arrhenius),"src")),name in files
        endswith(name,".jl") && push!(paths,joinpath(folder,name))
    end
    inventory()=Dict(path=>filehash(path) for path in paths)
    libraries()=Dict(basename(path)=>filehash(path) for path in Libdl.dllist()
        if isfile(path) && occursin(r"blas|mkl|accelerate|iomp|libgomp"i,basename(path)))
    before=inventory()
    threads_before=benchmark_julia_thread_settings(;enforce=true)
    gas=CreateSolution(joinpath(input,"gri30.yaml"));air=CreateSolution(joinpath(input,"air.yaml"))
    expected=NPZ.npzread(joinpath(reference,"native.npz"))
    meta=TOML.parsefile(joinpath(reference,"native.toml"))
    meta["converged"] && meta["inputs_unchanged"] || error("validated native calculation required")
    meta["native_sha256"]==filehash(joinpath(reference,"native.npz")) || error("native archive changed")
    meta["candidate_sha256"]==filehash(joinpath(source,"example","reactors","mixing_solver.jl")) || error("solver differs from validated calculation")
    for name in ("gri30.yaml","gri30.yaml.npz","air.yaml","air.yaml.npz")
        matches=[digest for (path,digest) in meta["inputs_before"] if basename(path)==name]
        length(matches)==1 && only(matches)==filehash(joinpath(input,name)) || error("mechanism differs from validated calculation")
    end
    function check(result)
        result.converged && result.physical_residual<=1e-9 || error("mixer failed convergence")
        result.iterations==meta["iterations"] || error("iteration count changed")
        isequal(result.state,expected["state"]) || error("stationary state changed")
        isequal(network_state(result.network),expected["initial_state"]) || error("initial state changed")
        diagnostic=network_diagnostics(result.network,result.state)
        isequal(diagnostic.mass_flow_rates,expected["flows"]) || error("boundary flow changed")
        isequal(diagnostic.nodes.mixer.mass_fractions,expected["Y"]) || error("species changed")
        nothing
    end
    # Specialize and warm the same complete call before timing.
    check(solve_mixing_network(gas,air))
    libraries_before=libraries()
    arrays=Dict{String,Any}();times=Float64[];iterations=Int[];residuals=Float64[]
    report=Dict{String,Any}("passed"=>false,"source_before"=>before,"threads_before"=>threads_before,
        "libraries_before"=>libraries_before,"scope"=>"fresh network construction and full steady solve; prepared mechanism loading and validation excluded",
        "julia_version"=>string(VERSION),"system"=>string(Sys.KERNEL),"cpu"=>Sys.cpu_info()[1].model,
        "kernel_release"=>Sys.iswindows() ? string(Sys.windows_version()) : readchomp(`uname -r`))
    try
        for i in 1:9
            start=time_ns()
            result=solve_mixing_network(gas,air)
            push!(times,(time_ns()-start)*1e-9)
            arrays["state_$i"]=copy(result.state)
            push!(iterations,result.iterations);push!(residuals,result.physical_residual)
            NPZ.npzwrite(joinpath(output,"states.npz"),arrays)
            check(result)
        end
        report["source_after"]=inventory()
        report["threads_after"]=benchmark_julia_thread_settings(;enforce=false)
        report["libraries_after"]=libraries()
        report["source_after"]==before || error("inputs changed during timing")
        report["libraries_after"]==libraries_before || error("loaded libraries changed during timing")
        any(path->occursin(r"cantera|coolprop|python"i,path),Libdl.dllist()) && error("foreign solver loaded")
        report["passed"]=true
    finally
        report["warm_seconds"]=times
        report["median_seconds"]=isempty(times) ? NaN : median(times)
        report["iterations"]=iterations;report["physical_residuals_per_s"]=residuals
        open(joinpath(output,"timing.toml"),"w") do io;TOML.print(io,report;sorted=true);end
    end
    println("Complete mixer median: ",median(times)," s; 9 exact output replays")
end
mixing_timing(ARGS...)
