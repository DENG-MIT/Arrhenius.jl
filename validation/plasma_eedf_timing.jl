
const script_started=time_ns()
using Arrhenius,LinearAlgebra,Libdl,SHA,TOML,NPZ,Statistics
source_root,model_path,reference,output=abspath.(ARGS)
ispath(output) && error("fresh output directory required")
mkpath(output)
realpath(dirname(dirname(pathof(Arrhenius))))==realpath(source_root) || error("wrong package")
include(joinpath(source_root,"validation","numerical_threads.jl"))
threads_before=benchmark_julia_thread_settings()
include(joinpath(source_root,"example","thermodynamics","plasma_eedf.jl"))
digest(p)=bytes2hex(sha256(read(p)))
source_paths=[joinpath(source_root,"src",name) for name in readdir(joinpath(source_root,"src")) if endswith(name,".jl")]
input_paths=[source_paths;model_path;reference;joinpath(source_root,"Manifest.toml");joinpath(source_root,"Project.toml");
    joinpath(source_root,"example","thermodynamics","plasma_eedf.jl");joinpath(source_root,"validation","numerical_threads.jl")]
input_hashes=Dict(path=>digest(path) for path in input_paths)
input_hashes[model_path]=="be36740e0db01d996668b2c3ff6cbf12a0f296155086a16b7335d26df43f3ce2" || error("wrong model")
input_hashes[reference]=="8282b42f1c715d77b8fa33e6cacb2622edc810aa2eb253b3cd45849ddb22d92c" || error("wrong reference")
function measure_complete_calls(path)
    timings=Float64[];saved=EEDFResult[];cold=0.0
    for rep in 0:9
        started=time_ns()
        result=Base.invokelatest(air_eedf,path) # fresh loader, state and complete solve
        elapsed=(time_ns()-started)*1e-9
        rep==0 && (cold=(time_ns()-script_started)*1e-9)
        push!(timings,elapsed);push!(saved,result)
        npzwrite(joinpath(output,"native-$rep.npz"),Dict("edges"=>result.edges,"edge_eedf"=>result.edge_eedf,
            "centers"=>result.centers,"center_eedf"=>result.center_eedf,"errors"=>result.errors,"deltas"=>result.deltas,"mobility"=>[result.mobility]))
        result.converged && result.iterations<=200 && last(result.errors)<1e-5 || error("convergence gate failed")
        result.edges==collect(0.:40.) && all(isfinite,result.edge_eedf) && all(result.edge_eedf.>=0) || error("shape/finite gate failed")
        println((;rep,elapsed,iterations=result.iterations));flush(stdout)
    end
    return timings,saved,cold
end
times,saved,cold=measure_complete_calls(model_path)
ref=npzread(reference);x=vec(ref["electron_energy_levels"]);f=vec(ref["electron_energy_distribution"])
trap(y,x)=sum((x[2:end].-x[1:end-1]).*(y[2:end].+y[1:end-1])./2)
weight=sqrt.(x);normref=trap(weight.*f,x);meanref=trap(x.*weight.*f,x)/normref
comparisons=Dict{String,Any}[]
for (rep,result) in enumerate(saved)
    weighted=trap(weight.*abs.(result.edge_eedf.-f),x)/normref
    meanvalue=trap(x.*weight.*result.edge_eedf,x)/trap(weight.*result.edge_eedf,x)
    pointwise=maximum(abs.(result.edge_eedf.-f)./(1e-10.+1e-4.*abs.(f)))
    normalization=Arrhenius._eedf_norm(result.center_eedf,result.centers)
    passed=weighted<=1e-4 && abs(meanvalue/meanref-1)<=1e-4 && pointwise<=1 && abs(normalization-1)<=1e-12
    push!(comparisons,Dict("repetition"=>rep-1,"passed"=>passed,"weighted_error"=>weighted,"pointwise_error_ratio"=>pointwise,"mean_energy"=>meanvalue,"normalization"=>normalization))
end
npzwrite(joinpath(output,"timings.npz"),Dict("seconds"=>times,"iterations"=>[r.iterations for r in saved]))
threads_after=benchmark_julia_thread_settings(;enforce=false)
unchanged=input_hashes==Dict(path=>digest(path) for path in input_paths)
passed=all(x["passed"] for x in comparisons) && unchanged
receipt=Dict("passed"=>passed,"qualified"=>false,"formal"=>false,"first_call_seconds"=>times[1],
    "script_start_to_first_result_seconds"=>cold,"warm_median_seconds"=>median(times[2:end]),"repetitions"=>10,
    "scope"=>"Fresh public air_eedf(path): YAML load, model and state construction, complete solve and returned arrays; output/validation outside timers",
    "seconds"=>times,"comparisons"=>comparisons,"source_inputs_unchanged"=>unchanged,
    "input_hashes"=>input_hashes,"threads_before"=>threads_before,"threads_after"=>threads_after)
open(joinpath(output,"receipt.toml"),"w") do io;TOML.print(io,receipt;sorted=true);end
passed || error("EEDF complete-call correctness failed")
println((;passed,warm_median=median(times[2:end]),first_call=times[1],cold));flush(stdout)
