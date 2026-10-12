# Complete published transport sequences; no reference states enter a solve.
using Arrhenius,LinearAlgebra,NPZ,Libdl,SHA,TOML
include("numerical_threads.jl")
thread_checks=Dict("before_first"=>benchmark_julia_thread_settings())
source_root=dirname(dirname(pathof(Arrhenius)))
function source_hashes(root)
    Dict(relpath(joinpath(folder,file),root)=>bytes2hex(sha256(read(joinpath(folder,file))))
        for (folder,_,files) in walkdir(joinpath(root,"src")) for file in files if endswith(file,".jl"))
end
initial_source_hashes=source_hashes(source_root)
parameters,output,case=ARGS[1:3]
case in ("free","burner","fixed","ion-free","ion-burner") || error("unknown flame case")
reps=parse(Int,ARGS[4]);reps>=5 || error("at least five warm repetitions required")
Threads.nthreads()==1 && BLAS.get_num_threads()==1 || error("single-thread benchmark required")
clock_diagnostics=length(ARGS)>=5 && ARGS[5]=="clock-diagnostics"
clock_diagnostics && !(Sys.islinux() && Sys.WORD_SIZE==64) && error("clock diagnostics require 64-bit Linux")
const on_linux64=Sys.islinux() && Sys.WORD_SIZE==64
const elapsed_clock=on_linux64 ? "CLOCK_MONOTONIC_RAW" : "monotonic"
function benchmark_elapsed_ns()::UInt64
    if Sys.islinux() && Sys.WORD_SIZE==64
        spec=Ref{NTuple{2,Clong}}()
        ccall(:clock_gettime,Cint,(Cint,Ref{NTuple{2,Clong}}),4,spec)==0 ||
            error("clock_gettime(CLOCK_MONOTONIC_RAW) failed")
        return UInt64(spec[][1])*UInt64(1_000_000_000)+UInt64(spec[][2])
    end
    return time_ns()
end
function flame_clock_observation()
    cpu=ccall(:sched_getcpu,Cint,())
    cpu>=0 || error("sched_getcpu failed")
    return (benchmark_elapsed_ns(),ccall(:clock,Clong,())/1e6,cpu)
end
mkpath(output)
ion_case=startswith(case,"ion-")
source_case=ion_case ? case[5:end] : case
if ion_case
    mechanism=joinpath(parameters,"gri30_ion.yaml")
    gas=CreateSolution(mechanism)
    shared_helper=joinpath(source_root,"example","flames","source_ion_flame_sequence.jl")
    include(shared_helper)
    realpath(String(which(run_ion_source_sequence,Tuple{typeof(gas),typeof(source_case)}).file))==realpath(shared_helper) ||
        error("ion source calculation must come from the public example helper")
    modes=["frozen","field"]
    fields=["grid","T","Y","X","E","velocity","rho","qdot","inlet_Y","state","P"]
else
    mechanism=joinpath(parameters,case=="fixed" ? "gri30.yaml" : "h2o2.yaml")
    gas=CreateSolution(mechanism);data=MultiTransportData(mechanism*".multicomponent.npz",gas;mechanism)
    profile=case=="fixed" ? npzread(joinpath(parameters,"fixed-profile.npz")) : nothing
    shared_helper=joinpath(source_root,"example","flames","source_flame_sequence.jl")
    include(shared_helper)
    realpath(String(which(run_source_sequence,Tuple{typeof(gas),typeof(data),typeof(case),typeof(profile)}).file))==realpath(shared_helper) ||
        error("source calculation must come from the public example helper")
    if case=="fixed"
        original_profile=source_flame_temperature_profile()
        all(profile[key]==original_profile[key] for key in ("positions","temperatures")) ||
            error("public prescribed temperature profile differs from pinned source data")
    end
    modes=case=="free" ? ["mass","mass-soret","multi","multi-soret"] : ["mole","multi"]
    fields=["grid","T","Y","velocity","inlet_Y","state","P"]
end
seconds=zeros(reps+1,length(modes));points=zeros(Int,size(seconds));snapshots=nothing;first_profiles=nothing
replay_errors=zeros(reps+1,length(modes),length(fields));hashes=String[]
automatic_gc_seconds=zeros(reps+1);warmup_gc_seconds=Ref(0.0)
clock_observations=zeros(clock_diagnostics ? reps+1 : 0,4)
for repetition in 0:reps
    gc_start_ns=Base.gc_time_ns()
    clock_before=clock_diagnostics ? flame_clock_observation() : nothing
    elapsed,nodes,saved=ion_case ? run_ion_source_sequence(gas,source_case;clock_ns=benchmark_elapsed_ns) :
        run_source_sequence(gas,data,case,profile;save_profiles=true,clock_ns=benchmark_elapsed_ns)
    if clock_diagnostics
        after=flame_clock_observation()
        clock_observations[repetition+1,:].=((after[1]-clock_before[1])/1e9,after[2]-clock_before[2],clock_before[3],after[3])
    end
    automatic_gc_seconds[repetition+1]=(Base.gc_time_ns()-gc_start_ns)/1e9
    thread_checks["after_repetition_"*string(repetition)]=benchmark_julia_thread_settings(;enforce=false)
    seconds[repetition+1,:].=elapsed;points[repetition+1,:].=nodes
    repetition==0 && (global first_profiles=saved)
    for (stage,mode) in enumerate(modes)
        buffer=IOBuffer()
        for (field,key) in enumerate(fields)
            current,baseline=saved[mode][key],first_profiles[mode][key]
            size(current)==size(baseline) || error("replay dimensions changed for $case $mode $key")
            all(isfinite,current) || error("nonfinite replay output")
            all(abs.(current.-baseline) .<= 1e-14 .+ 1e-12.*abs.(baseline)) || error("numerical replay failed for $case $mode $key")
            replay_errors[repetition+1,stage,field]=maximum(abs,current.-baseline)
            write(buffer,key);write(buffer,Int64.(collect(size(current))));write(buffer,vec(current))
        end
        push!(hashes,bytes2hex(SHA.sha256(take!(buffer))))
        repetition==0 && npzwrite(joinpath(output,"$case-$mode-first.npz"),saved[mode])
    end
    repetition==reps && (global snapshots=saved)
    println((;case,repetition,seconds=sum(elapsed),stages=elapsed,points=nodes));flush(stdout)
    # Exclude garbage left by the excluded warmup; automatic GC stays enabled
    # throughout all measured repetitions. No per-repetition collection.
    if repetition==0
        warmup_gc_seconds[]=@elapsed GC.gc()
    end
end
for (mode,values) in snapshots
    npzwrite(joinpath(output,"$case-$mode-0.npz"),values)
end
initial_source_hashes==source_hashes(source_root) || error("native source files changed during calculation")
tomlbytes(value)=collect(codeunits(sprint(io->TOML.print(io,value;sorted=true))))
library_hashes=Dict(realpath(path)=>bytes2hex(sha256(read(path))) for path in Libdl.dllist() if isfile(path))
npzwrite(joinpath(output,"timings.npz"),Dict("stage_seconds"=>seconds,"stage_points"=>points,
    "automatic_gc_seconds"=>automatic_gc_seconds,"warmup_gc_seconds"=>[warmup_gc_seconds[]],
    "clock_diagnostics"=>[clock_diagnostics],"clock_observations"=>clock_observations,
    "elapsed_clock_utf8"=>collect(codeunits(elapsed_clock)),
    "shared_calculation_path_utf8"=>collect(codeunits(realpath(shared_helper))),
    "shared_calculation_sha256_utf8"=>collect(codeunits(bytes2hex(sha256(read(shared_helper))))),
    "thread_checks_toml_utf8"=>tomlbytes(thread_checks),
    "source_hashes_toml_utf8"=>tomlbytes(initial_source_hashes),"source_hashes_unchanged"=>[true],
    "library_hashes_toml_utf8"=>tomlbytes(library_hashes),
    "julia_executable_sha256_utf8"=>collect(codeunits(bytes2hex(sha256(read(joinpath(Sys.BINDIR,"julia")))))),
    "replay_checked_stages"=>[(reps+1)*length(modes)],"replay_max_abs"=>replay_errors,
    "replay_hashes_utf8"=>collect(codeunits(join(hashes,"\n"))),"replay_fields_utf8"=>collect(codeunits(join(fields,"\n"))),
    "julia_version_utf8"=>collect(codeunits(string(VERSION))),"julia_threads"=>[Threads.nthreads()],"blas_threads"=>[BLAS.get_num_threads()],
    "accelerate_threading"=>[get(thread_checks["after_repetition_"*string(reps)],"accelerate_threading_mode",-1)],
    "kernel_utf8"=>collect(codeunits(readchomp(`uname -r`))),"machine_utf8"=>collect(codeunits(string(Sys.MACHINE))),
    "package_path_utf8"=>collect(codeunits(pathof(Arrhenius))),
    "blas_config_utf8"=>collect(codeunits(string(BLAS.get_config()))),
    "loaded_libraries_utf8"=>collect(codeunits(join(Libdl.dllist(),"\n")))))
