# julia --project=EOS_ENV validation/equations_of_state_timing.jl INPUT_DIR OUTPUT.npz
#       [REPETITIONS=9] [informational|controlled|validate-only]
# EOS_ENV contains Arrhenius and Clapeyron 0.6.28.
import_seconds = @elapsed @eval using Arrhenius, Clapeyron, SHA, Dates, Libdl
import_seconds += @elapsed @eval using Arrhenius.NPZ
import_seconds += @elapsed include(joinpath(@__DIR__,"..","example","thermodynamics","equations_of_state.jl"))
include("numerical_threads.jl")

filehash(path) = bytes2hex(sha256(read(path)))

function source_inventory(root)
    files = sort!([replace(relpath(joinpath(dir,name),root),'\\'=>'/')
        for (dir,_,names) in walkdir(root) for name in names])
    return files,[filehash(joinpath(root,file)) for file in files]
end

function co2_native_checks(result)
    result.T == 300.0 || error("source temperature changed")
    length(result.pressure) == 1000 || error("source pressure count changed")
    all(>(0),result.pressure) || error("invalid pressure")
    Set(result.helmholtz.phase) == Set((:vapour,:liquid)) || error("expected both source phases")
    isfinite(result.helmholtz.saturation_pressure) || error("invalid saturation pressure")
    for name in (:ideal,:rk,:helmholtz)
        state = getproperty(result,name)
        p = state.properties
        size(p) == (5,1000) && all(isfinite,p) || error("incomplete/nonfinite property output")
        all(>(0),state.density) && all(>(0),p[4:5,:]) || error("nonphysical density/heat capacity")
        all(iszero,p[1:3,1]) || error("reference-state offsets changed")
        pv = result.pressure ./ state.density ./ 1000
        all(isapprox.(p[1,:].-p[2,:],pv.-pv[1];rtol=1e-8,atol=1e-8)) ||
            error("enthalpy/internal-energy identity failed")
        name == :ideal || all(isapprox.(state.pressure,result.pressure;rtol=1e-8,atol=0.)) ||
            error("pressure closure failed")
    end
    any(path -> occursin(r"(?i)coolprop|cantera|python",path),Libdl.dllist()) &&
        error("foreign fluid calculation library loaded in Julia")
    return true
end

function co2_warm_samples(ideal,rk,helmholtz,first,repetitions,qualification)
    started = time_ns()
    warmup = equations_of_state_calculation(ideal,rk,helmholtz)
    warmup_seconds = (time_ns()-started)/1e9
    isequal(first,warmup) || error("native warmup changed output")
    samples = Float64[]
    batch_size = qualification == "validate-only" ? 0 : clamp(ceil(Int,.020/max(warmup_seconds,1e-9)),1,1000)
    GC.gc()
    if qualification != "validate-only"
        for _ in 1:repetitions
            elapsed = 0.0
            for _ in 1:batch_size
                started = time_ns()
                repeated = equations_of_state_calculation(ideal,rk,helmholtz)
                elapsed += (time_ns()-started)/1e9
                isequal(first,repeated) || error("native repetition changed output")
            end
            push!(samples,elapsed/batch_size)
            benchmark_julia_thread_settings(;enforce=false)
        end
    end
    return samples,batch_size
end

function co2_output_arrays(result)
    data = Dict{String,Any}("T"=>[result.T],"pressure"=>result.pressure,
        "saturation_pressure"=>[result.helmholtz.saturation_pressure],
        "helmholtz_phase"=>[phase==:vapour ? 1 : 2 for phase in result.helmholtz.phase],
        "rk_pressure"=>result.rk.pressure,"helmholtz_pressure"=>result.helmholtz.pressure)
    for name in (:ideal,:rk,:helmholtz)
        state = getproperty(result,name)
        data[String(name)] = state.properties
        data[String(name)*"_density"] = state.density
    end
    return data
end

function main(args)
    length(args) in 2:4 || error("supply INPUT_DIR OUTPUT [REPETITIONS] [QUALIFICATION]")
    directory,destination = args[1:2]
    repetitions = length(args)>=3 ? parse(Int,args[3]) : 9
    qualification = length(args)>=4 ? args[4] : "informational"
    repetitions>=9 || error("at least nine warm batches required")
    qualification in ("informational","controlled","validate-only") || error("invalid qualification")
    before = benchmark_julia_thread_settings()
    core_root = joinpath(pkgdir(Arrhenius),"src")
    inventory_before = source_inventory(core_root)
    paths = (mechanism=joinpath(directory,"co2-thermo.yaml"),
        sidecar=joinpath(directory,"co2-thermo.yaml.npz"),
        parameter=joinpath(directory,"carbon-dioxide.json"),
        thread_helper=joinpath(@__DIR__,"numerical_threads.jl"),
        example=joinpath(@__DIR__,"..","example","thermodynamics","equations_of_state.jl"),
        harness=@__FILE__)
    hashes_before = map(filehash,paths)
    started = time_ns()
    ideal = CreateSolution(paths.mechanism)
    rk = RedlichKwongThermo(paths.mechanism;phase="CO2-RK")
    helmholtz = SingleFluid(read(paths.parameter,String);coolprop_userlocations=false)
    preparation_seconds = (time_ns()-started)/1e9
    # Dynamic outer dispatch prevents compilation of the calculation from being
    # moved ahead of the cold timer by specialization of this wrapper.
    started = time_ns()
    result = Base.invokelatest(equations_of_state_calculation,ideal,rk,helmholtz)
    first_invocation_seconds = (time_ns()-started)/1e9
    data = co2_output_arrays(result)
    data["native_checks_pass"],data["strict_replay_checked"] = [false],[false]
    npzwrite(destination,data) # retain first output if any following check fails
    co2_native_checks(result)
    samples,batch_size = co2_warm_samples(ideal,rk,helmholtz,result,repetitions,qualification)
    after = benchmark_julia_thread_settings(;enforce=false)
    before == after || error("numerical thread settings changed")
    source_inventory(core_root) == inventory_before || error("loaded core source changed")
    map(filehash,paths) == hashes_before || error("example/harness/input changed")
    co2_native_checks(result)
    thread_keys = sort!(collect(keys(after)))
    utf8(value) = collect(codeunits(string(value)))
    cpu = Sys.isapple() ? readchomp(`sysctl -n machdep.cpu.brand_string`) :
        strip(split(first(filter(l -> startswith(l,"model name"),readlines("/proc/cpuinfo"))),':';limit=2)[2])
    merge!(data,Dict(
        "native_checks_pass"=>[true],"strict_replay_checked"=>[true],
        "seconds"=>samples,"batch_size"=>[batch_size],"warm_outputs_checked"=>[length(samples)*batch_size],
        "import_seconds"=>[import_seconds],"preparation_seconds"=>[preparation_seconds],
        "first_invocation_seconds"=>[first_invocation_seconds],
        "numerical_thread_names_utf8"=>utf8(join(thread_keys,"\n")),
        "numerical_threads"=>[after[key] for key in thread_keys],
        "julia_version_utf8"=>utf8(VERSION),"clapeyron_version_utf8"=>utf8(Base.pkgversion(Clapeyron)),
        "qualification_utf8"=>utf8(qualification),"cpu_utf8"=>utf8(cpu),
        "system_utf8"=>utf8(Sys.isapple() ? "Darwin" : "Linux"),
        "kernel_release_utf8"=>utf8(readchomp(`uname -r`)),"timestamp_utc_utf8"=>utf8(Dates.now(Dates.UTC)),
        "harness_sha256_utf8"=>utf8(hashes_before.harness),"example_sha256_utf8"=>utf8(hashes_before.example),
        "parameter_sha256_utf8"=>utf8(hashes_before.parameter),"mechanism_sha256_utf8"=>utf8(hashes_before.mechanism),
        "sidecar_sha256_utf8"=>utf8(hashes_before.sidecar),"thread_helper_sha256_utf8"=>utf8(hashes_before.thread_helper),
        "core_source_paths_utf8"=>utf8(join(inventory_before[1],"\n")),
        "core_source_sha256_utf8"=>utf8(join(inventory_before[2],"\n")),
        "first_call_scope_utf8"=>utf8("Outer invokelatest call includes specialization and complete first calculation; process startup, imports and model preparation excluded.")))
    npzwrite(destination,data)
    println("CO2: ",length(samples)," warm batches of ",batch_size," complete three-EOS sweeps; all outputs checked")
end
main(ARGS)
