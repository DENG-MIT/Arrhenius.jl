# Usage: julia --project=. validation/pure_water_timing.jl OUTPUT.npz
#        [REPETITIONS=9] [informational|controlled|validate-only]
# Run the paired Python harness afterward. A controlled result requires an
# otherwise idle target machine; first-call compilation is recorded separately.
import_seconds = @elapsed @eval using Arrhenius, NPZ, SHA, Dates
if !isdefined(Arrhenius,:PureWater)
    import_seconds += @elapsed Base.include(Arrhenius,joinpath(@__DIR__,"..","src","PureWater.jl"))
end
output_path = ARGS[1]
repetitions = length(ARGS)>=2 ? parse(Int,ARGS[2]) : 9
qualification = length(ARGS)>=3 ? ARGS[3] : "informational"
qualification in ("informational","controlled","validate-only") || error("invalid qualification")
repetitions >= 7 || error("at least seven warm repetitions required")

function rankine_calculation(T=300.,P=8e5)
    result = water_rankine(PureWater();inlet_temperature=T,boiler_pressure=P)
    return [result.pump_work,result.turbine_work,result.heat_added,result.efficiency]
end
function vapordome_calculation()
    water = PureWater()
    degc = [WATER_TMIN-273.15;4.;5.;6.;8.;collect(10.:36.);38.;collect(40.:5.:95.);
            collect(100.:10.:290.);collect(300.:20.:360.);collect(370.:373.);WATER_TC-273.15]
    table = zeros(length(degc),14)
    table[:,1] .= degc
    # Match the source's vapor pass, liquid pass, and reference conversion.
    for i in eachindex(degc)
        T = degc[i]+273.15
        state = water_state(water;T,Q=1.)
        table[i,2] = water_saturation(water,T).P/1e5
        table[i,5],table[i,8],table[i,11],table[i,14] = state.v,state.u/1000,state.h/1000,state.s/1000
    end
    for i in eachindex(degc)
        state = water_state(water;T=degc[i]+273.15,Q=0.)
        table[i,3],table[i,6],table[i,9],table[i,12] = state.v,state.u/1000,state.h/1000,state.s/1000
    end
    for (difference,gas,liquid) in ((4,5,3),(7,8,6),(10,11,9),(13,14,12))
        table[:,difference] .= table[:,gas].-table[:,liquid]
    end
    reference = water_state(water;T=WATER_TMIN,Q=0.)
    table[:,[6,8]] .-= reference.u/1000
    table[:,[9,11]] .-= (reference.h-reference.P*reference.v)/1000
    table[:,[12,14]] .-= reference.s/1000
    return table
end

utf8(value) = collect(codeunits(string(value)))
cpu = if Sys.isapple()
    readchomp(`sysctl -n machdep.cpu.brand_string`)
elseif Sys.islinux()
    strip(split(first(filter(l -> startswith(l,"model name"),readlines("/proc/cpuinfo"))),':';limit=2)[2])
else
    Sys.CPU_NAME
end
data = Dict{String,Any}(
    "cpu_utf8"=>utf8(cpu),"platform_utf8"=>utf8(Sys.MACHINE),"julia_version_utf8"=>utf8(VERSION),
    "system_utf8"=>utf8(Sys.isapple() ? "Darwin" : Sys.islinux() ? "Linux" : string(Sys.KERNEL)),
    "kernel_release_utf8"=>utf8(Sys.isunix() ? readchomp(`uname -r`) : "unknown"),
    "qualification_utf8"=>utf8(qualification),"timestamp_utc_utf8"=>utf8(Dates.now(Dates.UTC)),
    "source_sha256_utf8"=>utf8(bytes2hex(sha256(read(joinpath(@__DIR__,"..","src","PureWater.jl"))))),
    "harness_sha256_utf8"=>utf8(bytes2hex(sha256(read(@__FILE__)))),
    "gas_constants_sha256_utf8"=>utf8(bytes2hex(sha256(read(joinpath(@__DIR__,"..","src","Constants.jl"))))),
    "julia_threads"=>[Threads.nthreads()],"logical_cpus"=>[Sys.CPU_THREADS],
    "import_seconds"=>[import_seconds],"repetitions"=>[repetitions],
    "input_description_utf8"=>utf8("rankine:300K,800000Pa,.6,.8; units:80.33F,116.03psi,6894.757293168364Pa/psi; dome:74T,148satstates,1reference"),
)
for (name,run) in (("rankine",rankine_calculation),
        ("rankine_units",()->rankine_calculation((80.33-32)*5/9+273.15,116.03*6894.757293168364)),
        ("vapordome",vapordome_calculation))
    started = time_ns()
    result = run()
    first_seconds = (time_ns()-started)/1e9
    samples = Float64[]
    if qualification != "validate-only"
        run()
        GC.gc()
        for _ in 1:repetitions
            started = time_ns()
            run()
            push!(samples,(time_ns()-started)/1e9)
        end
    end
    data[name*"_output"] = result
    data[name*"_first_seconds"] = [first_seconds]
    data[name*"_seconds"] = samples
    println(name,": completed; first = ",first_seconds," s; warm samples = ",length(samples),"; ",qualification)
end
npzwrite(output_path,data)
