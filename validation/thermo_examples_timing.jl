# Paired native/Cantera thermo-example timing: exports native outputs and
# timings for validation/thermo_examples_timing.py, which independently
# recomputes the same workloads with Cantera 4 and compares. Workloads are the
# complete calculations of Cantera examples thermo/equivalenceRatio.py,
# thermo/isentropic.py, thermo/isentropic_units.py, thermo/sound_speed.py and
# thermo/sound_speed_units.py at pinned source commit
# 726522be4e2a13454d8415b7ef799d621f665cf3, obtained by calling the five
# completed example functions in example/thermodynamics (no native calculation
# is duplicated here).
#
# Usage: julia --project=. validation/thermo_examples_timing.jl OUTPUT.npz
#        [REPETITIONS=9] [informational|controlled|validate-only]
#        --mechanism-dir=PATH
#
# --mechanism-dir must contain the stock Cantera h2o2.yaml (10 species, with
# N2), gri30.yaml and gri30_highT.yaml (53 species each), each paired with its
# .yaml.npz sidecar exported from the same file by mechanism/export_sidecar.py
# (validation/isentropic_cases.py produces exactly this directory). Timed
# regions cover only the complete example calculations on prepared phases;
# mechanism loading, printing and I/O are excluded, and first-call compilation
# is recorded separately from at least nine warm batches. Each sample is the
# mean elapsed calculation time within a batch targeting 20 ms (at most 1000
# complete calls); setup and output checks remain outside every call's timer.
# A controlled
# result additionally requires an otherwise idle target machine.
import_seconds = @elapsed @eval using Arrhenius, NPZ, SHA, Dates, LinearAlgebra
BLAS.set_num_threads(1)
import_seconds += @elapsed @eval include(joinpath(@__DIR__,"..","example","thermodynamics","equivalence_ratio.jl"))
import_seconds += @elapsed @eval include(joinpath(@__DIR__,"..","example","thermodynamics","isentropic_units.jl"))
import_seconds += @elapsed @eval include(joinpath(@__DIR__,"..","example","thermodynamics","sound_speed_units.jl"))

positional = [a for a in ARGS if !startswith(a,"--")]
options = Dict{String,String}()
for a in ARGS
    startswith(a,"--") || continue
    kv = split(a[3:end],'=';limit=2)
    length(kv)==2 || error("expected --key=value, got ",a)
    options[kv[1]] = kv[2]
end
length(positional)>=1 || error("output NPZ path required; see header usage")
output_path = positional[1]
repetitions = length(positional)>=2 ? parse(Int,positional[2]) : 9
qualification = length(positional)>=3 ? positional[3] : "informational"
qualification in ("informational","controlled","validate-only") || error("invalid qualification")
repetitions>=9 || error("at least nine warm repetitions required")
haskey(options,"mechanism-dir") ||
    error("--mechanism-dir=PATH (with h2o2.yaml, gri30.yaml, gri30_highT.yaml and .yaml.npz sidecars) is required")
mechanism_dir = options["mechanism-dir"]
h2o2_path = joinpath(mechanism_dir,"h2o2.yaml")
gri_path = joinpath(mechanism_dir,"gri30.yaml")
hight_path = joinpath(mechanism_dir,"gri30_highT.yaml")
for path in (h2o2_path,gri_path,hight_path)
    isfile(path) || error("missing mechanism ",path)
    isfile(path*".npz") || error("missing native sidecar ",path,".npz (export it with mechanism/export_sidecar.py)")
end

utf8(value) = collect(codeunits(string(value)))
filehash(path) = bytes2hex(sha256(read(path)))
cpu = if Sys.isapple()
    readchomp(`sysctl -n machdep.cpu.brand_string`)
elseif Sys.islinux()
    strip(split(first(filter(l -> startswith(l,"model name"),readlines("/proc/cpuinfo"))),':';limit=2)[2])
else
    Sys.CPU_NAME
end

gas_h2o2 = CreateSolution(h2o2_path)
length(gas_h2o2.species_names)==10 && "N2" in gas_h2o2.species_names ||
    error("the nozzle source uses the stock 10-species h2o2.yaml (with N2); ",
          "the bundled 9-species mechanism/h2o2.yaml is not accepted")
gas_gri = CreateSolution(gri_path)
length(gas_gri.species_names)==53 ||
    error("expected the 53-species GRI-Mech 3.0 mechanism, got ",length(gas_gri.species_names))
gas_hight = CreateSolution(hight_path)
length(gas_hight.species_names)==53 ||
    error("expected the 53-species gri30_highT.yaml mechanism, got ",length(gas_hight.species_names))

data = Dict{String,Any}(
    "cpu_utf8"=>utf8(cpu),"platform_utf8"=>utf8(Sys.MACHINE),"julia_version_utf8"=>utf8(VERSION),
    "system_utf8"=>utf8(Sys.isapple() ? "Darwin" : Sys.islinux() ? "Linux" : string(Sys.KERNEL)),
    "kernel_release_utf8"=>utf8(Sys.isunix() ? readchomp(`uname -r`) : "unknown"),
    "qualification_utf8"=>utf8(qualification),"timestamp_utc_utf8"=>utf8(Dates.now(Dates.UTC)),
    "harness_sha256_utf8"=>utf8(filehash(@__FILE__)),
    "constants_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","Constants.jl"))),
    "solution_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","Solution.jl"))),
    "thermo_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","Thermo","IdealGasThermo.jl"))),
    "thermo_api_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","Thermo.jl"))),
    "equilibrium_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","Equilibrium.jl"))),
    "ideal_gas_states_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","IdealGasStates.jl"))),
    "equivalence_ratio_example_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","example","thermodynamics","equivalence_ratio.jl"))),
    "isentropic_example_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","example","thermodynamics","isentropic.jl"))),
    "isentropic_units_example_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","example","thermodynamics","isentropic_units.jl"))),
    "sound_speed_example_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","example","thermodynamics","sound_speed.jl"))),
    "sound_speed_units_example_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","example","thermodynamics","sound_speed_units.jl"))),
    "h2o2_mechanism_sha256_utf8"=>utf8(filehash(h2o2_path)),
    "gri_mechanism_sha256_utf8"=>utf8(filehash(gri_path)),
    "gri_highT_mechanism_sha256_utf8"=>utf8(filehash(hight_path)),
    "h2o2_sidecar_sha256_utf8"=>utf8(filehash(h2o2_path*".npz")),
    "gri_sidecar_sha256_utf8"=>utf8(filehash(gri_path*".npz")),
    "gri_highT_sidecar_sha256_utf8"=>utf8(filehash(hight_path*".npz")),
    "h2o2_species"=>[length(gas_h2o2.species_names)],
    "gri_species"=>[length(gas_gri.species_names)],
    "gri_highT_species"=>[length(gas_hight.species_names)],
    "julia_threads"=>[Threads.nthreads()],"logical_cpus"=>[Sys.CPU_THREADS],
    "blas_threads"=>[BLAS.get_num_threads()],
    "import_seconds"=>[import_seconds],"repetitions"=>[repetitions],
    "input_description_utf8"=>utf8("equivalence_ratio: 9 mixtures + 16 scalars from thermo/equivalenceRatio.py on 53 GRI species (HP equilibration with source defaults); isentropic: 200-point H2/N2 nozzle from 1200 K/10 atm on stock 10-species h2o2; isentropic_units: the 10-point gri30 variant with explicit K/m^2 conversions; sound_speed: 24 temperatures 300:200:4900 K on 53-species gri30_highT, source rtol=1e-8; sound_speed_units: 25 points 80..4880 degF on gri30, source rtol=1e-6, speeds in ft/s via explicit SI conversions (no Pint)"),
    "reset_note_utf8"=>utf8("Native calculations are stateless: every input (T, P, composition) is passed explicitly to each call, so there is no mutable phase state to reset; the carried-composition behavior of the sound_speed source is part of the example function itself. Every warm output is checked for equality outside the measured region. Native acoustic SP equilibration uses property_rtol=1e-13 (src/IdealGasStates.jl)."),
)

function measured_example(setup,calculation,repetitions,qualification)
    state = setup()
    started = time_ns()
    output = calculation(state)
    first_seconds = (time_ns()-started)/1e9
    samples = Float64[]
    batch_size = 0
    if qualification != "validate-only"
        state = setup()
        started = time_ns()
        warmup = calculation(state)
        warmup_seconds = (time_ns()-started)/1e9
        warmup == output || error("warmup changed the calculation output")
        batch_size = clamp(ceil(Int, 0.020/max(warmup_seconds, 1e-9)), 1, 1000)
        GC.gc()
        for _ in 1:repetitions
            elapsed = 0.0
            for _ in 1:batch_size
                state = setup()
                started = time_ns()
                repeated_output = calculation(state)
                elapsed += (time_ns()-started)/1e9
                repeated_output == output || error("repetition changed the calculation output")
            end
            push!(samples,elapsed/batch_size)
        end
    end
    return output,first_seconds,samples,batch_size
end

function store!(data,prefix,result,first_seconds,samples,batch_size)
    for key in keys(result)
        value = result[key]
        name = prefix*"_"*string(key)
        if value isa AbstractVector{<:AbstractString}
            data[name*"_utf8"] = utf8(join(value,"\n"))
        elseif value isa Number
            data[name] = [Float64(value)]
        else
            data[name] = Array{Float64}(value)
        end
    end
    data[prefix*"_first_seconds"] = [first_seconds]
    data[prefix*"_seconds"] = samples
    data[prefix*"_batch_size"] = [batch_size]
    data[prefix*"_warm_outputs_checked"] = [length(samples)*batch_size]
    println(prefix,": completed; first = ",first_seconds," s; warm batches = ",length(samples),
            "; calculations per batch = ",batch_size)
    return data
end

cases = (
    ("equivalence_ratio",() -> gas_gri,equivalence_ratio_calculation),
    ("isentropic",() -> gas_h2o2,isentropic_calculation),
    ("isentropic_units",() -> gas_gri,isentropic_units_calculation),
    ("sound_speed",() -> gas_hight,sound_speed_calculation),
    ("sound_speed_units",() -> gas_gri,sound_speed_units_calculation),
)
for (prefix,setup,calculation) in cases
    output,first_seconds,samples,batch_size = measured_example(setup,calculation,repetitions,qualification)
    store!(data,prefix,output,first_seconds,samples,batch_size)
end
println("qualification: ",qualification)

npzwrite(output_path,data)
