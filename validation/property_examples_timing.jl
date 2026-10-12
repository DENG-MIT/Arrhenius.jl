# Paired native/Cantera property-example timing: exports native outputs and
# timings for validation/property_examples_timing.py, which independently
# recomputes the same workloads with Cantera 4 and compares. Workloads are the
# complete calculations of Cantera examples kinetics/blowers_masel.py and
# transport/dusty_gas.py at pinned source commit
# 726522be4e2a13454d8415b7ef799d621f665cf3.
#
# Usage: julia --project=. validation/property_examples_timing.jl OUTPUT.npz
#        [REPETITIONS=9] [informational|controlled|validate-only]
#        [--gri-mechanism=PATH] --h2o2-mechanism=PATH --h2o2-sidecar=PATH
#
# --h2o2-mechanism must be the stock 10-species h2o2.yaml (with N2) used by the
# dusty_gas.py source, not the bundled 9-species mechanism/h2o2.yaml, and
# --h2o2-sidecar must be exported from that same file with
# mechanism/export_multicomponent.py. Timed regions cover only the complete
# example calculations on prepared phases/workspaces; mechanism loading,
# workspace construction, printing, plotting, and I/O are excluded, and
# first-call compilation is recorded separately from at least nine warm
# repetitions. A controlled result additionally requires an otherwise idle
# target machine.
import_seconds = @elapsed @eval using Arrhenius, NPZ, SHA, Dates, LinearAlgebra
BLAS.set_num_threads(1)
import_seconds += @elapsed @eval include(joinpath(@__DIR__,"..","example","kinetics","blowers_masel.jl"))
import_seconds += @elapsed @eval include(joinpath(@__DIR__,"..","example","transport","dusty_gas.jl"))

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
gri_path = get(options,"gri-mechanism",joinpath(@__DIR__,"..","mechanism","gri30.yaml"))
haskey(options,"h2o2-mechanism") ||
    error("--h2o2-mechanism=PATH (stock 10-species h2o2.yaml with N2) is required")
haskey(options,"h2o2-sidecar") ||
    error("--h2o2-sidecar=PATH (multicomponent sidecar exported from the same h2o2.yaml) is required")
h2o2_path = options["h2o2-mechanism"]
sidecar_path = options["h2o2-sidecar"]

utf8(value) = collect(codeunits(string(value)))
filehash(path) = bytes2hex(sha256(read(path)))
cpu = if Sys.isapple()
    readchomp(`sysctl -n machdep.cpu.brand_string`)
elseif Sys.islinux()
    strip(split(first(filter(l -> startswith(l,"model name"),readlines("/proc/cpuinfo"))),':';limit=2)[2])
else
    Sys.CPU_NAME
end

gas_bm = CreateSolution(gri_path)
length(gas_bm.species_names)==53 ||
    error("expected the 53-species GRI-Mech 3.0 mechanism, got ",length(gas_bm.species_names))
gas_dg = CreateSolution(h2o2_path)
length(gas_dg.species_names)==10 && "N2" in gas_dg.species_names ||
    error("dusty_gas.py uses the stock 10-species h2o2.yaml (with N2); ",
          "the bundled 9-species mechanism/h2o2.yaml is not accepted")
X_dg = mole_fractions(gas_dg,Dict("OH"=>1.,"H"=>2.,"O2"=>3.,"O"=>1e-8,"H2"=>1e-8,
    "H2O"=>1e-8,"H2O2"=>1e-8,"HO2"=>1e-8,"AR"=>1e-8))

data = Dict{String,Any}(
    "cpu_utf8"=>utf8(cpu),"platform_utf8"=>utf8(Sys.MACHINE),"julia_version_utf8"=>utf8(VERSION),
    "system_utf8"=>utf8(Sys.isapple() ? "Darwin" : Sys.islinux() ? "Linux" : string(Sys.KERNEL)),
    "kernel_release_utf8"=>utf8(Sys.isunix() ? readchomp(`uname -r`) : "unknown"),
    "qualification_utf8"=>utf8(qualification),"timestamp_utc_utf8"=>utf8(Dates.now(Dates.UTC)),
    "harness_sha256_utf8"=>utf8(filehash(@__FILE__)),
    "constants_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","Constants.jl"))),
    "blowers_masel_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","BlowersMasel.jl"))),
    "dusty_gas_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","DustyGasTransport.jl"))),
    "multicomponent_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","MulticomponentTransport.jl"))),
    "transport_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","Transport.jl"))),
    "thermo_module_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","src","Thermo","IdealGasThermo.jl"))),
    "blowers_masel_example_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","example","kinetics","blowers_masel.jl"))),
    "dusty_gas_example_sha256_utf8"=>utf8(filehash(joinpath(@__DIR__,"..","example","transport","dusty_gas.jl"))),
    "gri_mechanism_sha256_utf8"=>utf8(filehash(gri_path)),
    "h2o2_mechanism_sha256_utf8"=>utf8(filehash(h2o2_path)),
    "h2o2_sidecar_sha256_utf8"=>utf8(filehash(sidecar_path)),
    "gri_species"=>[length(gas_bm.species_names)],"h2o2_species"=>[length(gas_dg.species_names)],
    "julia_threads"=>[Threads.nthreads()],"logical_cpus"=>[Sys.CPU_THREADS],
    "blas_threads"=>[BLAS.get_num_threads()],
    "import_seconds"=>[import_seconds],"repetitions"=>[repetitions],
    "input_description_utf8"=>utf8("blowers_masel: 3 example reactions (Arrhenius + 2 BlowersMasel, A=38.7 b=2.7 Ea0=26191840 J/kmol W=1e9), 53 GRI species, T=300:100:3400, 100 enthalpy shifts within +/-5x the intrinsic barrier; dusty_gas: stock 10-species h2o2, T=500 K, P=1 and 1.2 atm, porosity .2, tortuosity 4, pore radius 1.5e-7 m, particle diameter 1.5e-6 m, delta 1e-3 m"),
    "reset_note_utf8"=>utf8("blowers_masel_calculation uses a stateless Solution and evaluates the prescribed reaction-enthalpy shift directly (physically equivalent to Cantera's H-thermo rewrite); dusty_gas_calculation runs on a freshly constructed workspace per repetition so transport caches stay cold, and every returned array is copied out of workspace buffers. Every warm output is checked outside the measured region."),
)

function measured_example(setup,calculation,repetitions,qualification)
    state = setup()
    started = time_ns()
    output = calculation(state)
    first_seconds = (time_ns()-started)/1e9
    samples = Float64[]
    if qualification != "validate-only"
        calculation(setup()) == output || error("warmup changed the calculation output")
        GC.gc()
        for _ in 1:repetitions
            state = setup()
            started = time_ns()
            repeated_output = calculation(state)
            elapsed = (time_ns()-started)/1e9
            repeated_output == output || error("repetition changed the calculation output")
            push!(samples,elapsed)
        end
    end
    return output,first_seconds,samples
end

bm,bm_first,bm_samples = measured_example(() -> gas_bm,blowers_masel_calculation,repetitions,qualification)
data["blowers_masel_temperatures"] = copy(bm.temperatures)
data["blowers_masel_rates"] = copy(bm.rates)
data["blowers_masel_enthalpies"] = copy(bm.enthalpies)
data["blowers_masel_barriers"] = copy(bm.barriers)
data["blowers_masel_first_seconds"] = [bm_first]
data["blowers_masel_seconds"] = bm_samples
println("blowers_masel: completed; first = ",bm_first," s; warm samples = ",length(bm_samples),"; ",qualification)

# Dusty gas: fresh workspace per repetition; dusty_gas_calculation copies the
# diffusion matrix out of the workspace and returns fresh flux vectors.
multi_data = MultiTransportData(sidecar_path,gas_dg)
make_dusty_workspace() = DustyGasTransport(gas_dg;porosity=.2,tortuosity=4.,
    mean_pore_radius=1.5e-7,mean_particle_diameter=1.5e-6,
    multicomponent_data=multi_data)
dg,dg_first,dg_samples = measured_example(make_dusty_workspace,
    w -> dusty_gas_calculation(w,X_dg),repetitions,qualification)
data["dusty_gas_diffusion"] = copy(dg.diffusion)
data["dusty_gas_conductivity"] = [dg.conductivity]
data["dusty_gas_zero_flux"] = copy(dg.uniform)
data["dusty_gas_pressure_flux"] = copy(dg.pressure_gradient)
data["dusty_gas_first_seconds"] = [dg_first]
data["dusty_gas_seconds"] = dg_samples
println("dusty_gas: completed; first = ",dg_first," s; warm samples = ",length(dg_samples),"; ",qualification)

npzwrite(output_path,data)
