# julia --project=. validation/condensed_equilibrium_timing.jl GAS.yaml PHASE.json OUTPUT.npz [REPETITIONS=9]
# Complete public 50-point adiabatic calculation on prepared input models.
import_seconds = @elapsed @eval using Arrhenius, NPZ, SHA, TOML, Dates, LinearAlgebra
import_seconds += @elapsed include(joinpath(@__DIR__, "..", "example", "thermodynamics", "adiabatic.jl"))
include("numerical_threads.jl")

file_sha256(path) = bytes2hex(sha256(read(path)))
function source_hashes(root)
    files = [relpath(joinpath(dir, name), root) for (dir, _, names) in walkdir(joinpath(root, "src"))
             for name in names if endswith(name, ".jl")]
    append!(files, ["Project.toml", "example/thermodynamics/adiabatic.jl",
                    "validation/condensed_equilibrium_timing.jl", "validation/numerical_threads.jl"])
    return Dict(replace(file, '\\'=>'/') => file_sha256(joinpath(root, file)) for file in files)
end

function check_output(result)
    length(result.phi) == length(result.T) == 50 || error("expected all 50 equivalence ratios")
    size(result.species_moles) == (54, 50) || error("expected all 54 gas/graphite species amounts")
    all(isfinite, result.phi) && all(isfinite, result.T) && all(isfinite, result.species_moles) ||
        error("nonfinite equilibrium output")
    all(>(0), result.T) && all(>=(0), result.species_moles) || error("nonphysical equilibrium output")
    result.phi ≈ collect(range(.3, 3.5; length=50)) || error("equivalence-ratio sweep differs")
    return nothing
end

function measured_sweep(gas, phase, repetitions)
    started = time_ns()
    result = adiabatic_calculation(gas, phase)
    first_seconds = (time_ns()-started)/1e9
    check_output(result)
    reference = deepcopy(result)
    samples = Float64[]
    GC.gc()
    for _ in 1:repetitions
        started = time_ns()
        repeated = adiabatic_calculation(gas, phase)
        push!(samples, (time_ns()-started)/1e9)
        check_output(repeated)
        isequal(reference, repeated) || error("complete native sweep changed on repetition")
        benchmark_julia_thread_settings(;enforce=false)
    end
    all(t -> isfinite(t) && t > 0, samples) || error("invalid elapsed time")
    return reference, first_seconds, samples
end

function main(args)
    length(args) in (3, 4) || error("supply GAS.yaml PHASE.json OUTPUT.npz [REPETITIONS]")
    gas_path, phase_path, output_path = abspath.(args[1:3])
    repetitions = length(args) == 4 ? parse(Int, args[4]) : 9
    repetitions >= 9 || error("at least nine complete warm sweeps required")
    root = normpath(joinpath(@__DIR__, ".."))
    realpath(pathof(Arrhenius)) == realpath(joinpath(root, "src", "Arrhenius.jl")) ||
        error("loaded a different Arrhenius checkout")
    sources_before = source_hashes(root)
    inputs_before = Dict("gas"=>file_sha256(gas_path), "sidecar"=>file_sha256(gas_path*".npz"),
                         "phase"=>file_sha256(phase_path))
    threads_before = benchmark_julia_thread_settings()
    started = time_ns()
    gas = CreateSolution(gas_path)
    phase = StoichiometricCondensedPhase(phase_path)
    preparation_seconds = (time_ns()-started)/1e9
    gas.n_species == 53 && gas.n_reactions == 325 || error("expected the stock GRI-Mech 3.0 mechanism")
    method = which(adiabatic_calculation, (typeof(gas), typeof(phase)))
    realpath(String(method.file)) == realpath(joinpath(root, "example", "thermodynamics", "adiabatic.jl")) ||
        error("benchmark must call the public example calculation")
    wrapper_started = time_ns()
    result, first_seconds, samples = measured_sweep(gas, phase, repetitions)
    wrapper_seconds = (time_ns()-wrapper_started)/1e9
    threads_after = benchmark_julia_thread_settings(;enforce=false)
    sources_before == source_hashes(root) || error("native source changed during measurement")
    inputs_before == Dict("gas"=>file_sha256(gas_path), "sidecar"=>file_sha256(gas_path*".npz"),
                         "phase"=>file_sha256(phase_path)) || error("input changed during measurement")
    threads_before == threads_after || error("thread settings changed during measurement")
    cpu = Sys.isapple() ? readchomp(`sysctl -n machdep.cpu.brand_string`) :
        strip(split(first(filter(l -> startswith(l, "model name"), readlines("/proc/cpuinfo"))), ':'; limit=2)[2])
    data = Dict("phi"=>result.phi, "T"=>result.T, "species_moles"=>result.species_moles,
                "seconds"=>samples, "first_seconds"=>[first_seconds])
    npzwrite(output_path, data)
    metadata = Dict(
        "source_hashes"=>sources_before, "input_hashes"=>inputs_before,
        "source_unchanged"=>true, "inputs_unchanged"=>true,
        "numerical_threads"=>threads_after, "warm_outputs_checked"=>repetitions,
        "julia_version"=>string(VERSION), "cpu"=>cpu,
        "system"=>(Sys.isapple() ? "Darwin" : "Linux"), "kernel_release"=>readchomp(`uname -r`),
        "import_seconds"=>import_seconds, "preparation_seconds"=>preparation_seconds,
        "measurement_wrapper_seconds"=>wrapper_seconds,
        "first_call_scope"=>"First call inside the measurement wrapper; wrapper specialization, imports and process startup excluded. This is not whole-process cold latency.",
        "timestamp_utc"=>string(Dates.now(Dates.UTC)), "native_artifact_sha256"=>file_sha256(output_path))
    open(output_path*".toml", "w") do stream
        TOML.print(stream, metadata)
    end
    println("Completed first plus ", repetitions, " warm sweeps, each with all 50 equilibrium states")
end

main(ARGS)
