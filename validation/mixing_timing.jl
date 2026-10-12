# julia --project=. validation/mixing_timing.jl MECHANISM OUTPUT.npz
#       [REPETITIONS=9] [informational|controlled|validate-only]
# Full native example on a prepared phase. Each normalized warm timing sample
# averages complete calculations; result checks remain outside every timer.
import_seconds = @elapsed @eval using Arrhenius, NPZ, SHA, Dates, LinearAlgebra
import_seconds += @elapsed include(joinpath(@__DIR__, "..", "example", "thermodynamics", "mixing.jl"))
BLAS.set_num_threads(1)

function measure_mixing(gas, repetitions, qualification)
    started = time_ns()
    output = mixing_calculation(gas)
    first_seconds = (time_ns()-started)/1e9
    samples, batch_size = Float64[], 0
    if qualification != "validate-only"
        started = time_ns()
        warmup = mixing_calculation(gas)
        warmup_seconds = (time_ns()-started)/1e9
        isequal(warmup, output) || error("native mixing warmup changed output")
        batch_size = clamp(ceil(Int, .020/max(warmup_seconds, 1e-9)), 1, 1000)
        GC.gc()
        for _ in 1:repetitions
            elapsed = 0.0
            for _ in 1:batch_size
                started = time_ns()
                repeated = mixing_calculation(gas)
                elapsed += (time_ns()-started)/1e9
                isequal(repeated, output) || error("native mixing repetition changed output")
            end
            push!(samples, elapsed/batch_size)
        end
    end
    return output, first_seconds, samples, batch_size
end

function main(args)
    length(args) in 2:4 || error("supply MECHANISM OUTPUT [REPETITIONS] [QUALIFICATION]")
    mechanism, destination = args[1:2]
    repetitions = length(args) >= 3 ? parse(Int, args[3]) : 9
    qualification = length(args) >= 4 ? args[4] : "informational"
    repetitions >= 9 || error("at least nine warm batches required")
    qualification in ("informational", "controlled", "validate-only") || error("invalid qualification")
    Threads.nthreads() == 1 || error("launch with JULIA_NUM_THREADS=1")
    gas = CreateSolution(mechanism)
    gas.n_species == 53 || error("the source example uses GRI-Mech with 53 species")
    result, first_seconds, samples, batch_size = measure_mixing(gas, repetitions, qualification)
    states = (result.before, result.after)
    fields = (:T, :P, :density, :meanMW, :h, :u, :s, :g, :cp, :cv, :moles, :mass, :enthalpy)
    labels = ("T", "P", "density", "mean_molecular_weight", "enthalpy_mole", "int_energy_mole",
              "entropy_mole", "gibbs_mole", "cp_mole", "cv_mole", "moles", "mass", "enthalpy")
    utf8(s) = collect(codeunits(string(s)))
    filehash(path) = bytes2hex(sha256(read(path)))
    cpu = if Sys.isapple()
        readchomp(`sysctl -n machdep.cpu.brand_string`)
    else
        strip(split(first(filter(l -> startswith(l, "model name"), readlines("/proc/cpuinfo"))), ':'; limit=2)[2])
    end
    data = Dict{String,Any}(
        "scalars" => [getproperty(state, key) for key in fields, state in states],
        "X" => hcat((state.X for state in states)...), "Y" => hcat((state.Y for state in states)...),
        "mu_RT" => hcat((state.mu_RT for state in states)...),
        "scalar_names_utf8" => utf8(join(labels, "\n")), "species_names_utf8" => utf8(join(gas.species_names, "\n")),
        "first_seconds" => [first_seconds], "seconds" => samples, "batch_size" => [batch_size],
        "warm_outputs_checked" => [length(samples)*batch_size], "import_seconds" => [import_seconds],
        "julia_threads" => [Threads.nthreads()], "blas_threads" => [BLAS.get_num_threads()],
        "julia_version_utf8" => utf8(VERSION), "cpu_utf8" => utf8(cpu),
        "system_utf8" => utf8(Sys.isapple() ? "Darwin" : "Linux"), "kernel_release_utf8" => utf8(readchomp(`uname -r`)),
        "qualification_utf8" => utf8(qualification), "timestamp_utc_utf8" => utf8(Dates.now(Dates.UTC)),
        "harness_sha256_utf8" => utf8(filehash(@__FILE__)),
        "mechanism_sha256_utf8" => utf8(filehash(mechanism)), "sidecar_sha256_utf8" => utf8(filehash(mechanism*".npz")),
        "first_call_scope_utf8" => utf8("First elapsed call inside the measurement wrapper, including only compilation triggered after its timer starts; wrapper specialization, process startup and imports are outside this timer. This is not whole-process cold latency."),
    )
    for file in ("src/Constants.jl", "src/Solution.jl", "src/Thermo.jl", "src/Thermo/IdealGasThermo.jl",
                 "src/Equilibrium.jl", "src/IdealGasMixing.jl", "example/thermodynamics/mixing.jl")
        data[file*"_sha256_utf8"] = utf8(filehash(joinpath(@__DIR__, "..", file)))
    end
    npzwrite(destination, data)
    println("mixing: ", length(samples), " warm batches of ", batch_size, " full calculations; first measured call ", first_seconds, " s")
end

main(ARGS)
