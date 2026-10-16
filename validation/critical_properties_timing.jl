# julia --project=. validation/critical_properties_timing.jl OUTPUT.npz
#       [REPETITIONS=9] [informational|controlled|validate-only]
import_seconds = @elapsed @eval using Arrhenius, NPZ, SHA, Dates, LinearAlgebra
import_seconds += @elapsed include(joinpath(@__DIR__, "..", "example", "thermodynamics", "critical_properties.jl"))
include("numerical_threads.jl")

function measure_critical_properties(fluids, repetitions, qualification)
    started = time_ns()
    output = critical_properties_calculation(fluids)
    first_seconds = (time_ns()-started)/1e9
    samples, batch_size = Float64[], 0
    if qualification != "validate-only"
        started = time_ns()
        warmup = critical_properties_calculation(fluids)
        warmup_seconds = (time_ns()-started)/1e9
        isequal(warmup, output) || error("native warmup changed output")
        batch_size = clamp(ceil(Int, .020/max(warmup_seconds, 1e-9)), 1, 1000)
        GC.gc()
        for _ in 1:repetitions
            elapsed = 0.0
            for _ in 1:batch_size
                started = time_ns()
                repeated = critical_properties_calculation(fluids)
                elapsed += (time_ns()-started)/1e9
                isequal(repeated, output) || error("native repetition changed output")
            end
            push!(samples, elapsed/batch_size)
            benchmark_julia_thread_settings(;enforce=false)
        end
    end
    return output, first_seconds, samples, batch_size
end

function main(args)
    length(args) in 1:3 || error("supply OUTPUT [REPETITIONS] [QUALIFICATION]")
    destination = args[1]
    repetitions = length(args) >= 2 ? parse(Int, args[2]) : 9
    qualification = length(args) >= 3 ? args[3] : "informational"
    repetitions >= 9 || error("at least nine warm batches required")
    qualification in ("informational", "controlled", "validate-only") || error("invalid qualification")
    before = benchmark_julia_thread_settings()
    result, first_seconds, samples, batch_size = measure_critical_properties(CRITICAL_FLUID_NAMES, repetitions, qualification)
    after = benchmark_julia_thread_settings(;enforce=false)
    before == after || error("numerical thread configuration changed")
    utf8(s) = collect(codeunits(string(s)))
    filehash(path) = bytes2hex(sha256(read(path)))
    cpu = if Sys.isapple()
        readchomp(`sysctl -n machdep.cpu.brand_string`)
    else
        strip(split(first(filter(l -> startswith(l, "model name"), readlines("/proc/cpuinfo"))), ':'; limit=2)[2])
    end
    thread_keys = sort!(collect(keys(after)))
    data = Dict{String,Any}(
        "properties" => result, "fluid_names_utf8" => utf8(join(CRITICAL_FLUID_NAMES, "\n")),
        "property_names_utf8" => utf8("critical_temperature\ncritical_pressure\ncritical_density\nmean_molecular_weight\ncritical_compressibility"),
        "first_seconds" => [first_seconds], "seconds" => samples, "batch_size" => [batch_size],
        "warm_outputs_checked" => [length(samples)*batch_size], "import_seconds" => [import_seconds],
        "numerical_thread_names_utf8" => utf8(join(thread_keys, "\n")), "numerical_threads" => [after[key] for key in thread_keys],
        "julia_version_utf8" => utf8(VERSION), "cpu_utf8" => utf8(cpu),
        "system_utf8" => utf8(Sys.isapple() ? "Darwin" : "Linux"), "kernel_release_utf8" => utf8(readchomp(`uname -r`)),
        "qualification_utf8" => utf8(qualification), "timestamp_utc_utf8" => utf8(Dates.now(Dates.UTC)),
        "harness_sha256_utf8" => utf8(filehash(@__FILE__)),
        "first_call_scope_utf8" => utf8("First elapsed call inside the measurement wrapper; wrapper specialization, process startup and imports are excluded. This is not whole-process cold latency."),
    )
    for file in ("src/Constants.jl", "src/CriticalProperties.jl", "example/thermodynamics/critical_properties.jl", "validation/numerical_threads.jl")
        data[file*"_sha256_utf8"] = utf8(filehash(joinpath(@__DIR__, "..", file)))
    end
    npzwrite(destination, data)
    println("critical properties: ", length(samples), " warm batches of ", batch_size, " complete eight-fluid calculations")
end

main(ARGS)
