# julia --project=. validation/coverage_thermo_timing.jl PARAMETER_DIRECTORY OUTPUT.npz
#       [REPETITIONS=9] [informational|controlled|validate-only]
# Complete coverage sweep on prepared models: four 101-state enthalpy/entropy
# curves plus the 5151-state triangular cross-interaction map (5555 states).
# Model preparation is timed separately; result checks stay outside every timer.
import_seconds = @elapsed @eval using Arrhenius, NPZ, SHA, Dates, LinearAlgebra
import_seconds += @elapsed include(joinpath(@__DIR__, "..", "example", "thermodynamics", "coverage_dependent_surf.jl"))
include("numerical_threads.jl")

function verified_archives(directory)
    archives = Dict{String,Any}()
    source_sha256 = nothing
    for name in COVERAGE_PHASE_NAMES
        path = joinpath(directory, name*".coverage.json")
        isfile(path) || error("missing coverage parameter archive: $path")
        data = Arrhenius.YAML.load_file(path)
        data["format"] == "arrhenius-coverage-thermo-v1" ||
            error("unsupported coverage archive format in $path")
        data["phase_name"] == name ||
            error("archive $path holds phase $(data["phase_name"]), expected $name")
        if source_sha256 === nothing
            source_sha256 = data["source_sha256"]
        else
            data["source_sha256"] == source_sha256 ||
                error("coverage archives were exported from different mechanisms")
        end
        archives[name] = path
    end
    return archives, source_sha256
end

function measure_coverage_sweep(prepared, repetitions, qualification)
    started = time_ns()
    output = coverage_dependent_surf_calculation(prepared)
    first_seconds = (time_ns()-started)/1e9
    samples, batch_size = Float64[], 0
    if qualification != "validate-only"
        started = time_ns()
        warmup = coverage_dependent_surf_calculation(prepared)
        warmup_seconds = (time_ns()-started)/1e9
        isequal(warmup, output) || error("native warmup changed output")
        batch_size = clamp(ceil(Int, .020/max(warmup_seconds, 1e-9)), 1, 1000)
        GC.gc()
        for _ in 1:repetitions
            elapsed = 0.0
            for _ in 1:batch_size
                started = time_ns()
                repeated = coverage_dependent_surf_calculation(prepared)
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
    length(args) in 2:4 || error("supply PARAMETER_DIRECTORY OUTPUT [REPETITIONS] [QUALIFICATION]")
    directory, destination = args[1:2]
    repetitions = length(args) >= 3 ? parse(Int, args[3]) : 9
    qualification = length(args) >= 4 ? args[4] : "informational"
    repetitions >= 9 || error("at least nine warm batches required")
    qualification in ("informational", "controlled", "validate-only") || error("invalid qualification")
    archives, source_sha256 = verified_archives(directory)
    before = benchmark_julia_thread_settings()
    started = time_ns()
    prepared = prepare_coverage_sweep(directory)
    prepare_seconds = (time_ns()-started)/1e9
    result, first_seconds, samples, batch_size = measure_coverage_sweep(prepared, repetitions, qualification)
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
        "coverages" => result.coverages, "curves" => result.curves, "cross" => result.cross,
        "phase_names_utf8" => utf8(join(COVERAGE_PHASE_NAMES, "\n")),
        "first_seconds" => [first_seconds], "prepare_seconds" => [prepare_seconds],
        "seconds" => samples, "batch_size" => [batch_size],
        "warm_outputs_checked" => [length(samples)*batch_size], "import_seconds" => [import_seconds],
        "numerical_thread_names_utf8" => utf8(join(thread_keys, "\n")), "numerical_threads" => [after[key] for key in thread_keys],
        "julia_version_utf8" => utf8(VERSION), "cpu_utf8" => utf8(cpu),
        "system_utf8" => utf8(Sys.isapple() ? "Darwin" : "Linux"), "kernel_release_utf8" => utf8(readchomp(`uname -r`)),
        "qualification_utf8" => utf8(qualification), "timestamp_utc_utf8" => utf8(Dates.now(Dates.UTC)),
        "mechanism_source_sha256_utf8" => utf8(source_sha256),
        "harness_sha256_utf8" => utf8(filehash(@__FILE__)),
        "first_call_scope_utf8" => utf8("First elapsed call inside the measurement wrapper; model preparation, wrapper specialization, process startup and imports are excluded. This is not whole-process cold latency."),
    )
    for name in COVERAGE_PHASE_NAMES
        data["archive_"*name*"_sha256_utf8"] = utf8(filehash(archives[name]))
    end
    for file in ("src/Constants.jl", "src/SurfaceKinetics.jl", "src/CoverageDependentThermo.jl",
                 "example/thermodynamics/coverage_dependent_surf.jl", "validation/numerical_threads.jl")
        data[file*"_sha256_utf8"] = utf8(filehash(joinpath(@__DIR__, "..", file)))
    end
    npzwrite(destination, data)
    println("coverage sweep: ", length(samples), " warm batches of ", batch_size,
            " complete 5555-state sweeps; first measured call ", first_seconds, " s")
end

main(ARGS)
