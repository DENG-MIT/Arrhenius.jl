# Usage: julia --project=. example/flames/diffusion_flame_continuation.jl h2o2.yaml [output_dir]
# Follow stable and unstable H2/O2 counterflow flame branches with two-point
# temperature control. Supply the stock Cantera h2o2.yaml mechanism.

using Arrhenius

include(joinpath(@__DIR__, "..", "..", "validation", "counterflow_continuation.jl"))

function run_diffusion_flame_continuation(gas; capture=true)
    (gas.n_species, gas.n_reactions) == (10, 29) ||
        error("expected stock Cantera h2o2.yaml (10 species, 29 reactions)")
    result = CounterflowContinuation.calculate(gas;
        slope=.0125, curve=.025, prune=.0125,
        initial_slope=.00625, initial_curve=.0125, initial_prune=.00625,
        max_points=5000, capture)
    result.source_success ||
        error("continuation failed: $(result.termination_reason)")
    result
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: julia $(basename(@__FILE__)) h2o2.yaml [output_dir]")
    outdir = length(ARGS) >= 2 ? ARGS[2] : joinpath(pwd(), "diffusion-flame-continuation")
    gas = CreateSolution(ARGS[1])
    result = run_diffusion_flame_continuation(gas)
    mkpath(outdir)
    csv = joinpath(outdir, "continuation.csv")
    open(csv, "w") do io
        names = propertynames(first(result.data))
        println(io, join(names, ","))
        for row in result.data
            println(io, join((getproperty(row, n) for n in names), ","))
        end
    end
    println("maximum strain: $(result.maximum_strain) 1/s")
    println("final strain fraction: $(result.final_ratio)")
    println("wrote $csv")
end
