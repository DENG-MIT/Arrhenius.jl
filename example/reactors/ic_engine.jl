# Eight revolutions of the published n-dodecane ideal-gas engine calculation.
# Prepare dodecane_IG.yaml + .npz with validation/ic_engine_case.py once.
# julia --project=<SciML environment> example/reactors/ic_engine.jl <mechanism.yaml> [output.csv]
# Solver dependencies: SciMLBase, OrdinaryDiffEqBDF and ForwardDiff.
# CSV profiles use the temperature/crank-angle limits. The companion TOML
# reports integrals from all accepted states and the pressure-work ledger.
include("ic_engine_qndf_solver.jl")
using .NativeEngineQNDF
include("ic_engine_setup.jl")
using TOML
if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("provide the prepared dodecane_IG.yaml mechanism")
    path = get(ARGS,2,"ic_engine.csv")
    calculation = NativeEngineQNDF.engine_cli_failure_guard(path*".diagnostics") do
        solve_ic_engine_qndf(ARGS[1];progress=true)
    end
    result,output = calculation.result,calculation.output
    path = write_ic_engine_csv(path,output,result.gas)
    summary = calculation.summary
    open(path*".toml","w") do io
        TOML.print(io,summary;sorted=true)
    end
    println(result.integrals)
    println("Integrated ",result.integration_points," accepted states; exported ",length(result.times)," profile states.")
end
