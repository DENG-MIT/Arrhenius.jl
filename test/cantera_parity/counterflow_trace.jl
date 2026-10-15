# Twin premixed counterflow decision trace on the production initial grid.
# TRACE_SPEC: harness spec; TRACE_OUT: directory containing ct_state.npz.
# TRACE_CAP: path to run_capped; PY: Cantera Python executable.
# Run from the mechanism directory under the task wall-time cap.
using Revise
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)
empty!(ARGS)
include(joinpath(@__DIR__, "native_dump_cf.jl"))
spec=YAML.load_file(ENV["TRACE_SPEC"])
spec["kind"]=="twin" || error("trace currently supports twin premixed counterflow")
gas=CreateSolution(spec["mechanism"]*".yaml")
dump_state(ENV["TRACE_OUT"],spec,gas,"ct_state.npz")
run(`$(ENV["TRACE_CAP"]) 60 $(ENV["PY"]) $(joinpath(@__DIR__,"compare_cf.py")) $(ENV["TRACE_OUT"])`)
f=build_native(spec,gas,spec["width"].*[0.,.2,.4,.5,.6,.8,1.])
Arrhenius.reset!(Arrhenius.COUNTERFLOW_STATS)
try
    solve!(f;ratio=spec["ratio"],slope=spec["slope"],curve=spec["curve"],prune=spec["prune"],loglevel=3)
catch err
    showerror(stdout,err,catch_backtrace()); println()
finally
    println("TRACE_FINAL converged=$(f.converged) points=$(length(f.grid)) Tmax=$(maximum(temperature(f)))")
    s=Arrhenius.COUNTERFLOW_STATS
    for k in fieldnames(typeof(s)); println("STATS ",k,"=",getfield(s,k)); end
    npzwrite(joinpath(ENV["TRACE_OUT"],"native_final.npz"),Dict("grid"=>f.grid,"state"=>f.state))
end
