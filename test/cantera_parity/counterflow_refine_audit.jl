# Audit a saved counterflow refinement and final state without solving.
# TRACE_SPEC, TRACE_OUT; BEFORE_REFINE (native NPZ grid/state), CT_AFTER_REFINE
# (ct_dump NPZ), FINAL_STATE (native NPZ grid/state).
using Arrhenius, YAML, NPZ
empty!(ARGS)
include(joinpath(@__DIR__,"native_dump_cf.jl"))
spec=YAML.load_file(ENV["TRACE_SPEC"])
gas=CreateSolution(spec["mechanism"]*".yaml")
before=npzread(ENV["BEFORE_REFINE"])
f=build_native(spec,gas,vec(before["grid"]))
f.state=before["state"]
@assert Arrhenius._refine_axisymmetric!(f;ratio=spec["ratio"],slope=spec["slope"],
    curve=spec["curve"],prune=get(spec,"prune",0.),grid_min=1e-10,max_points=1200,interpolate_velocity=true)
ct=npzread(ENV["CT_AFTER_REFINE"])
@assert f.grid==vec(ct["grid"])
println("REFINE_AUDIT max_velocity_error=",maximum(abs.(velocity(f).-vec(ct["values/velocity"]))))
println("REFINE_AUDIT max_temperature_error=",maximum(abs.(temperature(f).-vec(ct["values/T"]))))
npzwrite(joinpath(ENV["TRACE_OUT"],"refine_corrected.npz"),Dict("grid"=>f.grid,"state"=>f.state,"velocity"=>velocity(f)))
final=npzread(ENV["FINAL_STATE"])
f=build_native(spec,gas,vec(final["grid"]))
f.state=final["state"]
w=Arrhenius.CounterflowWorkspace(f)
r=Arrhenius.counterflow_residual!(similar(f.state),f,f.state,w)
println("INVARIANTS sumY_max=",maximum(abs.(vec(sum(mass_fractions(f);dims=1)).-1)))
println("INVARIANTS continuity_max=",maximum(abs.(@view(r[end,:]))))
println("INVARIANTS residual_max=",maximum(abs.(r)))
npzwrite(joinpath(ENV["TRACE_OUT"],"final_properties.npz"),Dict("grid"=>f.grid,"state"=>f.state,
    "residual"=>r,"source"=>w.properties.source,"flux"=>w.properties.flux,"rho"=>w.properties.rho,
    "MW"=>gas.MW))
