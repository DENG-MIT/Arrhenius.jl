# Capture the opening fixed-temperature Arrhenius decision sequence.
using Arrhenius, YAML, LinearAlgebra
BLAS.set_num_threads(1)
spec = YAML.load_file(ARGS[1]); gas = CreateSolution(spec["mechanism"] * ".yaml")
X = Dict(String(k) => parse(Float64, v) for (k, v) in split.(split(spec["X"], ","), ":"))
f = BurnerFlame(gas; mdot=spec["mdot"], T=spec["T"], P=spec["P_atm"] * one_atm,
    X, width=spec["width"], transport_model=:mixture_averaged, discretization=:finite_difference)
set_temperature_profile!(f, Float64.(spec["profile_z"]), Float64.(spec["profile_T"]); relative=false)
println("INITIAL_GUESS grid=", f.grid, " T_K=", temperature(f), " u=", velocity(f))
const JCOUNT = Ref(0)
const NCALL = Ref(0)
function Arrhenius._flame_jacobian(f::Arrhenius.AbstractPremixedFlame, u, w, r; kwargs...)
    JCOUNT[] += 1
    println("EVENT jacobian n=", JCOUNT[])
    return invoke(Arrhenius._flame_jacobian, Tuple{Any,Any,Any,Any}, f, u, w, r; kwargs...)
end
function Arrhenius._flame_newton!(f::Arrhenius.AbstractPremixedFlame, w; kwargs...)
    NCALL[] += 1
    previous = get(kwargs, :previous, nothing)
    dt = get(kwargs, :dt, Inf)
    println("EVENT newton-enter n=", NCALL[], " transient=", previous !== nothing, " dt=", dt)
    result = invoke(Arrhenius._flame_newton!, Tuple{Any,Any}, f, w; kwargs...)
    iterations = get(kwargs, :iterations, nothing)
    println("EVENT newton-exit n=", NCALL[], " success=", result,
        " iterations=", iterations === nothing ? -1 : iterations[], " jacobians=", JCOUNT[])
    return result
end
try
    solve!(f; ratio=spec["ratio"], slope=spec["slope"], curve=spec["curve"],
        prune=0.0, max_time_steps=10, refine_grid=false, loglevel=2,
        time_step_callback=dt -> (println("EVENT transient-accepted dt=", dt); 0.0))
catch err
    println("STOP ", sprint(showerror, err))
end
stats=Arrhenius.solver_stats(f)
println("FINAL converged=", f.converged, " points=", length(f.grid), " jacobians=", stats.jacobians,
    " newton=", stats.newton_iterations, " accepted_steps=", stats.time_steps)