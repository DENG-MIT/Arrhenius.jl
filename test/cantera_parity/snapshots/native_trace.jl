#!/usr/bin/env julia
# Trace the actual FreeFlame pseudo-time sequence and snapshot the first failed attempt.
# Usage from the mechanism directory: julia --project=WORKTREE native_trace.jl spec.json outdir
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)

spec = YAML.load_file(ARGS[1])
outdir = ARGS[2]
mkpath(outdir)
gas = CreateSolution(spec["mechanism"] * ".yaml")
X = Dict(String(k) => parse(Float64, v) for (k, v) in split.(split(spec["X"], ","), ":"))
P = spec["P_atm"] * one_atm
f = FreeFlame(gas; T=spec["T"], P, X, width=spec["width"], transport_model=:mixture_averaged)
f.energy_enabled = true

const ATTEMPTS = NamedTuple[]
const TRANSIENT_SUCCESSES = Ref(0)
const TRANSIENT_FAILURES = Ref(0)
const STEADY_ATTEMPTS = Ref(0)
const STEADY_FAILURES = Ref(0)
const FIRST_FAILURE_SAVED = Ref(false)

function save_first_failure!(f, input_state, dt, iterations, jac_before, jac_after, age_before, age_after)
    FIRST_FAILURE_SAVED[] && return
    names = String.(gas.species_names)
    z = copy(f.grid)
    T = 1000 .* input_state[1, :]
    Y = permutedims(copy(input_state[2:end-1, :]))
    meanMW = 1 ./ (Y' * (1 ./ gas.MW))
    rho = P .* meanMW ./ (Arrhenius.R .* T)
    u = input_state[end, :] ./ rho
    npzwrite(joinpath(outdir, "native_first_failed_input.npz"), Dict(
        "grid" => z, "T" => T, "Y" => Y, "velocity" => u,
        "species_names_utf8" => collect(codeunits(join(names, "\n")))))
    open(joinpath(outdir, "native_first_failure.txt"), "w") do io
        println(io, "attempt_dt=", repr(dt))
        println(io, "newton_iterations=", iterations)
        println(io, "jacobian_before=", jac_before, " jacobian_after=", jac_after)
        println(io, "jacobian_age_before=", age_before, " jacobian_age_after=", age_after)
        println(io, "successful_steps_before=", TRANSIENT_SUCCESSES[])
        println(io, "failed_attempt_index=", TRANSIENT_FAILURES[])
        println(io, "steady_attempts_before=", STEADY_ATTEMPTS[])
        println(io, "steady_failures_before=", STEADY_FAILURES[])
    end
    FIRST_FAILURE_SAVED[] = true
end

function Arrhenius._flame_newton!(flame::FreeFlame, w; kwargs...)
    previous = get(kwargs, :previous, nothing)
    transient = previous !== nothing
    dt = transient ? Float64(get(kwargs, :dt, NaN)) : NaN
    input_state = copy(flame.state)
    jac_before = w.jacobian_evals[]
    age_before = hasproperty(w, :jacobian_age) ? w.jacobian_age[] : -1
    iterations = get(kwargs, :iterations, nothing)
    iterations === nothing && (iterations = Ref(0))
    callkw = merge((; kwargs...), (; iterations))
    ok = invoke(Arrhenius._flame_newton!, Tuple{Any,Any}, flame, w; callkw...)
    jac_after = w.jacobian_evals[]
    age_after = hasproperty(w, :jacobian_age) ? w.jacobian_age[] : -1
    if transient
        ok ? (TRANSIENT_SUCCESSES[] += 1) : (TRANSIENT_FAILURES[] += 1)
    else
        STEADY_ATTEMPTS[] += 1
        ok || (STEADY_FAILURES[] += 1)
    end
    push!(ATTEMPTS, (kind=transient ? "transient" : "steady",
        dt=dt, success=ok, iterations=iterations[], jac_before=jac_before,
        jac_after=jac_after, jac_delta=jac_after-jac_before,
        age_before=age_before, age_after=age_after,
        successful_steps=TRANSIENT_SUCCESSES[], failed_steps=TRANSIENT_FAILURES[],
        steady_attempts=STEADY_ATTEMPTS[], steady_failures=STEADY_FAILURES[]))
    if transient && !ok
        save_first_failure!(flame, input_state, dt, iterations[], jac_before,
            jac_after, age_before, age_after)
    end
    return ok
end

function on_step(dt_next)
    println("CALLBACK successful_steps=", TRANSIENT_SUCCESSES[],
        " dt_next=", repr(dt_next), " failures=", TRANSIENT_FAILURES[],
        " steady_attempts=", STEADY_ATTEMPTS[], " steady_failures=", STEADY_FAILURES[],
        " jacobians=", solver_stats(f).jacobians)
    flush(stdout)
    return 0.0
end

err = "none"
try
    solve!(f; ratio=spec["ratio"], slope=spec["slope"], curve=spec["curve"],
        time_step_callback=on_step, loglevel=0)
catch e
    err = first(sprint(showerror, e), 240)
end
open(joinpath(outdir, "native_attempts.csv"), "w") do io
    println(io, "kind,dt,success,iterations,jac_before,jac_after,jac_delta,age_before,age_after,successful_steps,failed_steps,steady_attempts,steady_failures")
    for a in ATTEMPTS
        println(io, join((a.kind, repr(a.dt), a.success, a.iterations, a.jac_before,
            a.jac_after, a.jac_delta, a.age_before, a.age_after, a.successful_steps,
            a.failed_steps, a.steady_attempts, a.steady_failures), ","))
    end
end
println("NATIVE_TRACE converged=", f.converged, " error=", err,
    " transient_successes=", TRANSIENT_SUCCESSES[],
    " transient_failures=", TRANSIENT_FAILURES[],
    " steady_attempts=", STEADY_ATTEMPTS[], " steady_failures=", STEADY_FAILURES[],
    " attempts=", length(ATTEMPTS), " first_failure_saved=", FIRST_FAILURE_SAVED[])