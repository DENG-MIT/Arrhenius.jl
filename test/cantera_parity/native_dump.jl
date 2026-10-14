# native_dump.jl — native-side dump tool mirroring ct_dump.py (test/cantera_parity).
#
# Usage (from the nh3-parity directory so mechanism sidecars resolve):
#   julia --project=<worktree> native_dump.jl spec.json outdir state <ct_state.npz>
#       Loads the Cantera state (grid, T, u, Y) into native flames with the same
#       grid and dumps per-row residuals, Jacobian-derived Newton step, and the
#       bound fraction for discretization=:conservative and :finite_difference:
#           outdir/native_state_cons.npz, outdir/native_state_fd.npz
#   julia --project=<worktree> native_dump.jl spec.json outdir trace
#       Full solve! with FLAMEDBG decision lines on stdout (redirect to capture)
#       and per-timestep states in outdir/native_trace.npz
#
# State exchange: native rows are [T/1000; Y_1..Y_n; mdot].  mdot is rebuilt
# from the Cantera velocity as rho(T,Y,P)*u, exactly as set_initial_guess! does.
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)

spec = YAML.load_file(ARGS[1]); outdir = ARGS[2]; mode = ARGS[3]
gas = CreateSolution(spec["mechanism"] * ".yaml")
X = Dict(String(a) => parse(Float64, b) for (a, b) in split.(split(spec["X"], ","), ":"))
P = spec["P_atm"] * one_atm
mkpath(outdir)

# Solver-decision counts, fd_full_v6.jl style: invoke-wrappers around the
# production Jacobian and Newton entry points.
const TRACE_COUNTS = Dict(:jacobians => 0, :newton_calls => 0,
    :newton_iterations => 0, :steps => 0, :failures => 0)

function Arrhenius._flame_jacobian(f::Arrhenius.AbstractPremixedFlame, u, w, r; kwargs...)
    TRACE_COUNTS[:jacobians] += 1
    return invoke(Arrhenius._flame_jacobian, Tuple{Any,Any,Any,Any}, f, u, w, r; kwargs...)
end

function Arrhenius._flame_newton!(f::Arrhenius.AbstractPremixedFlame, w; kwargs...)
    TRACE_COUNTS[:newton_calls] += 1
    result = invoke(Arrhenius._flame_newton!, Tuple{Any,Any}, f, w; kwargs...)
    iterations = get(kwargs, :iterations, nothing)
    iterations === nothing || (TRACE_COUNTS[:newton_iterations] += iterations[])
    if get(kwargs, :previous, nothing) !== nothing
        result ? (TRACE_COUNTS[:steps] += 1) : (TRACE_COUNTS[:failures] += 1)
    end
    return result
end

const steps = NamedTuple[]

function build_native(scheme, grid)
    kw = (; T=spec["T"], P, X, width=spec["width"], transport_model=:mixture_averaged,
        discretization=scheme, grid)
    spec["kind"] == "free" && return FreeFlame(gas; kw...)
    return BurnerFlame(gas; mdot=spec["mdot"], kw...)
end

"Overwrite the native state with the Cantera solution at identical points."
function load_ct_state!(f, ct, meta)
    z = vec(ct["grid"]); T = vec(ct["values/T"]); u = vec(ct["values/velocity"])
    names = String.(meta["species_names"])
    names == gas.species_names || error("species order mismatch: $names vs $(gas.species_names)")
    n = gas.n_species; N = length(z)
    Y = reduce(hcat, [vec(ct["values/$s"]) for s in names])
    meanMW = 1 ./ (Y * (1 ./ gas.MW))
    rho = P .* meanMW ./ (Arrhenius.R .* T)
    state = zeros(n + 2, N)
    state[1, :] = T ./ 1000
    state[2:n+1, :] = permutedims(Y)
    state[end, :] = rho .* u
    f.grid = z
    f.state = state
    # Reselect the fixed-temperature anchor exactly as the constructor does,
    # then pin the fixed temperature at the anchor value as Cantera does.
    f.anchor = clamp(argmin(abs.(state[1, :] .- (f.fixed_temperature / 1000))), 2, N - 1)
    f.fixed_temperature = 1000 * state[1, f.anchor]
    return f
end

function flame_jacobian_step(f, w, r)
    u = f.state; B = size(u, 1); bandwidth = 2B - 1
    J = Arrhenius._flame_jacobian(f, u, w, r; analytic=true)
    _, pivots = LinearAlgebra.LAPACK.gbtrf!(bandwidth, bandwidth, length(u), J)
    correction = -vec(copy(r))
    LinearAlgebra.LAPACK.gbtrs!('N', bandwidth, bandwidth, length(u), w.band, pivots, correction)
    return reshape(correction, size(u))
end

function bound_fraction(f, step)
    u = f.state; B = size(u, 1)
    alpha = Inf
    for j in axes(u, 2), k in axes(u, 1)
        low, high = Arrhenius._flame_bounds(f, k, B)
        if step[k, j] < 0
            alpha = min(alpha, (u[k, j] - low) / (-step[k, j]))
        elseif step[k, j] > 0
            alpha = min(alpha, (high - u[k, j]) / step[k, j])
        end
    end
    return alpha
end

function dump_state(scheme)
    ct = npzread(joinpath(outdir, "ct_state.npz"))
    meta = YAML.load_file(joinpath(outdir, "ct_meta.json"))
    f = build_native(scheme, vec(ct["grid"]))
    load_ct_state!(f, ct, meta)
    w = Arrhenius.FlameWorkspace(f)
    r = Arrhenius.flame_residual!(similar(f.state), f, f.state, w)
    step = flame_jacobian_step(f, w, r)
    alpha = bound_fraction(f, step)
    tag = scheme == :conservative ? "cons" : "fd"
    npzwrite(joinpath(outdir, "native_state_$tag.npz"), Dict(
        "grid" => f.grid, "state" => f.state, "residual" => r,
        "newton_step" => step, "fbound" => [alpha], "anchor" => [f.anchor]))
    open(joinpath(outdir, "native_meta_$tag.json"), "w") do io
        names = join(("\"" * n * "\"" for n in ["T"; gas.species_names; "M"]), ", ")
        write(io, "{\"row_names\": [$names]}")
    end
    println("native_dump state $scheme: N=$(length(f.grid)) fbound=$alpha")
end

function dump_trace()
    scheme = length(ARGS) > 3 ? Symbol(ARGS[4]) : :conservative
    f = build_native(scheme, spec["width"] .* [0, .2, .3, .35, .4, .5, .6, .8, 1])
    Arrhenius.DEBUG_FLAME[] = true
    cb = dt -> push!(steps, (dt=dt, grid=copy(f.grid), state=copy(f.state),
        jac=TRACE_COUNTS[:jacobians], iters=TRACE_COUNTS[:newton_iterations]))
    err = "none"
    try
        solve!(f; ratio=spec["ratio"], slope=spec["slope"], curve=spec["curve"],
            loglevel=2, time_step_callback=cb)
    catch e
        err = first(sprint(showerror, e), 120)
    finally
        Arrhenius.DEBUG_FLAME[] = false
    end
    out = Dict{String,Any}("converged" => [Float64(f.converged)],
        "n_steps" => [length(steps)],
        "jacobians" => [TRACE_COUNTS[:jacobians]],
        "newton_calls" => [TRACE_COUNTS[:newton_calls]],
        "newton_iterations" => [TRACE_COUNTS[:newton_iterations]],
        "ts_failures" => [TRACE_COUNTS[:failures]],
        "step_dt" => [s.dt for s in steps],
        "step_jac" => [s.jac for s in steps],
        "step_iters" => [s.iters for s in steps])
    for (i, s) in enumerate(steps)
        out["step_$(i-1)/grid"] = s.grid
        out["step_$(i-1)/state"] = s.state
    end
    tag = scheme == :conservative ? "" : "_fd"
    npzwrite(joinpath(outdir, "native_trace$tag.npz"), out)
    println("native_dump trace: converged=$(f.converged) err=$err steps=$(length(steps)) ",
        "jac=$(TRACE_COUNTS[:jacobians]) newton=$(TRACE_COUNTS[:newton_iterations])")
end

if mode == "state"
    dump_state(:conservative)
    dump_state(:finite_difference)
elseif mode == "trace"
    dump_trace()
else
    error("unknown mode $mode")
end
