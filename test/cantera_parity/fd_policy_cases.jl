# fd_policy_cases.jl — counts + solution dump for NH3 premixed validation cases
# (test/cantera_parity). One Julia session loops over every manifest case with
# kind free/burner/fixed, solving each with the requested discretization and
# recording solver-decision counts (Jacobians, Newton calls/iterations, time
# steps, failures) plus the converged profiles for the accuracy gate.
#
# Usage (from the nh3-parity directory so mechanism sidecars resolve):
#   julia --project=<worktree> fd_policy_cases.jl manifest.csv outdir finite_difference|conservative [case_id ...]
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)

manifest, outdir = ARGS[1], ARGS[2]
scheme = Symbol(ARGS[3])
only = Set(ARGS[4:end])
mkpath(outdir)

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

function build(spec)
    gas = CreateSolution(spec["mechanism"] * ".yaml")
    transport = spec["transport"] == "multicomponent" ? :multicomponent : :mixture_averaged
    data = transport == :multicomponent ?
        MultiTransportData(spec["mechanism"] * ".yaml.multicomponent.npz", gas; mechanism=spec["mechanism"] * ".yaml") : nothing
    X = Dict(String(a) => parse(Float64, b) for (a, b) in split.(split(spec["X"], ","), ":"))
    P = spec["P_atm"] * one_atm
    kw = (; T=spec["T"], P, X, width=spec["width"], transport_model=transport,
        multicomponent_data=data, discretization=scheme)
    spec["kind"] == "free" && return FreeFlame(gas; kw...)
    f = BurnerFlame(gas; mdot=spec["mdot"], kw...)
    spec["kind"] == "fixed" && set_temperature_profile!(f, Float64.(spec["profile_z"]),
        Float64.(spec["profile_T"]); relative=false)
    return f
end

for line in eachline(manifest)
    startswith(line, "id,") && continue
    id, specpath, relpath = split(line, ",")
    spec = YAML.load_file(specpath)
    spec["kind"] in ("free", "burner", "fixed") || continue
    isempty(only) || id in only || continue
    for k in keys(TRACE_COUNTS)
        TRACE_COUNTS[k] = 0
    end
    err = "none"
    t = NaN
    f = build(spec)
    try
        t = @elapsed begin
            kw = haskey(spec, "prune") ? (; prune=min(spec["prune"], spec["slope"], spec["curve"])) : (;)
            solve!(f; ratio=spec["ratio"], slope=spec["slope"], curve=spec["curve"], kw...)
        end
        f.converged || error("did not converge")
    catch e
        err = first(sprint(showerror, e), 200)
    end
    sl = f.converged ? (spec["kind"] == "free" ? flame_speed(f) : velocity(f)[1]) : NaN
    npzwrite(joinpath(outdir, id * ".npz"), Dict(
        "grid" => f.grid, "T" => temperature(f), "Y" => mass_fractions(f),
        "velocity" => velocity(f), "inlet_Y" => copy(f.inlet_Y), "SL" => [sl],
        "species_names_utf8" => join(f.gas.species_names, "\n")))
    counts = Dict(String(k) => v for (k, v) in TRACE_COUNTS)
    open(joinpath(outdir, id * ".json"), "w") do io
        print(io, "{\"id\": \"$id\", \"scheme\": \"$scheme\", \"converged\": $(f.converged), ",
            "\"error\": $(repr(err)), \"seconds\": $t, \"points\": $(length(f.grid)), ",
            "\"SL\": $sl, \"counts\": ")
        print(io, "{\"jacobians\": $(counts["jacobians"]), ",
            "\"newton_calls\": $(counts["newton_calls"]), ",
            "\"newton_iterations\": $(counts["newton_iterations"]), ",
            "\"steps\": $(counts["steps"]), \"failures\": $(counts["failures"])}}")
    end
    println("$id $scheme converged=$(f.converged) J=$(counts["jacobians"]) ",
        "newton=$(counts["newton_iterations"]) steps=$(counts["steps"]) ",
        "fail=$(counts["failures"]) pts=$(length(f.grid)) SL=$sl t=$(round(t, digits=1))s err=$err")
end
