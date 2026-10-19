#!/usr/bin/env julia
using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)

struct StopAfterFirstSteady <: Exception end
const CAPTURE_PATH = Ref("")
const CAPTURED = Ref(false)

function Arrhenius._flame_steady!(f::Arrhenius.AbstractPremixedFlame; kwargs...)
    success = invoke(Arrhenius._flame_steady!, Tuple{Any}, f; kwargs...)
    if success && !CAPTURED[]
        gas = f.gas
        names = join(gas.species_names, Char(10))
        NPZ.npzwrite(CAPTURE_PATH[], Dict("grid" => copy(f.grid), "T" => Arrhenius.temperature(f),
            "Y" => Arrhenius.mass_fractions(f), "velocity" => Arrhenius.velocity(f),
            "mass_flux" => copy(f.state[end, :]), "species_names_utf8" => Vector{UInt8}(names)))
        CAPTURED[] = true
        throw(StopAfterFirstSteady())
    end
    return success
end

function capture_case(spec_path, output_path)
    spec = YAML.load_file(spec_path)
    gas = CreateSolution(spec["mechanism"] * ".yaml")
    X = Dict(String(k) => parse(Float64, v) for (k, v) in split.(split(spec["X"], ","), ":"))
    flame = BurnerFlame(gas; mdot=spec["mdot"], T=spec["T"], P=spec["P_atm"] * one_atm,
        X, width=spec["width"], discretization=:finite_difference)
    if spec["kind"] == "fixed"
        set_temperature_profile!(flame, Float64.(spec["profile_z"]), Float64.(spec["profile_T"]); relative=false)
    end
    CAPTURE_PATH[] = output_path
    CAPTURED[] = false
    try
        solve!(flame; ratio=spec["ratio"], slope=spec["slope"], curve=spec["curve"],
            prune=min(get(spec, "prune", 0.0), spec["slope"], spec["curve"]),
            refine_grid=false, cantera_compatible=true, auto=(spec["kind"] == "burner"))
    catch err
        err isa StopAfterFirstSteady || rethrow()
    end
    CAPTURED[] || error("no successful steady solve captured for $(spec["id"])")
    println("ARR_FIRST_STEADY case=$(spec["id"]) points=$(length(flame.grid))")
    flush(stdout)
end

length(ARGS) >= 2 && iseven(length(ARGS)) || error("usage: round9_first_arr.jl spec.json output.npz [pairs ...]")
for i in 1:2:length(ARGS)
    capture_case(ARGS[i], ARGS[i+1])
end
