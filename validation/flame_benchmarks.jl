# Native Julia flame benchmark runner for validation/flame_benchmarks.py.
#
# Usage:
#   julia --project=<repo> validation/flame_benchmarks.jl <mechanism.yaml> <output.npz> \
#       [--reps=N] [--cases=freeflame,burner,burner-fixed]
#
# Runs the two published example configurations (the Cantera 1D premixed-flame
# tutorial cases) plus the optional prescribed-temperature burner case:
#   freeflame:    H2:1.1,O2:1,AR:5 at 300 K, 1 atm, width 0.03 m,
#                 refinement ratio=3, slope=0.06, curve=0.12
#   burner:       H2:1.5,O2:1,AR:7 at 373 K, 0.05 atm, width 0.5 m,
#                 mdot=0.06 kg/(m^2 s), ratio=3, slope=0.05, curve=0.1
#   burner-fixed: same burner with a selected prescribed temperature profile
#
# Every flame is initialized natively (constant-enthalpy equilibrium guess);
# no reference profiles seed the solver. Mechanism/sidecar loading happens once
# before any timing. For each case the first repetition is recorded separately
# as the JIT/compilation run, followed by `reps` warm repetitions. Each timed
# region covers the full calculation: flame construction, initialization, and
# the adaptive solve. Profiles from the final warm repetition are exported for
# the Python orchestrator's correctness comparison. All runtime physics is
# native Julia.

using Arrhenius, NPZ, LinearAlgebra

BLAS.set_num_threads(1)

mechanism = ARGS[1]
output = ARGS[2]
reps = 5
cases = ["freeflame", "burner", "burner-fixed"]
for arg in ARGS[3:end]
    if startswith(arg, "--reps=")
        global reps = parse(Int, split(arg, "=", limit=2)[2])
    elseif startswith(arg, "--cases=")
        global cases = String.(split(split(arg, "=", limit=2)[2], ","))
    else
        error("unknown argument: $arg")
    end
end
reps >= 5 || error("at least 5 warm repetitions are required (got --reps=$reps)")

# Mechanism load and sidecar verification are outside every timed region.
gas = CreateSolution(mechanism)

# A selected prescribed-temperature input for the species boundary-value test.
# It is independent of both solvers' computed solutions.
const FIXED_POSITIONS = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 1.0]
const FIXED_TEMPERATURES = [373.0, 650.0, 1000.0, 1350.0, 1650.0, 1750.0, 1750.0]

function construct_freeflame()
    f = FreeFlame(gas; T=300.0, P=one_atm,
        X=Dict("H2" => 1.1, "O2" => 1.0, "AR" => 5.0), width=0.03)
    solve!(f; ratio=3.0, slope=0.06, curve=0.12)
    return f
end

function construct_burner()
    f = BurnerFlame(gas; T=373.0, P=0.05 * one_atm, mdot=0.06,
        X=Dict("H2" => 1.5, "O2" => 1.0, "AR" => 7.0), width=0.5)
    solve!(f; ratio=3.0, slope=0.05, curve=0.1)
    return f
end

function construct_burner_fixed()
    f = BurnerFlame(gas; T=373.0, P=0.05 * one_atm, mdot=0.06,
        X=Dict("H2" => 1.5, "O2" => 1.0, "AR" => 7.0), width=0.5)
    set_temperature_profile!(f, FIXED_POSITIONS, FIXED_TEMPERATURES)
    solve!(f; ratio=3.0, slope=0.05, curve=0.1)
    return f
end

const CONSTRUCTORS = Dict(
    "freeflame" => construct_freeflame,
    "burner" => construct_burner,
    "burner-fixed" => construct_burner_fixed,
)

results = Dict{String,Any}(
    "species_names_utf8" => Vector{UInt8}(join(gas.species_names, "\n")),
    "system_utf8" => collect(codeunits(Sys.isapple() ? "Darwin" : Sys.islinux() ? "Linux" : string(Sys.KERNEL))),
    "kernel_release_utf8" => collect(codeunits(Sys.isunix() ? readchomp(`uname -r`) : "unknown")),
    "cpu_utf8" => collect(codeunits(Sys.isapple() ? readchomp(`sysctl -n machdep.cpu.brand_string`) : Sys.CPU_NAME)),
)
println("mechanism=$mechanism reps=$reps cases=$(join(cases, ","))")

for case in cases
    construct = CONSTRUCTORS[case]
    seconds = Float64[]
    speeds = Float64[]
    tmax = Float64[]
    points = Int[]
    local f
    # Repetition 1 is the first/JIT run; repetitions 2..reps+1 are warm.
    for rep in 0:reps
        elapsed = @elapsed f = construct()
        f.converged || error("$case did not converge on repetition $rep")
        push!(seconds, elapsed)
        push!(speeds, flame_speed(f))
        push!(tmax, maximum(temperature(f)))
        push!(points, length(f.grid))
        println("$case rep=$rep seconds=", round(elapsed, digits=4),
            " speed=", flame_speed(f), " Tmax=", maximum(temperature(f)),
            " points=", length(f.grid))
    end
    key = replace(case, "-" => "_")
    results[key * "_seconds"] = seconds
    results[key * "_speed"] = speeds
    results[key * "_tmax"] = tmax
    results[key * "_points"] = points
    # Profiles of the final warm repetition for the correctness comparison.
    results[key * "_grid"] = f.grid
    results[key * "_T"] = temperature(f)
    results[key * "_Y"] = mass_fractions(f)
    results[key * "_velocity"] = velocity(f)
    results[key * "_inlet_Y"] = f.inlet_Y
end

npzwrite(output, results)
println("wrote $output")
