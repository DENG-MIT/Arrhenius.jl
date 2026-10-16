# Native half of the timestep-controls growth-strategy rows: the 8 BurnerFlame
# rows of flame_solver_timestep_controls.py (gri30 CH4/air phi=0.8, 300 K,
# 1 atm, mdot 0.09), one solve per (strategy, growth factor) pair, recording
# the accepted time-step count, the largest step, and the converged T_max and
# burner velocity. Prints the same summary row as the Cantera harness
# (tsc_growth_rows_ct.py); Jacobian evaluations are not exposed by the free-
# flame solver, so that column stays empty here.
# Usage: julia --project=<worktree> tsc_growth_rows.jl <mech.yaml> <tag>
using Arrhenius

const ROWS = [
    ("fixed-growth", 1.5),
    ("fixed-growth", 2.0),
    ("steady-norm", 2.0),
    ("transient-residual", 2.0),
    ("residual-ratio", 2.0),
    ("newton-iterations", 2.0),
    ("steady-norm", 6.0),
    ("newton-iterations", 6.0),
]

function build_flame(MECH)
    gas = CreateSolution(MECH)
    # set_equivalence_ratio(0.8, "CH4", {"O2": 1.0, "N2": 3.76})
    return BurnerFlame(gas; mdot=0.09, T=300.0, P=one_atm,
        X="CH4:0.4, O2:1, N2:3.76", width=0.02)
end

mechanism = abspath(ARGS[1])
tag = ARGS[2]
dts_io = open("$(tag)_growth.csv", "w")
println(dts_io, "row,step,dt")

println("row,ok,pts,Tmax,speed,steps,dt_max,wall_s")
for (strategy, growth) in ROWS
    f = build_flame(mechanism)
    steps = Float64[]
    t = @elapsed solve!(f; ratio=3.0, slope=0.25, curve=0.35, prune=0.05,
        time_step_growth_strategy=Symbol(replace(strategy, "-" => "_")),
        time_step_growth_factor=growth,
        time_step_callback=dt -> push!(steps, dt))
    ok = f.converged
    label = "$(strategy) g=$(growth)"
    tmax = ok ? maximum(temperature(f)) : NaN
    speed = ok ? flame_speed(f) : NaN
    println("$(label),$(ok),$(length(f.grid)),$(round(tmax; digits=2)),"
            * "$(round(speed; digits=6)),$(length(steps)),"
            * "$(isempty(steps) ? 0.0 : round(maximum(steps); sigdigits=4)),$(round(t; digits=1))")
    for (i, dt) in enumerate(steps)
        println(dts_io, "$(label),$(i),$(dt)")
    end
    flush(dts_io)
end
close(dts_io)
