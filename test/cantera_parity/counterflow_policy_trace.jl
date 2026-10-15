using Arrhenius
using LinearAlgebra
using YAML

BLAS.set_num_threads(1)

initial_time_step = parse(Float64, ARGS[1])
max_time_steps = parse(Int, ARGS[2])

for argument_index in 3:2:length(ARGS)
    global initial_time_step, max_time_steps
    spec = YAML.load_file(ARGS[argument_index])
    kind = ARGS[argument_index + 1]
    gas = CreateSolution(spec["mechanism"] * ".yaml")

    flame = if kind == "premixed"
        CounterflowPremixedFlame(gas; reactants=spec["X"],
            T_reactants=spec["T"], P=spec["P_atm"] * one_atm,
            mdot_reactants=spec["mdot"], mdot_products=spec["mdot_products"],
            width=spec["width"])
    elseif kind == "twin"
        CounterflowTwinPremixedFlame(gas; reactants=spec["X"],
            mdot=spec["mdot"], T=spec["T"], P=spec["P_atm"] * one_atm,
            width=spec["width"])
    else
        error("unsupported kind: $kind")
    end

    Arrhenius.reset!(COUNTERFLOW_STATS)
    converged = Arrhenius._counterflow_steady!(flame;
        max_time_steps, timestep=Ref(initial_time_step), loglevel=2)
    stats = COUNTERFLOW_STATS
    println("native kind=", kind, " converged=", converged,
        " points=", length(flame.grid),
        " Tmax=", round(maximum(temperature(flame)), digits=8),
        " jacobians=", stats.jacobians,
        " newton_iterations=", stats.newton_iterations,
        " time_steps=", stats.time_steps,
        " failures=", stats.time_step_failures)
end
