using Arrhenius
using LinearAlgebra
using NPZ
using YAML

BLAS.set_num_threads(1)

summary_path = ARGS[1]
open(summary_path, "w") do summary
    println(summary, "case\tconverged\tpoints\tjacobians\tnewton_iterations\t",
        "time_steps\ttime_step_failures\tsteady_failures\tTmax\t",
        "sumY_max\tminY\tresidual_norm")
    global argument_index = 2
    while argument_index <= length(ARGS)
        spec = YAML.load_file(ARGS[argument_index])
        output_path = ARGS[argument_index + 1]
        flame_spec = spec
        gas = CreateSolution(flame_spec["mechanism"] * ".yaml")
        pressure = flame_spec["P_atm"] * one_atm
        kind = flame_spec["kind"]
        flame = if kind == "diffusion"
            CounterflowDiffusionFlame(gas; fuel=flame_spec["fuel"],
                oxidizer=flame_spec["oxidizer"], mdot_fuel=flame_spec["mdot_fuel"],
                mdot_oxidizer=flame_spec["mdot_oxidizer"], T_fuel=flame_spec["T"],
                T_oxidizer=flame_spec["T"], P=pressure, width=flame_spec["width"])
        elseif kind == "premixed"
            CounterflowPremixedFlame(gas; reactants=flame_spec["X"],
                T_reactants=flame_spec["T"], P=pressure,
                mdot_reactants=flame_spec["mdot"],
                mdot_products=flame_spec["mdot_products"], width=flame_spec["width"])
        elseif kind == "twin"
            CounterflowTwinPremixedFlame(gas; reactants=flame_spec["X"],
                mdot=flame_spec["mdot"], T=flame_spec["T"], P=pressure,
                width=flame_spec["width"])
        elseif kind == "stagnation"
            ImpingingJet(gas; reactants=flame_spec["X"], mdot=flame_spec["mdot"],
                T_inlet=flame_spec["T"], T_surface=flame_spec["T_surface"],
                P=pressure, width=flame_spec["width"])
        else
            error("unsupported counterflow kind: " * kind)
        end
        Arrhenius.reset!(COUNTERFLOW_STATS)
        solve!(flame; ratio=flame_spec["ratio"], slope=flame_spec["slope"],
            curve=flame_spec["curve"],
            prune=min(get(flame_spec, "prune", 0.0), flame_spec["slope"], flame_spec["curve"]))
        flame.converged || error(flame_spec["id"] * " did not converge")
        species = mass_fractions(flame)
        residual = similar(flame.state)
        counterflow_residual!(residual, flame)
        stats = COUNTERFLOW_STATS
        values = (flame_spec["id"], flame.converged, length(flame.grid),
            stats.jacobians, stats.newton_iterations, stats.time_steps,
            stats.time_step_failures, stats.steady_failures,
            maximum(temperature(flame)), maximum(abs.(sum(species; dims=1))),
            minimum(species), norm(residual))
        println(summary, join(values, "\t"))
        flush(summary)
        npzwrite(output_path, Dict("batch_means" => [0.0], "per_batch" => [1],
            "grid" => copy(flame.grid), "T" => temperature(flame),
            "Y" => species, "velocity" => velocity(flame),
            "inlet_Y" => zeros(gas.n_species),
            "SL" => [maximum(temperature(flame))],
            "species_names_utf8" => Vector{UInt8}(join(gas.species_names, "\n"))))
        global argument_index += 2
    end
end
