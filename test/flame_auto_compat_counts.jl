using Arrhenius, YAML, LinearAlgebra

BLAS.set_num_threads(1)

function make_flame(spec, spec_path)
    spec_root = dirname(dirname(abspath(spec_path)))
    mechanism = joinpath(spec_root, spec["mechanism"] * ".yaml")
    gas = CreateSolution(mechanism)
    transport = spec["transport"] == "multicomponent" ? :multicomponent : :mixture_averaged
    data = transport == :multicomponent ?
        MultiTransportData(mechanism * ".multicomponent.npz", gas; mechanism) : nothing
    X = Dict(String(k) => parse(Float64, v)
        for (k, v) in split.(split(spec["X"], ","), ":"))
    common = (; T=spec["T"], P=spec["P_atm"] * one_atm, X,
        width=spec["width"], transport_model=transport, multicomponent_data=data,
        discretization=:finite_difference)
    if spec["kind"] == "free"
        return FreeFlame(gas; common...)
    elseif spec["kind"] in ("burner", "fixed")
        flame = BurnerFlame(gas; mdot=spec["mdot"], common...)
        if spec["kind"] == "fixed"
            set_temperature_profile!(flame, Float64.(spec["profile_z"]),
                Float64.(spec["profile_T"]); relative=false)
        end
        return flame
    end
    error("only free, burner, and fixed specs are supported: $(spec["kind"])")
end

dry_run = any(arg -> strip(arg) == "--dry-run", ARGS)
spec_paths = filter(arg -> strip(arg) != "--dry-run", ARGS)
isempty(spec_paths) && error("provide one or more flame spec paths")
for path in spec_paths
    spec = YAML.load_file(path)
    modes = if spec["kind"] == "fixed"
        (false,) # Cantera-compatible auto staging rejects prescribed-temperature burners.
    else
        dry_run ? (true,) : (false, true)
    end
    for compatible in modes
        flame = make_flame(spec, path)
        stats = solver_stats(flame)
        solve!(flame; ratio=spec["ratio"], slope=spec["slope"],
            curve=spec["curve"], prune=min(get(spec, "prune", 0.0),
            spec["slope"], spec["curve"]), cantera_compatible=compatible,
            refine_grid=!dry_run, max_time_steps=dry_run ? 40 : 500, loglevel=0)
        flame.converged || error("$(spec["id"]) did not converge")
        value = spec["kind"] == "free" ? 100 * flame_speed(flame) : 100 * velocity(flame)[1]
        Y = mass_fractions(flame)
        sum_y_error = maximum(abs.(sum(Y; dims=1) .- 1))
        atomic_weights = flame.gas.ele_matrix' \ flame.gas.MW
        element_mass_matrix = Diagonal(atomic_weights) * flame.gas.ele_matrix *
            Diagonal(1 ./ flame.gas.MW)
        element_in = element_mass_matrix * flame.inlet_Y
        element_out = element_mass_matrix * @view(Y[:, end])
        element_endpoint_delta = maximum(abs.(element_out .- element_in)) /
            max(maximum(abs.(element_in)), 1e-300)
        mass_flux_error = maximum(abs.(flame.state[end, :] .- flame.state[end, 1])) /
            max(abs(flame.state[end, 1]), 1e-300)
        expected_ladder = get(Dict(
            "stagni_2023__free__phi1.2_h20.3_p1__mixture-averaged" => (9, 283),
            "alzueta_2023__free__phi1.2_h20.3_p1__mixture-averaged" => (9, 305)),
            spec["id"], nothing)
        if compatible && !dry_run && expected_ladder !== nothing
            expected_regrids, expected_points = expected_ladder
            @assert stats.regrids == expected_regrids
            @assert length(flame.grid) == expected_points
            @assert isapprox(flame.grid[end] - flame.grid[1], spec["width"]; atol=1e-14, rtol=0)
        end
        println("case=$(spec["id"]) mode=$(compatible ? "compat" : "plain") " *
            "jacobians=$(stats.jacobians) newton_calls=$(stats.newton_calls) " *
            "newton_iterations=$(stats.newton_iterations) steady_calls=$(length(stats.grids)) " *
            "grid_trace=$(join(stats.grids, ",")) steps=$(stats.time_steps) " *
            "steady_failures=$(stats.steady_failures) " *
            "time_step_failures=$(stats.time_step_failures) regrids=$(stats.regrids) " *
            "points=$(length(flame.grid)) speed_cm_s=$(value) sumY_error=$(sum_y_error) " *
            "element_endpoint_delta=$(element_endpoint_delta) mass_flux_error=$(mass_flux_error)")
    end
end
