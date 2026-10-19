using Arrhenius, YAML, NPZ, LinearAlgebra
BLAS.set_num_threads(1)

function run_case(spec_path, output_path, variant="compat")
    spec = YAML.load_file(spec_path)
    mechanism = spec["mechanism"] * ".yaml"
    gas = CreateSolution(mechanism)
    transport = spec["transport"] == "multicomponent" ? :multicomponent : :mixture_averaged
    data = transport == :multicomponent ? MultiTransportData(mechanism * ".multicomponent.npz", gas;
        mechanism) : nothing
    X = Dict(String(k) => parse(Float64, v) for (k, v) in split.(split(spec["X"], ","), ":"))
    flame = BurnerFlame(gas; mdot=spec["mdot"], T=spec["T"], P=spec["P_atm"] * one_atm,
        X, width=spec["width"], transport_model=transport, multicomponent_data=data,
        discretization=:finite_difference)
    if spec["kind"] == "fixed"
        set_temperature_profile!(flame, Float64.(spec["profile_z"]), Float64.(spec["profile_T"]); relative=false)
    elseif spec["kind"] != "burner"
        error("unsupported case kind $(spec["kind"])")
    end
    stats = Arrhenius.reset!(Arrhenius.solver_stats(flame))
    prune = min(get(spec, "prune", 0.0), spec["slope"], spec["curve"])
    refine = get(ENV, "ONE_STAGE", "0") != "1"
    if variant == "initial"
        Arrhenius._fresh_analytic_guess!(flame, flame.grid; cantera_compatible=true)
        names = join(gas.species_names, Char(10))
        NPZ.npzwrite(output_path, Dict("grid" => copy(flame.grid),
            "T" => temperature(flame), "Y" => mass_fractions(flame),
            "velocity" => velocity(flame), "species_names_utf8" => Vector{UInt8}(names)))
        println("ARR_INITIAL case=$(spec["id"]) points=$(length(flame.grid))")
        return
    end
    if variant == "default"
        solve!(flame; ratio=spec["ratio"], slope=spec["slope"], curve=spec["curve"], prune, refine_grid=refine, loglevel=0)
    else
        variant in ("compat", "no_matrix_omission", "no_newton_preset", "no_initial_guess_ramp") ||
            error("unknown variant $variant")
        newton = variant == "no_newton_preset" ?
            Arrhenius._newton_policy(nothing, false) : Arrhenius._newton_policy(nothing, true)
        controls = (;growth_factor=1.5, growth_strategy=:fixed_growth, step_callback=nothing,
            analytic=true, cantera_omission=variant != "no_matrix_omission", newton_options=newton)
        timestep = Ref(1e-5)
        guess_compatible = variant != "no_initial_guess_ramp"
        if spec["kind"] == "fixed"
            Arrhenius._fresh_analytic_guess!(flame, flame.grid; cantera_compatible=guess_compatible)
            Arrhenius._auto_refine_loop!(flame; refine_grid=refine, ratio=spec["ratio"],
                slope=spec["slope"], curve=spec["curve"], prune, max_points=1000,
                max_time_steps=500, loglevel=0, domain_check=false, timestep, controls,
                initial_time_step=1e-5, ts_regrid=3, cantera_compatible=true)
        else
            controls = merge(controls, (;steady_callback=()->Arrhenius._burner_blowoff_check!(flame)))
            Arrhenius._auto_staged!(flame; refine_grid=refine, ratio=spec["ratio"],
                slope=spec["slope"], curve=spec["curve"], prune, max_points=1000,
                max_time_steps=500, loglevel=0, domain_check=false, timestep, controls,
                ts_regrid=3, cantera_compatible=true, guess_compatible)
        end
        flame.converged = true
    end
    !refine && (flame.converged = true)
    flame.converged || error("$(spec["id"]) did not converge")
    names = join(gas.species_names, Char(10))
    workspace = Arrhenius.FlameWorkspace(flame)
    residual = similar(flame.state)
    Arrhenius.flame_residual!(residual, flame, flame.state, workspace)
    NPZ.npzwrite(output_path, Dict("grid" => copy(flame.grid), "T" => temperature(flame),
        "Y" => mass_fractions(flame), "inlet_Y" => copy(flame.inlet_Y),
        "velocity" => velocity(flame), "mass_flux" => copy(flame.state[end, :]),
        "species_diffusive_flux" => copy(workspace.flux),
        "species_total_flux_available" => [workspace.conservative !== nothing],
        "species_total_flux" => workspace.conservative === nothing ? zeros(0, 0) :
            copy(workspace.conservative.species_flux),
        "species_names_utf8" => Vector{UInt8}(names),
        "jacobians" => [stats.jacobians], "newton_iterations" => [stats.newton_iterations],
        "time_steps" => [stats.time_steps], "regrids" => [stats.regrids],
        "failures" => [stats.time_step_failures + stats.steady_failures]))
    sumY = maximum(abs.(sum(mass_fractions(flame); dims=1) .- 1))
    mdot_spread = maximum(abs.(flame.state[end, :] .- flame.mass_flux))
    println("ARR variant=$variant case=$(spec["id"]) kind=$(spec["kind"]) points=$(length(flame.grid)) " *
        "jacobians=$(stats.jacobians) newton=$(stats.newton_iterations) steps=$(stats.time_steps) " *
        "regrids=$(stats.regrids) failures=$(stats.time_step_failures + stats.steady_failures) " *
        "sumY_inf=$(sumY) mdot_spread=$(mdot_spread)")
end

length(ARGS) >= 3 && iseven(length(ARGS) - 1) ||
    error("usage: fixed_t_six_counts.jl variant spec.json output.npz [spec.json output.npz ...]")
variant = ARGS[1]
for i in 2:2:length(ARGS)
    run_case(ARGS[i], ARGS[i + 1], variant)
end
