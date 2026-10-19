using Arrhenius, Test

@testset "Cantera-compatible burner auto solve" begin
    mechanism = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = mole_fractions(gas, "H2:1.1,O2:1,AR:5")
    grid = 0.03 .* [0.0, 0.2, 0.3, 0.6, 0.8, 1.0]
    burner = BurnerFlame(gas; mdot=0.06, T=300.0, P=one_atm, X, grid,
        discretization=:finite_difference)
    solve!(burner; cantera_compatible=true, refine_grid=false)
    @test burner.converged
    sum_y_error = maximum(abs.(sum(burner.state[2:end-1, :]; dims=1) .- 1))
    mass_flux_error = maximum(abs.(burner.state[end, :] ./ burner.mass_flux .- 1))
    @info "burner invariants" sum_y_error mass_flux_error
    @test sum_y_error < 1e-10
    @test mass_flux_error < 1e-12
end

@testset "Cantera-compatible burner velocity guess" begin
    mechanism = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")
    grid = 0.03 .* [0.0, 0.05, 0.1, 0.15, 0.2, 0.4, 1.0]
    for fixed_temperature in (false, true)
        gas = CreateSolution(mechanism)
        burner = BurnerFlame(gas; mdot=0.06, T=300.0, P=one_atm,
            X="H2:1.1,O2:1,AR:5", grid, discretization=:finite_difference)
        fixed_temperature && set_temperature_profile!(burner,
            [0.0, 0.2, 1.0], [300.0, 1500.0, 700.0])
        Arrhenius._fresh_analytic_guess!(burner, grid; cantera_compatible=true)
        u = velocity(burner)
        initial_T = temperature(burner)
        x = burner.inlet_Y ./ gas.MW
        x ./= sum(x)
        eq = equilibrate(gas; T=burner.inlet_temperature, P=burner.pressure, X=x, mode=:HP)
        u0 = burner.mass_flux / Arrhenius._flame_density(burner,
            burner.inlet_temperature, burner.inlet_Y).density
        u1 = burner.mass_flux / Arrhenius._flame_density(burner, eq.T, eq.Y).density
        for j in eachindex(grid)
            fraction = clamp((grid[j] - grid[1]) / (0.2 * (grid[end] - grid[1])), 0.0, 1.0)
            @test initial_T[j] ≈ burner.inlet_temperature + fraction * (eq.T - burner.inlet_temperature) rtol=1e-12
            @test u[j] ≈ u0 + fraction * (u1 - u0) rtol=1e-10
        end
        if fixed_temperature
            @test maximum(abs.(initial_T .- burner.imposed_temperature)) > 100
        end
        compatible_state = copy(burner.state)
        Arrhenius._fresh_analytic_guess!(burner, grid; cantera_compatible=false)
        @test maximum(abs.(burner.state[end, :] .- burner.mass_flux)) < 1e-14
        @test maximum(abs.(compatible_state[end, :] .- burner.state[end, :])) > 1e-6
    end
end


@testset "Cantera-compatible fixed-temperature burner solve" begin
    gas = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))
    grid = collect(range(0.0, 0.03; length=5))
    burner = BurnerFlame(gas; mdot=0.07, T=300.0, P=one_atm,
        X="H2:1.1,O2:1,AR:5", grid, discretization=:finite_difference)
    set_temperature_profile!(burner, [0.0, 0.25, 1.0], [300.0, 1600.0, 300.0])
    prescribed = copy(burner.imposed_temperature)
    Arrhenius._fresh_analytic_guess!(burner, burner.grid; cantera_compatible=true)
    initial_T = copy(temperature(burner))
    initial_velocity = copy(velocity(burner))
    x = burner.inlet_Y ./ gas.MW
    x ./= sum(x)
    eq = equilibrate(gas; T=burner.inlet_temperature, P=burner.pressure, X=x, mode=:HP)
    u0 = burner.mass_flux / Arrhenius._flame_density(burner,
        burner.inlet_temperature, burner.inlet_Y).density
    u1 = burner.mass_flux / Arrhenius._flame_density(burner, eq.T, eq.Y).density
    for j in eachindex(burner.grid)
        fraction = clamp((burner.grid[j] - burner.grid[1]) /
            (0.2 * (burner.grid[end] - burner.grid[1])), 0.0, 1.0)
        @test initial_T[j] ≈ burner.inlet_temperature + fraction * (eq.T - burner.inlet_temperature) rtol=1e-12
        @test initial_velocity[j] ≈ u0 + fraction * (u1 - u0) rtol=1e-10
    end
    @test maximum(abs.(initial_T .- prescribed)) > 100
    solve!(burner; cantera_compatible=true, refine_grid=false)
    @test burner.converged
    @test temperature(burner) ≈ prescribed
    @test maximum(abs.(sum(mass_fractions(burner); dims=1) .- 1)) < 1e-10
    @test solver_stats(burner).jacobians > 0
end


@testset "staged burner guess compatibility control" begin
    mechanism = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")
    grid = 0.03 .* [0.0, 0.05, 0.1, 0.2, 0.4, 1.0]
    controls = (;growth_factor=1.5, growth_strategy=:fixed_growth, step_callback=nothing,
        analytic=true, cantera_omission=true,
        newton_options=Arrhenius._newton_policy(nothing, true))

    function stop_after_guess!(guess_compatible)
        gas = CreateSolution(mechanism)
        flame = BurnerFlame(gas; mdot=0.06, T=300.0, P=one_atm,
            X="H2:1.1,O2:1,AR:5", grid, discretization=:finite_difference)
        stop_solve = (f; kwargs...) -> error("stop after staged guess")
        @test_throws ErrorException Arrhenius._auto_staged!(flame; refine_grid=false,
            ratio=3.0, slope=0.06, curve=0.12, prune=0.0, max_points=100,
            max_time_steps=1, loglevel=0, domain_check=false, timestep=Ref(1e-5),
            controls, ts_regrid=3, cantera_compatible=true, guess_compatible,
            steady_solve=stop_solve)
        return copy(flame.state)
    end

    expected_compat = BurnerFlame(CreateSolution(mechanism); mdot=0.06, T=300.0,
        P=one_atm, X="H2:1.1,O2:1,AR:5", grid, discretization=:finite_difference)
    Arrhenius._fresh_analytic_guess!(expected_compat, grid; cantera_compatible=true)
    expected_legacy = BurnerFlame(CreateSolution(mechanism); mdot=0.06, T=300.0,
        P=one_atm, X="H2:1.1,O2:1,AR:5", grid, discretization=:finite_difference)
    Arrhenius._fresh_analytic_guess!(expected_legacy, grid; cantera_compatible=false)

    @test stop_after_guess!(true) ≈ expected_compat.state rtol=1e-14 atol=1e-14
    @test stop_after_guess!(false) ≈ expected_legacy.state rtol=1e-14 atol=1e-14
end
