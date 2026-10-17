using Arrhenius, Test

@testset "Cantera-compatible free-flame fixed-temperature anchor" begin
    mechanism = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = mole_fractions(gas, "H2:1.1,O2:1,AR:5")
    eq = equilibrate(gas; T=300.0, P=one_atm, X, mode=:HP)
    grid = 0.03 .* [0.0, 0.2, 0.3, 0.6, 0.8, 1.0]
    f = FreeFlame(gas; T=300.0, P=one_atm, X, grid,
        discretization=:finite_difference)

    f.state[end, :] .= 0.42
    Arrhenius._fresh_analytic_guess!(f, grid; cantera_compatible=true)
    @test all(isapprox.(f.state[end, :], 0.42; rtol=4eps()))  # 1-ulp platform (FMA) differences
    @test length(f.grid) == length(grid) + 1
    @test isapprox(f.grid[f.anchor], 0.01125; atol=1e-14)
    @test isapprox(f.fixed_temperature, 0.75 * 300.0 + 0.25 * eq.T; atol=1e-10)
    @test isapprox(f.state[1, f.anchor] * 1000, f.fixed_temperature; atol=1e-10)

    Arrhenius._fresh_analytic_guess!(f, grid; cantera_compatible=false)
    @test f.grid == grid
    @test length(f.state[1, :]) == length(grid)
    @test f.fixed_temperature == 1000 * f.state[1, f.anchor]

    temperature_lock = copy(f.state[1, :]) .+ 0.01
    Arrhenius.FLAME_TEMPERATURE_LOCKS[f] = temperature_lock
    locked_residual = similar(f.state)
    try
        flame_residual!(locked_residual, f)
        @test locked_residual[1, :] == f.state[1, :] .- temperature_lock
    finally
        delete!(Arrhenius.FLAME_TEMPERATURE_LOCKS, f)
    end

    burner = BurnerFlame(gas; mdot=0.06, T=300.0, P=one_atm, X, grid,
        discretization=:finite_difference)
    Arrhenius._fresh_analytic_guess!(burner, grid; cantera_compatible=true)
    @test burner.grid == grid
    @test burner.state[end, :] == fill(burner.mass_flux, length(grid))
    @test isapprox(burner.state[1, end] * 1000, eq.T; atol=1e-8)
    @test !Arrhenius._burner_blown_off(burner)
    flat = copy(burner.state)
    flat[1, :] .= burner.inlet_temperature / 1000
    burner.state .= flat
    @test Arrhenius._burner_blown_off(burner)
    blown = copy(burner.state)
    blown[1, :] .= 300.0 / 1000
    blown[1, 4] += 1e-8
    blown[1, end] = 1500.0 / 1000
    burner.state .= blown
    @test Arrhenius._burner_blown_off(burner)
end
