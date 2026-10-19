using Arrhenius, Test

@testset "failed auto stages restore the pre-attempt flame" begin
    mechanism = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = mole_fractions(gas, "H2:1.1,O2:1,AR:5")
    grid = 0.03 .* [0.0, 0.2, 0.3, 0.6, 0.8, 1.0]
    burner = BurnerFlame(gas; mdot=0.06, T=300.0, P=one_atm, X, grid,
        discretization=:finite_difference)
    original_grid = copy(burner.grid)
    original_state = copy(burner.state)
    original_anchor = burner.anchor
    original_fixed_temperature = burner.fixed_temperature
    original_converged = burner.converged
    attempted_grids = Int[]
    attempt = Ref(0)
    fail_after_mutating!(f; kwargs...) = begin
        attempt[] += 1
        push!(attempted_grids, length(f.grid))
        f.state .+= 0.01
        mod(attempt[], 3) == 2
    end

    @test_throws ErrorException Arrhenius._auto_staged!(burner;
        refine_grid=true, ratio=3.0, slope=0.06, curve=0.10, prune=0.05,
        max_points=48, max_time_steps=0, loglevel=0, domain_check=false,
        timestep=Ref(1e-6), controls=(;), ts_regrid=3,
        steady_solve=fail_after_mutating!)

    @test attempted_grids == [6, 6, 6, 12, 12, 12, 24, 24, 24, 48, 48, 48]
    @test burner.grid == original_grid
    @test burner.state == original_state
    @test burner.anchor == original_anchor
    @test burner.fixed_temperature == original_fixed_temperature
    @test burner.converged == original_converged
    @test isempty(Arrhenius.FLAME_TEMPERATURE_LOCKS)
end
