using Arrhenius, LinearAlgebra, Test

@testset "staged auto solve (Cantera FreeFlame.solve(auto=True))" begin
    mechanism = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = mole_fractions(gas, "H2:1.1,O2:1,AR:5")
    transport = MultiTransportData(mechanism * ".multicomponent.npz", gas;
        mechanism)
    grid0 = 0.03 .* [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]

    # The staged ladder's fresh analytic guess rebuilds the constructor's
    # initial state on the target grid.
    flame = FreeFlame(gas; T=300.0, P=one_atm, X,
        multicomponent_data=transport, soret=true,
        discretization=:finite_difference, grid=grid0)
    state0 = copy(flame.state)
    Arrhenius._fresh_analytic_guess!(flame, copy(grid0))
    @test flame.grid == grid0
    @test flame.state ≈ state0 rtol=0 atol=1e-12
    @test flame.fixed_temperature == 1000 * flame.state[1, flame.anchor]
    @test !flame.converged

    # auto staging mirrors Cantera's FreeFlame.solve(auto=True) and is
    # FreeFlame-only.
    burner = BurnerFlame(gas; mdot=0.06, T=300.0, P=one_atm, X,
        discretization=:finite_difference)
    @test_throws ArgumentError solve!(burner; auto=true)

    # H2/air at 300 K, 1 atm: thermal diffusion (Soret) shifts the flame speed
    # by >10%, so the final stage with Soret restored must re-solve. Before
    # the solve-first refinement loop, the final stage returned at the
    # already-converged no-Soret state and both runs agreed bit-for-bit.
    # Cantera 4.0.0a2, commit 726522be4e2a13454d8415b7ef799d621f665cf3,
    # stock FreeFlame.solve(auto=True) on the same mechanism, composition and
    # starting grid: Su = 62.63803417388627 cm/s on 153 points (steady-state
    # callback trace), after the same 0.03 -> 0.06 m domain expansion.
    staged(mechanism, transport, X, soret) = solve!(FreeFlame(gas; T=300.0,
        P=one_atm, X, multicomponent_data=transport, soret=soret,
        discretization=:finite_difference, grid=grid0); ratio=3.0, slope=0.06,
        curve=0.10, prune=0.05, max_points=10000, loglevel=0, auto=true)
    with_soret = staged(mechanism, transport, X, true)
    @test with_soret.converged
    su_soret = flame_speed(with_soret) * 100
    su_nosoret = flame_speed(staged(mechanism, transport, X, false)) * 100
    @test su_soret ≈ 61.90167652524685 rtol=1e-4
    @test su_nosoret ≈ 70.65494896801594 rtol=1e-4
    @test (su_nosoret - su_soret) / su_soret > 0.1
    # Native lands on 138 points against Cantera's 153 (the guide-grid
    # tolerance sequencing and the per-state-agreed refine/prune rules
    # compound differently along the trajectory); with prune disabled the
    # same solve reaches 151 points and Cantera's value to 0.09%.
    @test abs(su_soret - 62.63803417388627) / 62.63803417388627 < 0.02
    # Both runs expanded the too-narrow 0.03 m domain like Cantera.
    @test with_soret.grid[end] - with_soret.grid[1] ≈ 0.06 rtol=1e-12
    @test 62 < length(with_soret.grid) < 200
    @test maximum(abs.(sum(mass_fractions(with_soret); dims=1) .- 1)) < 1e-10
    @test minimum(mass_fractions(with_soret)) > -1e-12
    @test norm(flame_residual!(similar(with_soret.state), with_soret), Inf) < 1e-6
end
