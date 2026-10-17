using Arrhenius, LinearAlgebra, Test

@testset "staged auto solve (Cantera FreeFlame.solve(auto=True))" begin
    mechanism = joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = mole_fractions(gas, "H2:1.1,O2:1,AR:5")
    transport = MultiTransportData(mechanism * ".multicomponent.npz", gas;
        mechanism)
    grid0 = 0.03 .* [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]

    # The staged ladder's fresh analytic guess rebuilds the constructor state
    # and inserts the Cantera Tmid anchor when the adjacent temperatures fall
    # outside the fixed-temperature snap band.
    flame = FreeFlame(gas; T=300.0, P=one_atm, X,
        multicomponent_data=transport, soret=true,
        discretization=:finite_difference, grid=grid0)
    state0 = copy(flame.state)
    Arrhenius._fresh_analytic_guess!(flame, copy(grid0))
    @test length(flame.grid) == length(grid0) + 1
    @test flame.grid ≈ vcat(grid0[1:3], 0.0105, grid0[4:end]) atol=1e-14
    @test flame.state[:, [1, 2, 3, 5, 6, 7, 8, 9]] ≈ state0 rtol=0 atol=1e-12
    @test flame.fixed_temperature == 1000 * flame.state[1, flame.anchor]
    @test !flame.converged

    # Cantera inserts the Tmid crossing on a coarse stage grid when
    # neither adjacent grid temperature falls within its 20% snap band.
    eq = equilibrate(gas; T=300.0, P=one_atm, X, mode=:HP)
    uniform12 = 0.03 .* [0.0, 0.2, 0.3, 0.6, 0.8, 1.0]
    Arrhenius._fresh_analytic_guess!(flame, uniform12)
    @test length(flame.grid) == 7
    @test flame.fixed_temperature ≈ 0.75 * 300.0 + 0.25 * eq.T atol=1e-10
    Arrhenius._fresh_analytic_guess!(flame, uniform12; cantera_compatible=false)
    @test length(flame.grid) == 6

    # The Cantera-compatible initial guess supports burner flames as well.
    burner = BurnerFlame(gas; mdot=0.06, T=300.0, P=one_atm, X,
        discretization=:finite_difference)
    @test Arrhenius._fresh_analytic_guess!(burner, copy(burner.grid)) === burner
    @test burner.state[end, :] == fill(burner.mass_flux, length(burner.grid))
    @test isapprox(burner.state[1, end] * 1000, eq.T; atol=1e-8)

    # H2/air at 300 K, 1 atm: thermal diffusion (Soret) shifts the flame speed
    # by >10%, so the final stage with Soret restored must re-solve. Before
    # the solve-first refinement loop, the final stage returned at the
    # already-converged no-Soret state and both runs agreed bit-for-bit.
    # Cantera 4.0.0a2, commit 726522be4e2a13454d8415b7ef799d621f665cf3,
    # stock FreeFlame.solve(auto=True) on the same mechanism, composition and
    # starting grid: Su = 64.777928779356 cm/s on 148 points (steady-state
    # callback trace), after the same 0.03 -> 0.06 m domain expansion.
    staged(mechanism, transport, X, soret) = solve!(FreeFlame(gas; T=300.0,
        P=one_atm, X, multicomponent_data=transport, soret=soret,
        transport_model=:multicomponent, discretization=:finite_difference, grid=grid0);
        ratio=3.0, slope=0.06, curve=0.10, prune=0.05, max_points=10000,
        loglevel=0, auto=true)
    with_soret = staged(mechanism, transport, X, true)
    @test with_soret.converged
    su_soret = flame_speed(with_soret) * 100
    no_soret_trace = Symbol[]
    no_soret = FreeFlame(gas; T=300.0, P=one_atm, X,
        multicomponent_data=transport, transport_model=:multicomponent,
        discretization=:finite_difference, grid=grid0)
    observe_transport = () -> begin
        model = no_soret.transport_model
        (isempty(no_soret_trace) || last(no_soret_trace) != model) && push!(no_soret_trace, model)
        nothing
    end
    Arrhenius._solve_auto!(no_soret; refine_grid=true, ratio=3.0, slope=0.06,
        curve=0.10, prune=0.05, max_points=10000, max_time_steps=500,
        loglevel=0, domain_check=true, timestep=Ref(1e-5),
        controls=(;growth_factor=1.5, growth_strategy=:fixed_growth,
            step_callback=nothing, analytic=true, steady_callback=observe_transport),
        ts_regrid=3, cantera_compatible=true)
    @test no_soret.converged
    @test no_soret_trace == [:mixture_averaged, :multicomponent]
    @test no_soret.transport_model == :multicomponent
    su_nosoret = flame_speed(no_soret) * 100
    @test solver_stats(no_soret).newton_calls < solver_stats(with_soret).newton_calls
    @test su_soret ≈ 64.777928779356 rtol=2e-4
    @test su_nosoret ≈ 71.557165413082 rtol=1e-4
    @test (su_nosoret - su_soret) / su_soret > 0.1
    # The updated fixed-temperature anchor trajectory lands near the Cantera
    # flame-speed reference while using the same configured refinement rules.
    @test abs(su_soret - 64.777928779356) / 64.777928779356 < 0.02
    # Both runs expanded the too-narrow 0.03 m domain like Cantera.
    @test with_soret.grid[end] - with_soret.grid[1] ≈ 0.06 rtol=1e-12
    @test 62 < length(with_soret.grid) < 200
    @test maximum(abs.(sum(mass_fractions(with_soret); dims=1) .- 1)) < 1e-10
    @test minimum(mass_fractions(with_soret)) > -1e-12
    # The compatibility preset stops on Cantera's weighted Newton correction,
    # so this residual is a bounded diagnostic rather than its stop criterion.
    @test norm(flame_residual!(similar(with_soret.state), with_soret), Inf) < 2e-4
end
