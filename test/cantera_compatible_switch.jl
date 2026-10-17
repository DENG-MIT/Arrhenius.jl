using Arrhenius, Test

const MECHANISM = joinpath(dirname(pathof(Arrhenius)), "..", "mechanism", "h2o2.yaml")
newgas() = CreateSolution(MECHANISM)

const NEWTON_FIELDS = (:damping_trials, :damping_factor, :damping_acceptance,
    :step_bounds, :norm_denominator, :flow_weight, :steady_convergence,
    :premixed_steady_maxiters, :counterflow_steady_maxiters,
    :premixed_counterflow_steady_maxiters, :premixed_transient_maxiters,
    :premixed_policy_transient_maxiters, :counterflow_transient_maxiters,
    :jacobian_age, :jacobian_reevaluations)
const CURRENT_NEWTON = (damping_trials=7, damping_factor=sqrt(2.0),
    damping_acceptance=:current, step_bounds=:current, norm_denominator=:current,
    flow_weight=:current, steady_convergence=:residual_inf,
    premixed_steady_maxiters=35, counterflow_steady_maxiters=45,
    premixed_counterflow_steady_maxiters=0, premixed_transient_maxiters=20,
    premixed_policy_transient_maxiters=200, counterflow_transient_maxiters=100,
    jacobian_age=20, jacobian_reevaluations=4)
const CANTERA_NEWTON = (damping_trials=7, damping_factor=sqrt(2.0),
    damping_acceptance=:step_or_reduction, step_bounds=:cantera,
    norm_denominator=:cantera, flow_weight=:velocity,
    steady_convergence=:weighted_step, premixed_steady_maxiters=0,
    counterflow_steady_maxiters=0, premixed_counterflow_steady_maxiters=0,
    premixed_transient_maxiters=0, premixed_policy_transient_maxiters=0,
    counterflow_transient_maxiters=0, jacobian_age=20, jacobian_reevaluations=4)

function check_solver_options(flame, compatible, expected)
    compatible ? solve!(flame; cantera_compatible=true) : solve!(flame)
    @test flame.converged
    settings = solver_options(flame)
    @test settings !== nothing
    @test settings.jacobian_mode == expected.jacobian_mode
    @test settings.staged_ladder == expected.staged_ladder
    @test settings.auto == expected.auto
    @test settings.newton_options.cantera_compatible == compatible
    newton_expected = compatible ? CANTERA_NEWTON : CURRENT_NEWTON
    @test all(getproperty(settings.newton_options, name) == getproperty(newton_expected, name)
        for name in NEWTON_FIELDS)
    @test all(getproperty(settings, name) == value for (name, value) in pairs(expected.defaults))
end

function check_solver_options_reset(flame, expected)
    solve!(flame; jacobian_mode=:finite_difference)
    @test flame.converged
    first = solver_options(flame)
    @test first !== nothing
    @test first.jacobian_mode == :finite_difference

    Arrhenius.reset!(solver_stats(flame))
    @test solver_options(flame) === nothing

    solve!(flame; jacobian_mode=:auto, newton_options=NewtonPolicyOptions(damping_trials=5))
    @test flame.converged
    latest = solver_options(flame)
    @test latest !== nothing
    @test latest.jacobian_mode == :auto
    @test latest.newton_options.damping_trials == 5
    @test latest !== first
    @test all(getproperty(latest, name) == value for (name, value) in pairs(expected.defaults))
end

@testset "cantera_compatible solve switch integration" begin
    X = "H2:2,O2:1,AR:3"
    premixed_defaults = (refine_grid=true, ratio=3.0, slope=.06, curve=.12,
        prune=0.0, max_points=1000, max_time_steps=500, loglevel=0,
        initial_time_step=1e-5, domain_check=true, time_step_growth_factor=1.5,
        time_step_growth_strategy=:fixed_growth, time_step_callback=nothing, time_step_regrid=3)
    premixed_expected(compatible, kind) = (jacobian_mode=compatible ? :auto : :finite_difference,
        staged_ladder=compatible, auto=false,
        defaults=merge(premixed_defaults, (domain_check=kind === FreeFlame,)))
    for kind in (FreeFlame, BurnerFlame), compatible in (false, true)
        gas = newgas()
        flame = kind === FreeFlame ? FreeFlame(gas; X,
            grid=collect(range(0.0, .01; length=17))) :
            BurnerFlame(gas; X, mdot=.1, grid=collect(range(0.0, .01; length=6)))
        check_solver_options(flame, compatible, premixed_expected(compatible, kind))
    end

    cf_defaults = (refine_grid=true, ratio=4.0, slope=.2, curve=.3, prune=0.0,
        max_points=1200, max_time_steps=800, loglevel=0,
        initial_time_step=Arrhenius._COUNTERFLOW_FD_INITIAL_TIME_STEP,
        time_step_growth_factor=1.5, time_step_growth_strategy=:fixed_growth,
        time_step_callback=nothing, time_step_regrid=0)
    cf_expected() = (jacobian_mode=:auto, staged_ladder=false, auto=true, defaults=cf_defaults)
    for compatible in (false, true)
        flame = CounterflowDiffusionFlame(newgas(); fuel="H2:1,AR:1",
            oxidizer="O2:.2,AR:.8", mdot_fuel=.24, mdot_oxidizer=.72,
            P=one_atm, width=.02, grid=collect(range(0.0,.02;length=5)))
        check_solver_options(flame, compatible, cf_expected())
    end

    pcf_defaults = (refine_grid=true, ratio=3.0, slope=.1, curve=.2, prune=0.0,
        grid_min=1e-10, max_points=1200, max_time_steps=800, loglevel=0,
        initial_time_step=Arrhenius._COUNTERFLOW_FD_INITIAL_TIME_STEP,
        time_step_growth_factor=1.5, time_step_growth_strategy=:fixed_growth,
        time_step_callback=nothing, time_step_regrid=0, refinement_interpolation=:velocity)
    pcf_expected() = (jacobian_mode=:auto, staged_ladder=false, auto=true, defaults=pcf_defaults)
    composition = "H2:1.6,O2:1,AR:7"
    for kind in (CounterflowPremixedFlame, CounterflowTwinPremixedFlame, ImpingingJet),
            compatible in (false, true)
        gas = newgas()
        flame = if kind === CounterflowPremixedFlame
            kind(gas; reactants=composition, mdot_reactants=.12,
                mdot_products=.06, T_reactants=373.0, P=.05*one_atm, width=.2)
        elseif kind === CounterflowTwinPremixedFlame
            kind(gas; reactants=composition, mdot=.12, T=373.0, P=.05*one_atm, width=.2)
        else
            kind(gas; reactants=composition, mdot=.06, T_inlet=373.0,
                T_surface=500.0, P=.05*one_atm, width=.2)
        end
        check_solver_options(flame, compatible, pcf_expected())
    end

    @testset "counterflow solver options reset and refresh" begin
        diffusion = CounterflowDiffusionFlame(newgas(); fuel="H2:1,AR:1",
            oxidizer="O2:.2,AR:.8", mdot_fuel=.24, mdot_oxidizer=.72,
            P=one_atm, width=.02, grid=collect(range(0.0,.02;length=5)))
        check_solver_options_reset(diffusion, cf_expected())

        composition = "H2:1.6,O2:1,AR:7"
        premixed = CounterflowPremixedFlame(newgas(); reactants=composition,
            mdot_reactants=.12, mdot_products=.06, T_reactants=373.0,
            P=.05*one_atm, width=.2)
        twin = CounterflowTwinPremixedFlame(newgas(); reactants=composition,
            mdot=.12, T=373.0, P=.05*one_atm, width=.2)
        jet = ImpingingJet(newgas(); reactants=composition, mdot=.06,
            T_inlet=373.0, T_surface=500.0, P=.05*one_atm, width=.2)
        for flame in (premixed, twin, jet)
            check_solver_options_reset(flame, pcf_expected())
        end
    end
end
