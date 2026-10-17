using Arrhenius, Test

@testset "Newton policy option resolution" begin
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    X=mole_fractions(gas,"H2:2,O2:1,AR:3")
    f=FreeFlame(gas;X,grid=collect(range(0.,1.;length=5)),discretization=:finite_difference)
    current=Arrhenius._newton_policy(nothing,false)
    compatible=Arrhenius._newton_policy(nothing,true)

    @test current.cantera_compatible == false
    @test (current.damping_trials,current.damping_factor,current.damping_acceptance)==(7,sqrt(2.),:current)
    @test (current.step_bounds,current.norm_denominator,current.flow_weight,current.steady_convergence)==
        (:current,:current,:current,:residual_inf)
    @test (current.premixed_steady_maxiters,current.counterflow_steady_maxiters,
        current.premixed_counterflow_steady_maxiters,current.premixed_transient_maxiters,
        current.premixed_policy_transient_maxiters,current.counterflow_transient_maxiters)==(35,45,0,20,200,100)
    @test current.jacobian_age==20 && current.jacobian_reevaluations==4

    expected_cantera=(damping_trials=7,damping_factor=sqrt(2.),
        damping_acceptance=:step_or_reduction,step_bounds=:cantera,
        norm_denominator=:cantera,flow_weight=:velocity,
        steady_convergence=:weighted_step,premixed_steady_maxiters=0,
        counterflow_steady_maxiters=0,premixed_counterflow_steady_maxiters=0,
        premixed_transient_maxiters=0,premixed_policy_transient_maxiters=0,
        counterflow_transient_maxiters=0,jacobian_age=20,jacobian_reevaluations=4)
    @test compatible.cantera_compatible
    @test all(getproperty(compatible,k)==v for (k,v) in pairs(expected_cantera))

    # Every field supplied by the caller overlays the compatible preset alone.
    overrides=(damping_trials=3,damping_factor=2.0,
        damping_acceptance=:correction_reduction,step_bounds=:current,
        norm_denominator=:current,flow_weight=:current,
        steady_convergence=:residual_inf,premixed_steady_maxiters=1,
        counterflow_steady_maxiters=2,premixed_counterflow_steady_maxiters=3,
        premixed_transient_maxiters=4,premixed_policy_transient_maxiters=5,
        counterflow_transient_maxiters=6,jacobian_age=9,jacobian_reevaluations=2)
    for (field,value) in pairs(overrides)
        option=NewtonPolicyOptions(cantera_compatible=true; (field=>value,)...)
        resolved=Arrhenius._newton_policy(option,false)
        @test resolved.cantera_compatible
        @test getproperty(resolved,field)==value
        for (other,preset_value) in pairs(expected_cantera)
            other==field || @test getproperty(resolved,other)==preset_value
        end
    end
    @test Arrhenius._newton_policy(NewtonPolicyOptions(cantera_compatible=true,damping_trials=3),false).damping_trials==3
    @test Arrhenius._newton_policy(NewtonPolicyOptions(damping_trials=3),true).damping_trials==3

    @test !Arrhenius._newton_merit_accepts(current,0.95,0.9,false)
    @test Arrhenius._newton_merit_accepts(compatible,0.95,0.9,false)
    reduction=NewtonPolicyOptions(damping_acceptance=:correction_reduction)
    @test !Arrhenius._newton_merit_accepts(Arrhenius._newton_policy(reduction,false),0.95,0.9,true)
    # Construct the reported edge: all seven damping trials reject, although
    # the last trial merit is <= 1. Neither solver may infer success from it.
    rejected_trial_sequence=[(false,1.4),(false,1.3),(false,1.2),(false,1.1),(false,1.0),(false,0.9),(false,0.75)]
    @test all(!accepted for (accepted,_) in rejected_trial_sequence)
    @test last(rejected_trial_sequence)[2] <= 1
    @test !Arrhenius._damping_failure_fallback(last(rejected_trial_sequence)...,true)
    @test !Arrhenius._damping_failure_fallback(false,0.75,true)
    @test Arrhenius._damping_failure_fallback(false,1.01,true)
    @test Arrhenius._cantera_damped_step_converged(true,1.0)
    @test Arrhenius._cantera_damped_step_converged(true,0.9)
    @test !Arrhenius._cantera_damped_step_converged(false,0.9)
    @test !Arrhenius._cantera_damped_step_converged(true,1.000001)

    # A five-point startup can fail to converge, but its bound-rejection path
    # must return the solver's ordinary failure instead of reading an undefined
    # correction merit (the regression that arose while preparing these traces).
    bound_failure_flame=FreeFlame(gas;X,grid=collect(range(0.,1.;length=5)),
        discretization=:finite_difference)
    @test_throws ErrorException solve!(bound_failure_flame;refine_grid=false,auto=false,
        max_time_steps=100,cantera_compatible=true)

    # Cantera 4.0.0a2 FreeFlame bounds probe (Flow1D::setBounds):
    #                    Cantera T       Arrhenius T (stored in kK)
    # h2o2.yaml          (200, 7000)      (0.2, 7.0)
    # glarborg_2018      (200, 4000)      (0.2, 4.0)
    # creck_nh3          (200, 7000)      (0.2, 7.0)
    # Every probe also returned Y=(-1e-7,1e5), velocity=(-1e20,1e20),
    # and Lambda=(-1e20,1e20). Those are the bounds used by the compatible
    # preset for species, flow, and counterflow Lambda, respectively.
    cantera_probe=(
        ("h2o2.yaml",3500.0,(200.0,7000.0),(0.2,7.0)),
        ("glarborg_2018.yaml",2000.0,(200.0,4000.0),(0.2,4.0)),
        ("creck_nh3.yaml",3500.0,(200.0,7000.0),(0.2,7.0)),
    )
    for (_,minimum_species_Tmax,cantera_T,arrhenius_T) in cantera_probe
        @test cantera_T==(200.0,2minimum_species_Tmax)
        @test arrhenius_T==(0.2,2minimum_species_Tmax/1000)
        @test cantera_T[1]/1000==arrhenius_T[1]
    end
    @test Arrhenius._cantera_flame_temperature_bounds(f)==(0.2,7.0)
    cantera_Y=(-1e-7,1e5)
    cantera_velocity=(-1e20,1e20)
    cantera_lambda=(-1e20,1e20)
    @test Arrhenius._cantera_flame_species_bounds()==cantera_Y
    @test Arrhenius._cantera_flame_velocity_bounds()==cantera_velocity
    velocity_scale=Arrhenius._flame_state_density(f,f.state,1)
    preset_flow_bounds=Arrhenius._cantera_flame_flow_state_bounds(f,f.state,1)
    @test collect(preset_flow_bounds./velocity_scale) ≈ collect(cantera_velocity) rtol=eps(Float64)
    burner=BurnerFlame(gas;X="H2:2,O2:1,AR:3",mdot=.1,
        grid=collect(range(0.,.5;length=5)),discretization=:finite_difference)
    burner_density=Arrhenius._flame_state_density(burner,burner.state,1)
    burner_flow_bounds=Arrhenius._cantera_flame_flow_state_bounds(burner,burner.state,1)
    @test collect(burner_flow_bounds./burner_density) ≈ collect(cantera_velocity) rtol=eps(Float64)
    @test (cantera_velocity[1]/Arrhenius._counterflow_Vscale,
        cantera_velocity[2]/Arrhenius._counterflow_Vscale)==(-1e18,1e18)
    @test (cantera_lambda[1]/Arrhenius._counterflow_Lscale,
        cantera_lambda[2]/Arrhenius._counterflow_Lscale)==(-1e16,1e16)
    # At s1 == 1, Cantera accepts only when the correction also decreases;
    # convergence then follows from the accepted damping outcome.
    @test Arrhenius._newton_merit_accepts(compatible,1.0,1.1,true)
    @test !Arrhenius._newton_merit_accepts(compatible,1.0,0.9,true)

    # Start each hand trace from a converged premixed/counterflow state, then perturb its
    # center temperature by 2% so the Newton call must accept a damped step and
    # subsequently report convergence under Cantera's weighted-step rule.
    premixed=FreeFlame(gas;X,grid=collect(range(0.,1.;length=17)),discretization=:finite_difference)
    solve!(premixed;refine_grid=false,auto=false,max_time_steps=100,cantera_compatible=true)
    @test premixed.converged
    premixed.state[1,9] *= 1.02
    premixed_trace=Tuple{Bool,Float64,Float64,Bool}[]
    premixed_converged=Arrhenius._flame_newton!(premixed,Arrhenius.FlameWorkspace(premixed);maxiters=20,
        newton_options=compatible,step_trace=premixed_trace)
    @test !isempty(premixed_trace)
    # Cantera accepts the first step because 3.067909053575063 <
    # 23.39533432218312, but it does not converge because s1 >= 1.
    # It accepts and converges the second because 0.3086583470757785 < 1.
    expected_premixed_sequence=[(true,false),(true,true)]
    @test [(accepted,converged) for (accepted,_,_,converged) in premixed_trace]==expected_premixed_sequence
    @test premixed_converged==last(expected_premixed_sequence)[2]

    cf=CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1,AR:3",
        mdot_fuel=.12,mdot_oxidizer=.12,P=.05*one_atm,width=.02,
        grid=collect(range(0.,.02;length=5)))
    solve!(cf;refine_grid=false,auto=false,max_time_steps=100,cantera_compatible=true)
    @test cf.converged
    cf.state[1,3] *= 1.02
    counterflow_trace=Tuple{Bool,Float64,Float64,Bool}[]
    counterflow_converged=Arrhenius._counterflow_newton!(cf,Arrhenius.CounterflowWorkspace(cf);maxiters=20,
        analytic=false,newton_options=compatible,step_trace=counterflow_trace)
    @test !isempty(counterflow_trace)
    # Cantera accepts the first step because 8.315338762561373 <
    # 35.61845052974506, but s1 >= 1 so it does not converge.
    # The next step converges because 0.25362760981519045 < 1.
    expected_counterflow_sequence=[(true,false),(true,true)]
    @test [(accepted,converged) for (accepted,_,_,converged) in counterflow_trace]==expected_counterflow_sequence
    @test counterflow_converged==last(expected_counterflow_sequence)[2]

    wcurrent=Arrhenius._flame_correction_weights(f,f.state,false,current)
    wvelocity=Arrhenius._flame_correction_weights(f,f.state,false,compatible)
    @test wcurrent[end] != wvelocity[end]
    step=ones(size(f.state))
    ncurrent=Arrhenius._flame_policy_correction_norm(step,wcurrent,f,current)
    ncompatible=Arrhenius._flame_policy_correction_norm(step,wvelocity,f,compatible)
    @test ncurrent != ncompatible
    @test_throws ArgumentError solve!(f;newton_options=:invalid)
end
