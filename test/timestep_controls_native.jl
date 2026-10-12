using Arrhenius,Test

@testset "Cantera time-step controls" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    X="H2:1.1,O2:1,AR:5"
    defaults=(;time_step_growth_factor=1.5,time_step_growth_strategy=:fixed_growth,
        time_step_callback=nothing,time_step_regrid=0,jacobian_mode=:auto)
    reference=solve!(FreeFlame(gas;X))
    explicit=solve!(FreeFlame(gas;X);defaults...)
    @test explicit.grid==reference.grid && explicit.state==reference.state
    speed=flame_speed(reference)
    for strategy in Arrhenius._timestep_growth_strategies
        steps=Float64[]
        f=solve!(FreeFlame(gas;X);time_step_growth_strategy=strategy,time_step_growth_factor=2.,
            time_step_callback=dt->push!(steps,dt))
        @test f.converged && isapprox(flame_speed(f),speed;rtol=1e-6)
        @test !isempty(steps) && all(>(0),steps)
    end
    fd=solve!(FreeFlame(gas;X);jacobian_mode=:finite_difference)
    @test fd.converged && isapprox(flame_speed(fd),speed;rtol=1e-6)
    @test solve!(FreeFlame(gas;X);time_step_regrid=3).state==reference.state
    for bad in ((;time_step_growth_factor=.5),(;time_step_growth_factor=Inf),
            (;time_step_growth_strategy=:invalid),(;time_step_regrid=-1),(;jacobian_mode=:analytic))
        @test_throws ArgumentError solve!(FreeFlame(gas;X);bad...)
        @test_throws ArgumentError solve!(CounterflowDiffusionFlame(gas;fuel="H2:1",oxidizer="O2:1",
            mdot_fuel=.24,mdot_oxidizer=.72);bad...)
    end

    kwargs=(fuel="H2:1,AR:1",oxidizer="O2:.2,AR:.8",mdot_fuel=.24,mdot_oxidizer=.72)
    counterflow=solve!(CounterflowDiffusionFlame(gas;kwargs...);slope=.2,curve=.3)
    explicit=solve!(CounterflowDiffusionFlame(gas;kwargs...);slope=.2,curve=.3,defaults...)
    @test explicit.grid==counterflow.grid && explicit.state==counterflow.state
    Tmax=maximum(temperature(counterflow))
    for options in ((;time_step_growth_strategy=:residual_ratio),(;time_step_growth_strategy=:newton_iterations),
            (;jacobian_mode=:finite_difference))
        steps=Float64[]
        c=solve!(CounterflowDiffusionFlame(gas;kwargs...);slope=.2,curve=.3,options...,
            time_step_callback=dt->push!(steps,dt))
        @test c.converged && isapprox(maximum(temperature(c)),Tmax;rtol=1e-6)
        @test !isempty(steps) && all(>(0),steps)
    end
end
