@testset "Ideal-gas isentropes and sound speed" begin
    gas = CreateSolution(joinpath(@__DIR__,"..","mechanism","gri30.yaml"))
    x = mole_fractions(gas,"AR:1")
    mw = dot(gas.MW,x)
    for ratio in (.01,.2,1.,2.,10.)
        state = isentropic_state(gas;T=1000.,P=one_atm,X=x,pressure=ratio*one_atm)
        @test state.T ≈ 1000ratio^(2/5) rtol=2e-12
        @test state.X == x
        @test cal_smass_mean(gas,state.T,state.P,state.X) ≈ cal_smass_mean(gas,1000.,one_atm,x) rtol=2e-13
        @test frozen_sound_speed(gas;T=state.T,P=state.P,X=x) ≈ sqrt((5/3)*R*state.T/mw) rtol=2e-13
    end
    for composition in ("H2:1,N2:0.1","CH4:1,O2:2,N2:7.52")
        state = isentropic_state(gas;T=1200.,P=10one_atm,X=composition,pressure=.2one_atm)
        reverse = isentropic_state(gas;T=state.T,P=state.P,X=state.X,pressure=10one_atm)
        @test reverse.T ≈ 1200. rtol=2e-12
    end
    eq = equilibrium_sound_speeds(gas;T=1000.,X="AR:1")
    @test eq.equilibrium ≈ eq.frozen rtol=2e-8
    @test eq.equilibrium ≈ frozen_sound_speed(gas;T=1000.,X="AR:1") rtol=3e-5
    @test_throws ArgumentError frozen_sound_speed(gas;T=-1.,X=x)
    @test_throws ArgumentError isentropic_state(gas;T=300.,X=x,pressure=-1.)
    @test_throws ArgumentError isentropic_state(gas;T=300.,X=x,pressure=1e-5,temperature_bounds=(200.,6000.))
    @test_throws ArgumentError equilibrium_sound_speeds(gas;T=300.,X=x,pressure_step=0.)
    @test_throws ArgumentError equilibrate(gas;T=300.,X=x,property_rtol=0.)
end
