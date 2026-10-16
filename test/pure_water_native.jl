using Arrhenius, Test, ForwardDiff
isdefined(Arrhenius,:PureWater) || Base.include(Arrhenius,joinpath(@__DIR__,"..","src","PureWater.jl"))

@testset "native water thermodynamic invariants" begin
    water = PureWater()
    for T in (WATER_TMIN,300.,373.15,500.,600.,640.,647.,WATER_TC)
        sat = water_saturation(water,T)
        @test sat.P > 0
        @test sat.rhof >= sat.rhog > 0
        @test water_saturation_temperature(water,sat.P) ≈ T atol=2e-6
        liquid = water_state(water;T,Q=0.)
        vapor = water_state(water;T,Q=1.)
        @test liquid.g ≈ vapor.g atol=.002 rtol=0
        for Q in (.1,.5,.9)
            state = water_state(water;T,Q)
            @test state.v ≈ (1-Q)*liquid.v+Q*vapor.v rtol=1e-14
            if T < WATER_TC
                @test state.Q ≈ Q atol=1e-13
                @test state.h ≈ (1-Q)*liquid.h+Q*vapor.h rtol=1e-13
                @test state.u ≈ (1-Q)*liquid.u+Q*vapor.u rtol=1e-13
                @test state.s ≈ (1-Q)*liquid.s+Q*vapor.s rtol=1e-13
                @test state.cp == Inf
                @test isnan(state.cv)
                @test state.h-state.u ≈ state.P*state.v rtol=2e-7 atol=1e-4
                inverse = water_state(water;P=state.P,h=state.h)
                @test inverse.T ≈ state.T atol=2e-6
                @test inverse.Q ≈ state.Q atol=2e-7
            end
        end
    end
    for (T,P) in ((300.,one_atm),(450.,5e6),(500.,one_atm),(700.,25e6),(1200.,10e6))
        state = water_state(water;T,P)
        @test state.rho > 0
        @test state.P ≈ P rtol=2e-7
        @test state.h-state.u ≈ state.P*state.v rtol=1e-10 atol=1e-7
        @test state.g ≈ state.h-T*state.s
        @test state.a ≈ state.u-T*state.s
        @test state.cp >= state.cv > 0
        cv = ForwardDiff.derivative(t -> Arrhenius._water_raw(t,state.rho).u,T)
        @test state.cv ≈ cv rtol=5e-6
        pressure_T = ForwardDiff.derivative(t -> Arrhenius._water_pressure(t,state.rho),T)
        pressure_rho = ForwardDiff.derivative(rho -> Arrhenius._water_pressure(T,rho),state.rho)
        @test pressure_rho > 0
        @test state.cp-state.cv ≈ T*pressure_T^2/(state.rho^2*pressure_rho) rtol=3e-4
        for key in (:h,:s,:u)
            inverse = water_state(water;P,NamedTuple{(key,)}((getproperty(state,key),))...)
            @test inverse.T ≈ state.T atol=2e-5
            @test inverse.rho ≈ state.rho rtol=1e-6
        end
        inverse = water_state(water;T,v=state.v)
        @test inverse.h ≈ state.h rtol=1e-13
    end
    # The Clapeyron relation provides an independent saturation-curve check.
    for T in (300.,450.,600.)
        liquid,vapor = water_state(water;T,Q=0.),water_state(water;T,Q=1.)
        delta = .01
        slope = (water_saturation(water,T+delta).P-water_saturation(water,T-delta).P)/(2delta)
        @test slope ≈ (vapor.h-liquid.h)/(T*(vapor.v-liquid.v)) rtol=5e-5
    end
    cycle = water_rankine(water)
    @test 0 < cycle.efficiency < 1-cycle.states[1].T/cycle.states[3].T
    @test cycle.pump_work > 0
    @test cycle.turbine_work > cycle.pump_work
    @test cycle.heat_added-(cycle.states[4].h-cycle.states[1].h) ≈ cycle.turbine_work-cycle.pump_work atol=1e-3
    @test cycle.ideal_states[1].s ≈ cycle.states[1].s atol=1e-5
    @test cycle.ideal_states[2].s ≈ cycle.states[3].s atol=1e-5
end

@testset "native water input validation" begin
    water = PureWater()
    @test_throws ArgumentError water_state(water;T=270.,P=one_atm)
    @test_throws ArgumentError water_state(water;T=1700.,P=one_atm)
    @test_throws ArgumentError water_state(water;T=NaN,P=one_atm)
    @test_throws ArgumentError water_state(water;T=300.,P=-1.)
    @test_throws ArgumentError water_state(water;T=300.,Q=-.1)
    @test_throws ArgumentError water_state(water;T=300.,Q=1.1)
    @test_throws ArgumentError water_state(water;T=700.,Q=.5)
    @test_throws ArgumentError water_state(water;T=300.,v=-1.)
    @test_throws ArgumentError water_state(water;T=300.)
    @test_throws ArgumentError water_state(water;T=300.,P=one_atm,Q=0.)
    @test_throws ArgumentError water_state(water;s=0.,h=0.)
    @test_throws ArgumentError water_state(water;P=one_atm,h=NaN)
    @test_throws ArgumentError water_state(water;P=one_atm,h=1e12)
    @test_throws ArgumentError water_state(water;T=300.,P=water_saturation(water,300.).P)
    @test_throws ArgumentError water_saturation_temperature(water,0.)
    @test_throws ArgumentError water_saturation_temperature(water,2WATER_PC)
    @test_throws ArgumentError water_rankine(water;pump_efficiency=0.)
    @test_throws ArgumentError water_rankine(water;turbine_efficiency=2.)
    @test_throws ArgumentError water_rankine(water;boiler_pressure=100.)
end
