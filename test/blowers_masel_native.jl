using Arrhenius, Test, ForwardDiff
@testset "Blowers–Masel barrier and temperature response" begin
    rate = BlowersMaselRate(38.7,2.7,26_191_840.,1e9)
    @test activation_energy(rate,0.) ≈ rate.Ea0 rtol=1e-14
    @test activation_energy(rate,-5rate.Ea0) == 0
    @test activation_energy(rate,5rate.Ea0) == 5rate.Ea0
    @test activation_energy(rate,-4rate.Ea0) == 0
    @test activation_energy(rate,nextfloat(4rate.Ea0)) == nextfloat(4rate.Ea0)
    @test ForwardDiff.derivative(h->activation_energy(rate,h),0.) ≈ .5 rtol=1e-14
    for T in (300.,1000.,3400.)
        value = rate_constant(rate,T,0.)
        @test value ≈ rate.A*T^rate.b*exp(-rate.Ea0/(R*T)) rtol=1e-13
        @test ForwardDiff.derivative(t->rate_constant(rate,t,0.),T) ≈
            value*(rate.b/T+rate.Ea0/(R*T^2)) rtol=1e-13
    end
    @test_throws ArgumentError BlowersMaselRate(1.,0.,1.,1.)
    @test_throws ArgumentError BlowersMaselRate(1.,0.,-1.,2.)
    @test_throws ArgumentError rate_constant(rate,0.,1.)
    @test_throws ArgumentError rate_constant(rate,300.,NaN)
end
