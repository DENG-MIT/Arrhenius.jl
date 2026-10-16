using Arrhenius, NPZ, Test
const CF=Arrhenius
gas=CreateSolution(joinpath(ARGS[1],"mechanism","gri30.yaml"))
@testset "Cantera H2O/CO2 radiation source, isolated property evaluation" begin
    for path in ARGS[2:end]
        ref=npzread(path)
        # Property evaluation only. No nonlinear solve follows reference loading.
        probe=CF.CounterflowDiffusionFlame(gas;fuel="C2H6:1",oxidizer="O2:.21,N2:.78,AR:.01",
            mdot_fuel=.24,mdot_oxidizer=.72,grid=ref["grid"],radiation=true,
            boundary_emissivities=Tuple(ref["boundary_emissivities"]))
        probe.state[1,:] .= ref["T"]./1000
        probe.state[2:gas.n_species+1,:] .= ref["Y"]
        q=CF.radiative_heat_loss(probe)
        error=maximum(abs.(q.-ref["radiative_heat_loss"]))/maximum(abs.(ref["radiative_heat_loss"]))
        @test error<1e-12
        @test q[end]==0
        @test maximum(q)>0
        @test_throws ArgumentError CF.set_radiation!(probe,true;boundary_emissivities=(NaN,0))
        CF.set_radiation!(probe,false)
        @test all(iszero,CF.radiative_heat_loss(probe))
        println(basename(path)," radiation source normalized error: ",error)
    end
end
