module TimestepExampleFixture
include("cantera_parity/fstc_rows.jl")
end
using Test, Arrhenius
@testset "Cantera timestep example effective inlet settings" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    for pressure in (7e6,15e6)
        f=TimestepExampleFixture.build_flame(mechanism,pressure)
        # Cantera 4.0.0a2 _onedim.py Domain1D.phase returns self.gas for
        # both inlets. The example reads equal densities and therefore sets
        # mdot_o=.3/rho*rho*3=.9, rather than separate pure-inlet densities.
        @test f.fuel_mass_flux≈.3
        @test f.oxidizer_mass_flux≈.9
        @test f.fuel_temperature==800
        @test f.oxidizer_temperature==711
        @test f.pressure==pressure
        @test f.grid==collect(range(0.,.03;length=20))
        @test f.fuel_Y[findfirst(==("H2"),f.gas.species_names)]≈1
        @test f.oxidizer_Y[findfirst(==("O2"),f.gas.species_names)]≈1
    end
end
@testset "Cantera residual-ratio growth on a small state" begin
    # Cantera SteadyStateSystem.cpp:253-263; oracle norms 32 -> 1 give
    # uncapped factor 2, equal/increasing/zero-after norms give factor 1.
    before=[32.,0.]; after=[1.,0.]; residual=zeros(2)
    residual_at!(r,x,previous,dt)=copyto!(r,x)
    @test Arrhenius._timestep_growth(residual_at!,residual,:residual_ratio,1.5,before,after,1e-5,1)==1.5
    @test Arrhenius._timestep_growth(residual_at!,residual,:residual_ratio,3.,before,after,1e-5,1)≈2.
    @test Arrhenius._timestep_growth(residual_at!,residual,:residual_ratio,3.,after,before,1e-5,1)==1.
    @test Arrhenius._timestep_growth(residual_at!,residual,:residual_ratio,3.,before,before,1e-5,1)==1.
    @test Arrhenius._timestep_growth(residual_at!,residual,:residual_ratio,3.,before,zeros(2),1e-5,1)==1.
end
