@testset "plasma module boundary and compatibility API" begin
    @test Arrhenius.Plasma.PlasmaMechanism === Arrhenius.PlasmaMechanism
    @test Arrhenius.Plasma.PlasmaState === Arrhenius.PlasmaState
    @test Arrhenius.Plasma.solve_eedf === Arrhenius.solve_eedf
    @test Arrhenius.Plasma.set_reduced_electric_field! === Arrhenius.set_reduced_electric_field!
    @test isdefined(Arrhenius, :PlasmaEnergyReactor)
    @test Arrhenius.Plasma.SurfaceMechanism === Arrhenius.SurfaceMechanism
    @test Arrhenius.Plasma.coverage_thermo === Arrhenius.coverage_thermo
end
