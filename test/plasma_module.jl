@testset "plasma module boundary and compatibility API" begin
    legacy_path = joinpath(@__DIR__, "data", "plasma_legacy_exports.txt")
    legacy_exports = Set(Symbol.(strip.(readlines(legacy_path))))
    for name in legacy_exports
        @test isdefined(Arrhenius, name)
    end

    # Every name exported by Plasma is the same binding from the legacy module.
    for name in names(Arrhenius.Plasma)
        @test isdefined(Arrhenius, name)
        @test getfield(Arrhenius.Plasma, name) === getfield(Arrhenius, name)
    end

    @test Arrhenius.Plasma.PlasmaMechanism === Arrhenius.PlasmaMechanism
    @test Arrhenius.Plasma.PlasmaState === Arrhenius.PlasmaState
    @test Arrhenius.Plasma.solve_eedf === Arrhenius.solve_eedf
    @test Arrhenius.Plasma.set_reduced_electric_field! === Arrhenius.set_reduced_electric_field!
    @test Arrhenius.Plasma.SurfaceMechanism === Arrhenius.SurfaceMechanism
    @test Arrhenius.Plasma.coverage_thermo === Arrhenius.coverage_thermo
end
