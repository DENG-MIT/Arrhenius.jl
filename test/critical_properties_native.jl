@testset "pure-fluid critical properties" begin
    include(joinpath(@__DIR__, "..", "example", "thermodynamics", "critical_properties.jl"))
    output = critical_properties_calculation()
    @test size(output) == (5, 8)
    @test all(isfinite, output)
    @test all(>(0), output)
    @test all(0 .< output[5, :] .< 1)
    @test critical_properties("HFC-134a") == critical_properties("hfc-134a")
    @test critical_properties("WATER").T == Arrhenius.WATER_TC
    @test critical_properties("water").P == Arrhenius.WATER_PC
    @test critical_properties("water").density == Arrhenius.WATER_RHOC
    @test critical_properties("water").molecular_weight == Arrhenius.WATER_MW
    @test_throws ArgumentError critical_properties("argon")
    @test_throws ArgumentError critical_properties("")
    # Preserve the critical point, not HFC-134a's EOS reducing parameters.
    refrigerant = critical_properties("HFC-134a")
    @test refrigerant.T == 374.21
    @test refrigerant.P == 4059280.0
    @test refrigerant.density == 511.95
    # Inputs determine column order; a returned table cannot mutate model data.
    @test critical_properties_calculation(reverse(CRITICAL_FLUID_NAMES)) == reverse(output; dims=2)
    output[1, 1] = 0
    @test critical_properties("water").T == Arrhenius.WATER_TC
end
