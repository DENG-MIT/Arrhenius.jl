module IdealGasMixingTests

using Arrhenius, LinearAlgebra, Test

const mix_h2o2 = CreateSolution(joinpath(@__DIR__, "..", "mechanism", "h2o2.yaml"))

# Independent totals via species-vector contractions, not the mixture means
# the implementation accumulates.
function stream_enthalpy(gas, T, P, X, moles)
    x = mole_fractions(gas, X)
    return moles * dot(x, cal_h(gas, T, P, x))
end

@testset "mass, species and enthalpy conservation" begin
    gas = mix_h2o2
    streams = ((T=900.0, P=one_atm, X="H2:2, O2:1", moles=0.5),
               (T=300.0, P=one_atm, X="AR:1, H2O:0.1", moles=1.5))
    result = mix_constant_pressure(gas, streams)
    xa = mole_fractions(gas, streams[1].X)
    xb = mole_fractions(gas, streams[2].X)
    @test result.P === one_atm
    @test result.moles == 2.0
    @test result.mass ≈ 0.5 * dot(gas.MW, xa) + 1.5 * dot(gas.MW, xb)
    @test result.X ≈ (0.5 .* xa .+ 1.5 .* xb) ./ 2.0
    @test result.Y ≈ result.X .* gas.MW ./ dot(result.X, gas.MW)
    @test sum(result.X) ≈ 1
    @test sum(result.Y) ≈ 1
    @test 300.0 < result.T < 900.0
    @test result.enthalpy ≈ stream_enthalpy(gas, 900.0, one_atm, xa, 0.5) +
                            stream_enthalpy(gas, 300.0, one_atm, xb, 1.5) rtol=1e-12
    @test result.enthalpy ≈ result.moles * dot(result.X, cal_h(gas, result.T, result.P, result.X)) rtol=1e-9
end

@testset "equal inlet temperatures are returned exactly" begin
    gas = mix_h2o2
    result = mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="H2:1", moles=0.3),
                                         (T=300.0, P=one_atm, X="O2:1, AR:1", moles=0.7)))
    @test result.T === 300.0
    xa = mole_fractions(gas, "H2:1")
    xb = mole_fractions(gas, "O2:1, AR:1")
    @test result.X ≈ 0.3 .* xa .+ 0.7 .* xb
    # A stream with zero amount does not affect the mixed state.
    with_inert = mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="H2:1", moles=0.3),
                                             (T=300.0, P=one_atm, X="O2:1, AR:1", moles=0.7),
                                             (T=1500.0, P=one_atm, X="H2O:1", moles=0.0)))
    @test with_inert.T === 300.0
    @test with_inert.X == result.X
    @test with_inert.moles == result.moles
end

@testset "single stream returns its own state" begin
    gas = mix_h2o2
    result = mix_constant_pressure(gas, [(T=500.0, P=2one_atm, X="H2:1, O2:0.5", moles=2.0)])
    @test result.T === 500.0
    @test result.P === 2one_atm
    @test result.X == mole_fractions(gas, "H2:1, O2:0.5")
    @test result.moles == 2.0
    @test result.enthalpy ≈ stream_enthalpy(gas, 500.0, 2one_atm, "H2:1, O2:0.5", 2.0) rtol=1e-13
end

@testset "stream order independence" begin
    gas = mix_h2o2
    streams = ((T=900.0, P=one_atm, X="H2:2, O2:1", moles=0.5),
               (T=300.0, P=one_atm, X="AR:1, H2O:0.1", moles=1.5),
               (T=1200.0, P=one_atm, X="O2:1, H2O2:0.2", moles=0.8))
    forward = mix_constant_pressure(gas, streams)
    backward = mix_constant_pressure(gas, reverse(collect(streams)))
    @test backward.T ≈ forward.T rtol=1e-12
    @test backward.X ≈ forward.X
    @test backward.moles == forward.moles
    @test backward.mass ≈ forward.mass
    @test backward.enthalpy ≈ forward.enthalpy rtol=1e-12
end

@testset "scaling all amounts scales only the totals" begin
    gas = mix_h2o2
    streams = ((T=900.0, P=one_atm, X="H2:2, O2:1", moles=0.5),
               (T=300.0, P=one_atm, X="AR:1, H2O:0.1", moles=1.5))
    base = mix_constant_pressure(gas, streams)
    scaled = mix_constant_pressure(gas, ((T=s.T, P=s.P, X=s.X, moles=5 * s.moles) for s in streams))
    @test scaled.T ≈ base.T rtol=1e-12
    @test scaled.X ≈ base.X
    @test scaled.moles ≈ 5 * base.moles
    @test scaled.mass ≈ 5 * base.mass
    @test scaled.enthalpy ≈ 5 * base.enthalpy rtol=1e-12
end

@testset "inputs are not modified" begin
    gas = mix_h2o2
    xa = mole_fractions(gas, "H2:2, O2:1")
    xb = mole_fractions(gas, "AR:1, H2O:0.5")
    xa_before, xb_before = copy(xa), copy(xb)
    mix_constant_pressure(gas, ((T=800.0, P=one_atm, X=xa, moles=0.4),
                                (T=350.0, P=one_atm, X=xb, moles=1.1)))
    @test xa == xa_before
    @test xb == xb_before
end

@testset "invalid streams are rejected" begin
    gas = mix_h2o2
    valid = (T=300.0, P=one_atm, X="H2:1", moles=1.0)
    @test_throws ArgumentError mix_constant_pressure(gas, ())
    @test_throws ArgumentError mix_constant_pressure(gas, NamedTuple[])
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=-1.0, P=one_atm, X="H2:1", moles=1.0),))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=NaN, P=one_atm, X="H2:1", moles=1.0),))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=0.0, X="H2:1", moles=1.0),))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=Inf, X="H2:1", moles=1.0),))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="H2:1", moles=-1.0),))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="H2:1", moles=NaN),))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="H2:1", moles=1e308),
                                                          (T=300.0, P=one_atm, X="H2:1", moles=1e308)))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="H2:1", moles=0.0),))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="H2:1"),))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="H2:1", moles="1"),))
    @test_throws ArgumentError mix_constant_pressure(gas, ([300.0, one_atm],))
    @test_throws ArgumentError mix_constant_pressure(gas, (valid, (T=300.0, P=1.000000001one_atm, X="O2:1", moles=1.0)))
    @test_throws ArgumentError mix_constant_pressure(gas, ((T=300.0, P=one_atm, X="bogus:1", moles=1.0),))
    @test_throws DimensionMismatch mix_constant_pressure(gas, ((T=300.0, P=one_atm, X=[1.0, 2.0], moles=1.0),))
    # A pressure matched within 1e-10 relative is accepted.
    matched = mix_constant_pressure(gas, (valid, (T=300.0, P=one_atm * (1 + 1e-11), X="O2:1", moles=1.0)))
    @test matched.P === one_atm
end

gri30_path = joinpath(@__DIR__, "..", "mechanism", "gri30.yaml")
if isfile(gri30_path)
    @testset "gri30 stoichiometric mixing and TP equilibration" begin
        gas = CreateSolution(gri30_path)
        air = mole_fractions(gas, "O2:0.21, N2:0.78, AR:0.01")
        nO2 = air[findfirst(==("O2"), gas.species_names)]
        mixed = mix_constant_pressure(gas, ((T=300.0, P=one_atm, X=air, moles=1.0),
                                            (T=300.0, P=one_atm, X="CH4:1", moles=0.5 * nO2)))
        @test mixed.T === 300.0
        @test mixed.P === one_atm
        @test mixed.moles ≈ 1 + 0.5 * nO2
        # Stoichiometric proportions: O2 exactly matches the CH4 demand.
        xch4 = mixed.X[findfirst(==("CH4"), gas.species_names)]
        xo2 = mixed.X[findfirst(==("O2"), gas.species_names)]
        @test xo2 ≈ 2 * xch4 rtol=1e-12
        # The unequal-temperature solve on the full 53-species mechanism.
        hot = mix_constant_pressure(gas, ((T=1400.0, P=one_atm, X="N2:3.76, O2:1", moles=2.0),
                                          (T=400.0, P=one_atm, X="CH4:1", moles=0.3)))
        @test 400.0 < hot.T < 1400.0
        @test hot.enthalpy ≈ stream_enthalpy(gas, 1400.0, one_atm, "N2:3.76, O2:1", 2.0) +
                             stream_enthalpy(gas, 400.0, one_atm, "CH4:1", 0.3) rtol=1e-12
        @test hot.enthalpy ≈ hot.moles * dot(hot.X, cal_h(gas, hot.T, hot.P, hot.X)) rtol=1e-9
        # TP equilibration of the frozen mixture burns the methane away.
        eq = equilibrate(gas; T=mixed.T, P=mixed.P, X=mixed.X, mode=:TP)
        @test eq.T == 300.0
        @test eq.P == one_atm
        @test dot(gas.MW, eq.X) ≈ dot(gas.MW, mixed.X)
        @test eq.X[findfirst(==("CH4"), gas.species_names)] < 1e-6
        # Independent Cantera TP reference is approximately 7.4e-20 O2.
        # Intermediate equilibrium residuals must not leave excess oxidizer.
        @test eq.X[findfirst(==("O2"), gas.species_names)] < 1e-12
    end
end

end
