using Arrhenius, LinearAlgebra, Test

@testset "native composition validation" begin
    gas = CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    X = mole_fractions(gas,"H2:2, O2:1, AR:3")
    @test X == mole_fractions(gas,Dict(:H2=>2,:O2=>1,:AR=>3))
    @test X == mole_fractions(gas,"H2:1,H2:1,O2:1,AR:3")
    @test X ≈ mole_fractions(gas,mass_fractions(gas,X);basis=:mass)
    @test X ≈ mole_fractions(gas,"H2:2e307,O2:1e307,AR:3e307")
    before = copy(X)
    mole_fractions(gas,X)
    @test X == before
    @test sum(mole_fractions(gas,"AR")) == 1
    for invalid in ("", "H2:-1", "H2:2,H2:-1", "H2:NaN", "H2:Inf", "H2:bad", "missing:1", "H2:0")
        @test_throws ArgumentError mole_fractions(gas,invalid)
    end
    @test_throws DimensionMismatch mole_fractions(gas,[1.,2.])
    @test_throws ArgumentError mole_fractions(gas,zeros(gas.n_species))
    @test_throws ArgumentError mole_fractions(gas,X;basis=:volume)
    @test_throws ArgumentError set_equivalence_ratio(gas,-1;fuel="H2",oxidizer="O2")
    @test_throws ArgumentError set_equivalence_ratio(gas,1;fuel="O2",oxidizer="H2")
    @test_throws ArgumentError set_equivalence_ratio(gas,1;fuel="AR",oxidizer="O2")
    @test_throws ArgumentError set_equivalence_ratio(gas,1;fuel="H2",oxidizer="O2",fraction=(diluent=.5,))
    @test_throws ArgumentError set_equivalence_ratio(gas,1;fuel="H2",oxidizer="O2",diluent="AR")
    @test_throws ArgumentError set_equivalence_ratio(gas,1;fuel="H2",oxidizer="O2",diluent="AR",fraction=(fuel=.99,))
    @test_throws ArgumentError set_mixture_fraction(gas,1.1;fuel="H2",oxidizer="O2")
    @test_throws ArgumentError equivalence_ratio(gas,X;fuel="H2")
    @test_throws ArgumentError mixture_fraction(gas,X;fuel="H2",oxidizer="O2",element="C")
    @test_throws ArgumentError mixture_fraction(gas,X;fuel="AR",oxidizer="AR")
    @test equivalence_ratio(gas,mole_fractions(gas,"O2")) == 0
    @test equivalence_ratio(gas,mole_fractions(gas,"H2")) == Inf
    @test isnan(equivalence_ratio(gas,mole_fractions(gas,"AR")))
    @test set_equivalence_ratio(gas,0;fuel="H2",oxidizer="O2") == mole_fractions(gas,"O2")
    @test set_mixture_fraction(gas,0;fuel="H2",oxidizer="O2") == mole_fractions(gas,"O2")
    @test set_mixture_fraction(gas,1;fuel="H2",oxidizer="O2") == mole_fractions(gas,"H2")
    for basis in (:mole,:mass), phi in (.1,1.,10.)
        X = set_equivalence_ratio(gas,phi;fuel="H2:1,H2O:.2",oxidizer="O2:1,AR:3",basis)
        @test equivalence_ratio(gas,X;fuel="H2:1,H2O:.2",oxidizer="O2:1,AR:3",basis) ≈ phi rtol=1e-12
    end
end

@testset "native equilibrium invariants and degenerate element systems" begin
    gas = CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    X0 = mole_fractions(gas,"H2:2,O2:1,AR:3")
    initial = copy(X0)
    baseline = equilibrate(gas;T=2500.,P=one_atm,X=X0)
    @test X0 == initial
    @test sum(baseline.X) ≈ 1 atol=1e-14
    @test sum(baseline.Y) ≈ 1 atol=1e-14
    @test minimum(baseline.X) >= 0
    @test baseline == equilibrate(gas;T=2500.,P=one_atm,X=X0,mode="TP")
    # Chemical potentials must lie in the span of elemental potentials.
    mu = cal_g_RT(gas,baseline.T,baseline.P,baseline.X)
    potentials = transpose(gas.ele_matrix) \ mu
    @test norm(transpose(gas.ele_matrix)*potentials-mu,Inf) < 1e-8
    second = equilibrate(gas;T=baseline.T,P=baseline.P,X=baseline.X)
    @test second.X ≈ baseline.X rtol=1e-8 atol=1e-13

    # A redundant element row must not create an extra conservation constraint.
    redundant = Arrhenius.Solution(gas.n_species,gas.n_reactions,gas.MW,
        gas.species_names,[gas.elements;"duplicate"],
        vcat(gas.ele_matrix,2 .* gas.ele_matrix[1:1,:]),gas.thermo,gas.trans,gas.reaction)
    result = equilibrate(redundant;T=2500.,P=one_atm,X=X0)
    @test result.X ≈ baseline.X rtol=1e-8 atol=1e-13
    @test size(Arrhenius._equilibrium_system(redundant,X0).A,1) == size(gas.ele_matrix,1)

    for species in ("AR","H2","O2")
        X = mole_fractions(gas,species)
        result = equilibrate(gas;T=1500.,X)
        absent = findall(==(0),gas.ele_matrix*X)
        for k in eachindex(result.X)
            any(gas.ele_matrix[e,k] != 0 for e in absent) && @test result.X[k] == 0
        end
    end
    for mode in (:TP,:TV,:HP,:UV,:SP,:SV)
        result = equilibrate(gas;T=1200.,P=3one_atm,X="AR",mode)
        @test result.T ≈ 1200. atol=1e-7
        @test result.P ≈ 3one_atm rtol=1e-12
        @test result.X == mole_fractions(gas,"AR")
    end
    # Single molecular species with multiple dependent element rows.
    water = findfirst(==("H2O"),gas.species_names)
    inds = [water]
    thermo = Arrhenius.IdealGasThermo(gas.thermo.nasa_low[inds,:],
        gas.thermo.nasa_high[inds,:],gas.thermo.Trange[inds,:],true)
    single = Arrhenius.Solution(1,gas.n_reactions,gas.MW[inds],gas.species_names[inds],
        gas.elements,gas.ele_matrix[:,inds],thermo,gas.trans,gas.reaction)
    for mode in (:TP,:TV,:HP,:UV,:SP,:SV)
        result = equilibrate(single;T=800.,P=one_atm,X=[1.],mode)
        @test result.T ≈ 800. atol=1e-7
        @test result.X == [1.]
    end
    # A boundary of feasible composition can force a species to vanish even
    # though no element is absent: pure H2O in an H2O/H2O2-only phase.
    inds = [water,findfirst(==("H2O2"),gas.species_names)]
    thermo = Arrhenius.IdealGasThermo(gas.thermo.nasa_low[inds,:],
        gas.thermo.nasa_high[inds,:],gas.thermo.Trange[inds,:],true)
    boundary = Arrhenius.Solution(2,gas.n_reactions,gas.MW[inds],gas.species_names[inds],
        gas.elements,gas.ele_matrix[:,inds],thermo,gas.trans,gas.reaction)
    @test equilibrate(boundary;T=2500.,X=[1.,0.]).X ≈ [1.,0.] atol=1e-9

    @test_throws ArgumentError equilibrate(gas;T=0.,X=X0)
    @test_throws ArgumentError equilibrate(gas;T=Inf,X=X0)
    @test_throws ArgumentError equilibrate(gas;T=300.,P=-1.,X=X0)
    @test_throws ArgumentError equilibrate(gas;T=300.,X=X0,mode=:invalid)
    @test_throws ArgumentError equilibrate(gas;T=300.,X=X0,mode=:HP,temperature_bounds=(400.,300.))
    @test_throws ArgumentError equilibrate(gas;T=300.,X=X0,mode=:HP,temperature_bounds=(290.,400.))
end
