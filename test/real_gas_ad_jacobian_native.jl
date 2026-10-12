using Arrhenius, LinearAlgebra, YAML, Test
include(joinpath(@__DIR__,"..","example","reactors","real_gas_ad_jacobian.jl"))

@testset "Optional shock-tube automatic Jacobian" begin
    path=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(path)
    yaml=YAML.load_file(path)
    yaml["phases"][1]["thermo"]="Redlich-Kwong"
    for (i,species) in enumerate(yaml["species"])
        # Numerical derivative fixture; these coefficients are not an EOS fit.
        species["equation-of-state"]=Dict("model"=>"Redlich-Kwong",
            "a"=>string(1e5*(1+i/10))*" Pa*m^6*K^0.5/kmol^2",
            "b"=>string(.015+.001i)*" m^3/kmol")
    end
    model=RedlichKwongThermo(yaml;molecular_weights=gas.MW)
    for real in (false,true), absent in (false,true)
        composition=absent ? Dict("H2"=>2.,"O2"=>1.,"AR"=>4.) : collect(1.:gas.n_species)
        reactor=real ? RedlichKwongReactor(gas,model;temperature=1200.,pressure=40one_atm,mole_fractions=composition) :
            IdealGasReactor(gas;temperature=1200.,pressure=40one_atm,mole_fractions=composition,constraint=:constant_volume)
        u=reactor_state(reactor)
        original=copy(u)
        jac=shocktube_ad_jacobian(reactor)
        J=zeros(length(u),length(u))
        jac(J,u,nothing,0.)
        @test u==original
        @test all(isfinite,J)
        direction=vcat(0.01 .* abs.(sin.(collect(1.:gas.n_species))),200.)
        rhs=reactor_rhs(reactor)
        base,plus,minus=zero(u),zero(u),zero(u)
        step=1e-5
        rhs(base,u,nothing,0.)
        rhs(plus,u.+step.*direction,nothing,0.)
        if absent
            rhs(minus,u.+2step.*direction,nothing,0.)
            numerical=(-3base.+4plus.-minus)./(2step)
        else
            rhs(minus,u.-step.*direction,nothing,0.)
            numerical=(plus.-minus)./(2step)
        end
        @test maximum(abs.(J*direction-numerical)./max.(abs.(numerical),1))<2e-5
        species=view(J,1:gas.n_species,:)
        @test norm(sum(species,dims=1),Inf)<1e-10*max(norm(species,Inf),1)
        elemental=gas.ele_matrix*(species./gas.MW)
        @test norm(elemental,Inf)<1e-10*max(norm(species,Inf),1)
        # A temperature change must invalidate the cached Arrhenius/Kc factors.
        expected=copy(J)
        shifted=copy(u); shifted[end]+=75.
        jac(J,shifted,nothing,0.)
        @test !isapprox(J,expected;rtol=1e-3)
        jac(J,u,nothing,0.)
        @test J≈expected rtol=2e-13 atol=1e-12
        @test_throws DimensionMismatch jac(zeros(2,2),u,nothing,0.)
        @test_throws ArgumentError shocktube_ad_jacobian(reactor;chunk=0)
    end
    isothermal=IdealGasReactor(gas;temperature=1200.,mole_fractions=Dict("H2"=>1.),
        constraint=:constant_volume,energy=:isothermal)
    @test_throws ArgumentError shocktube_ad_jacobian(isothermal)
    old=gas.reaction
    fields=ntuple(i -> getfield(old,i),fieldcount(typeof(old))-1)
    bm=Arrhenius.BlowersMaselData([1],reshape([1.,0.,1e7,1e9],1,4))
    reaction=Arrhenius.Reaction(fields...,bm)
    badgas=Arrhenius.Solution(gas.n_species,gas.n_reactions,gas.MW,gas.species_names,
        gas.elements,gas.ele_matrix,gas.thermo,gas.trans,reaction)
    bm_reactor=IdealGasReactor(badgas;temperature=1200.,mole_fractions=Dict("H2"=>1.),constraint=:constant_volume)
    @test_throws ArgumentError shocktube_ad_jacobian(bm_reactor)
end
