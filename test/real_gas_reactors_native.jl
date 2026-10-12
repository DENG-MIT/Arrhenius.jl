using Arrhenius, LinearAlgebra, YAML, Test
for (name,file) in ((:RedlichKwongThermo,"RealGasThermo.jl"),(:RedlichKwongReactor,"RealGasReactors.jl"))
    isdefined(Arrhenius,name) || Base.include(Arrhenius,joinpath(@__DIR__,"..","src",file))
end

@testset "Native RK reactor conservation and callbacks" begin
    path = joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas = CreateSolution(path)
    function fixture(ideal=false)
        yaml = YAML.load_file(path)
        yaml["phases"][1]["thermo"] = "Redlich-Kwong"
        for (i,species) in enumerate(yaml["species"])
            # Numerical test coefficients, not a proposed hydrogen EOS fit.
            species["equation-of-state"] = Dict("model"=>"Redlich-Kwong",
                "a"=>string(ideal ? 0. : 1e5*(1+i/10))*" Pa*m^6*K^0.5/kmol^2",
                "b"=>string(ideal ? 0. : .015+.001i)*" m^3/kmol")
        end
        return RedlichKwongThermo(yaml;molecular_weights=gas.MW)
    end
    model = fixture()
    X = collect(1.:gas.n_species)
    X ./= sum(X)
    reactor = RedlichKwongReactor(gas,model;temperature=1400.,pressure=40one_atm,mole_fractions=X)
    u = reactor_state(reactor)
    rhs = reactor_rhs(reactor)
    du = similar(u)
    rhs(du,u,nothing,0.)
    props = reactor_properties(reactor)
    state = redlich_kwong_state(model;T=u[end],rho=reactor.density,X)
    @test abs(sum(du[1:end-1])) <= 1e-12*norm(du[1:end-1],1)
    @test norm(gas.ele_matrix*(du[1:end-1]./gas.MW),Inf) <= 1e-12*norm(du[1:end-1],1)
    source = dot(state.u_TV./gas.MW,du[1:end-1])
    @test abs(source+state.cv_mass*du[end]) <= 1e-12*abs(source)
    @test props.internal_energy ≈ state.u_mass
    @test props.density == reactor.density
    @test props.pressure ≈ 40one_atm
    @test props.mass_fraction_sum ≈ 1.
    @test props.elemental_inventory ≈ gas.ele_matrix*(reactor.mass_fractions./gas.MW)
    @test reactor_state(reactor) !== reactor.mass_fractions
    J = zeros(length(u),length(u))
    before = copy(u)
    reactor_jacobian!(J,u,rhs)
    @test u == before
    direction = collect(range(-.1,.1;length=length(u)))
    direction[end] = 100.
    plus,minus = similar(u),similar(u)
    step = 1e-5
    rhs(plus,u+step*direction,nothing,0.)
    rhs(minus,u-step*direction,nothing,0.)
    @test J*direction ≈ (plus-minus)/(2step) rtol=2e-6
    problem = reactor_problem(reactor,(0.,.005))
    @test problem.u0 == reactor_state(reactor)
    @test problem.u0 !== u
    problem.tgrad(du,u,nothing,0.)
    @test all(iszero,du)
    @test solve_reactor(reactor,(0.,.005);integrator=p -> p.u0) == reactor_state(reactor)
    inactive = RedlichKwongReactor(gas,model;temperature=1400.,mole_fractions=X,rate_multipliers=zeros(gas.n_reactions))
    reactor_rhs(inactive)(du,reactor_state(inactive),nothing,0.)
    @test all(iszero,du)
    from_mass = RedlichKwongReactor(gas,model;temperature=1400.,pressure=40one_atm,mass_fractions=reactor.mass_fractions)
    @test reactor_state(from_mass) ≈ reactor_state(reactor)
    @test from_mass.density ≈ reactor.density
    for kw in ((temperature=-1.,mole_fractions=X),(temperature=1400.,),
        (temperature=1400.,mole_fractions=X,mass_fractions=X),
        (temperature=1400.,mole_fractions=X,rate_multipliers=[1.]),
        (temperature=1400.,mole_fractions=X,rate_multipliers=fill(-1.,gas.n_reactions)))
        @test_throws ArgumentError RedlichKwongReactor(gas,model;kw...)
    end
    @test_throws ArgumentError reactor_problem(reactor,(1.,0.))
    @test_throws DimensionMismatch rhs(zeros(2),u,nothing,0.)
    # Zero attraction/covolume recovers the established ideal-gas reactor RHS.
    ideal_rk = RedlichKwongReactor(gas,fixture(true);temperature=1400.,pressure=40one_atm,mole_fractions=X)
    ideal = IdealGasReactor(gas;temperature=1400.,pressure=40one_atm,mole_fractions=X,constraint=:constant_volume)
    reactor_rhs(ideal_rk)(plus,reactor_state(ideal_rk),nothing,0.)
    reactor_rhs(ideal)(minus,reactor_state(ideal),nothing,0.)
    @test plus ≈ minus rtol=2e-12
    @test_throws ArgumentError Arrhenius._check_real_gas_rate_models((blowers_masel=(reaction_indices=[1],),))
    @test Arrhenius._check_real_gas_rate_models((blowers_masel=(reaction_indices=Int[],),)) === nothing
    if hasproperty(gas.reaction,:blowers_masel)
        # Exercise actual public constructors on the combined core snapshot.
        old = gas.reaction
        fields = ntuple(i -> getfield(old,i),fieldcount(typeof(old))-1)
        bm = Arrhenius.BlowersMaselData([1],reshape([1.,0.,1e7,1e9],1,4))
        reaction = Arrhenius.Reaction(fields...,bm)
        badgas = Arrhenius.Solution(gas.n_species,gas.n_reactions,gas.MW,gas.species_names,
            gas.elements,gas.ele_matrix,gas.thermo,gas.trans,reaction)
        @test_throws ArgumentError RedlichKwongReactor(badgas,model;temperature=1400.,mole_fractions=X)
        @test_throws ArgumentError RealGasKineticsWorkspace(badgas,model)
    end
end
