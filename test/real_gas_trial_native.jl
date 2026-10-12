module RealGasTrialContracts

using Arrhenius, ForwardDiff, LinearAlgebra, NPZ, SparseArrays, Test, YAML
include(joinpath(@__DIR__,"..","example","reactors","real_gas_ad_jacobian.jl"))

@testset "Signed mass-action column polynomials" begin
    # The columns represent 2x, 2x^2 and (2x)y. At negative trials, two
    # negative factors suppress the product; a single negative factor remains.
    A=sparse([1,1,1,2],[1,2,3,3],[1.,2.,1.,1.],2,3)
    product=(i,c)->_integer_trial_product(2.,c,rowvals(A),nonzeros(A),nzrange(A,i))
    for (c,values,gradients) in (
            ([2.,3.],[4.,8.,12.],([2.,0.],[8.,0.],[6.,4.])),
            ([0.,3.],[0.,0.,0.],([2.,0.],[0.,0.],[6.,0.])),
            ([0.,0.],[0.,0.,0.],([2.,0.],[0.,0.],[0.,0.])),
            ([-2.,3.],[-4.,0.,-12.],([2.,0.],[0.,0.],[6.,-4.])),
            ([2.,-3.],[4.,8.,-12.],([2.,0.],[8.,0.],[-6.,4.])),
            ([-2.,-3.],[-4.,0.,0.],([2.,0.],[0.,0.],[0.,0.])))
        for i in 1:3
            @test product(i,c)==values[i]
            @test ForwardDiff.gradient(z->product(i,z),c)==gradients[i]
        end
    end
    # Preserve the stored-order grouping; k*(x*y) rounds differently here.
    @test _integer_trial_product(.1,[.2,.3],rowvals(A),nonzeros(A),nzrange(A,3))===.006000000000000001
    for T in (Float32,Float64)
        orders=T.(A)
        for (i,expected) in enumerate(T[4,8,12])
            value=_integer_trial_product(2f0,Float32[2,3],rowvals(orders),nonzeros(orders),nzrange(orders,i))
            @test value===expected
        end
    end

    # Mutations must select the live arity/order without prepared metadata.
    live=sparse([1],[1],[1.],2,1)
    evaluate=c->_integer_trial_product(2.,c,rowvals(live),nonzeros(live),nzrange(live,1))
    @test evaluate([2.,3.])==4.
    live[2,1]=1. # (1) -> (1,1)
    @test evaluate([2.,3.])==12.
    @test ForwardDiff.gradient(evaluate,[2.,3.])==[6.,4.]
    live[1,1]=2. # (2,1) exercises the unchanged generic fallback: 2x^2 y.
    @test evaluate([2.,3.])==24.
    @test ForwardDiff.gradient(evaluate,[2.,3.])==[24.,8.]
    @test evaluate([-2.,3.])==0.
    live[2,1]=0.;dropzeros!(live) # (2,1) -> (2)
    @test evaluate([2.,3.])==8.
    @test ForwardDiff.gradient(evaluate,[2.,3.])==[8.,0.]
end

@testset "Mechanism-bound signed shock-tube trials" begin
    path=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(path)
    make_reactor(g)=IdealGasReactor(g;temperature=1200.,pressure=40one_atm,
        mole_fractions=collect(1.:g.n_species),constraint=:constant_volume)
    reactor=make_reactor(gas)
    policy=SignedIntegerShockTubeTrials(gas,path)
    @test _validate_trial_policy(policy,reactor)===nothing
    @test _trial_data_equal(gas,policy.verified)
    @test _check_integer_trial_reactions(gas.reaction)===nothing
    @test_throws ArgumentError shocktube_ad_jacobian(make_reactor(deepcopy(gas));trial_policy=policy)

    # Same species and reaction counts do not establish numerical equivalence.
    wrong=deepcopy(gas)
    wrong.reaction.Arrhenius_coeffs[1,1]+=1.
    @test wrong.species_names==gas.species_names && wrong.n_reactions==gas.n_reactions
    @test_throws ArgumentError SignedIntegerShockTubeTrials(wrong,path)
    old=gas.thermo.nasa_high[1,1]
    try
        gas.thermo.nasa_high[1,1]+=1.
        @test_throws ArgumentError shocktube_ad_jacobian(reactor;trial_policy=policy)
    finally
        gas.thermo.nasa_high[1,1]=old
    end
    @test _validate_trial_policy(policy,reactor)===nothing

    mktempdir() do tmp
        input=YAML.load_file(path)
        input["reactions"][1]["orders"]=Dict(gas.species_names[k]=>v
            for (k,v) in enumerate(gas.reaction.reactant_stoich_coeffs[:,1]) if v!=0)
        custom=joinpath(tmp,"explicit-default-orders.yaml")
        YAML.write_file(custom,input)
        # Equal numerical orders still select Cantera's different Cn kernel.
        # Reject the explicit YAML metadata before requiring any sidecar.
        err=try SignedIntegerShockTubeTrials(gas,custom);nothing catch e;e end
        @test err isa ArgumentError
        @test occursin("explicit custom orders",sprint(showerror,err))
        delete!(input["reactions"][1],"orders")
        input["reactions"][1]["type"]="Blowers-Masel"
        YAML.write_file(custom,input)
        @test_throws ArgumentError SignedIntegerShockTubeTrials(gas,custom)

        source=joinpath(tmp,"unbound-sidecar.yaml")
        cp(path,source)
        sidecar=npzread(path*".npz")
        delete!(sidecar,"source_sha256_utf8")
        npzwrite(source*".npz",sidecar)
        @test_throws ArgumentError SignedIntegerShockTubeTrials(gas,source)
        cp(path*".npz",source*".npz";force=true)
        open(source,"a") do io;println(io,"# changed source bytes");end
        @test_throws ArgumentError SignedIntegerShockTubeTrials(gas,source)
    end

    fractional=deepcopy(gas.reaction)
    fractional.reactant_orders[1,1]=0.5
    fractional.reactant_stoich_coeffs[1,1]=0.5
    @test_throws ArgumentError _check_integer_trial_reactions(fractional)
    high_order=deepcopy(gas.reaction)
    high_order.reactant_orders[1,1]=4.
    high_order.reactant_stoich_coeffs[1,1]=4.
    @test_throws ArgumentError _check_integer_trial_reactions(high_order)
end

@testset "Signed trial Jacobian and unchanged input states" begin
    path=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(path)
    policy=SignedIntegerShockTubeTrials(gas,path)
    input=YAML.load_file(path)
    input["phases"][1]["thermo"]="Redlich-Kwong"
    for (i,species) in enumerate(input["species"])
        # Derivative fixture only; these coefficients are not an EOS fit.
        species["equation-of-state"]=Dict("model"=>"Redlich-Kwong",
            "a"=>string(1e5*(1+i/10))*" Pa*m^6*K^0.5/kmol^2",
            "b"=>string(.015+.001i)*" m^3/kmol")
    end
    model=RedlichKwongThermo(input;molecular_weights=gas.MW)
    for real in (false,true)
        reactor=real ? RedlichKwongReactor(gas,model;temperature=1200.,pressure=40one_atm,
            mole_fractions=collect(1.:gas.n_species)) :
            IdealGasReactor(gas;temperature=1200.,pressure=40one_atm,
                mole_fractions=collect(1.:gas.n_species),constraint=:constant_volume)
        jac=shocktube_ad_jacobian(reactor;trial_policy=policy)
        J=zeros(gas.n_species+1,gas.n_species+1)
        for branch in (:positive,:zero,:one_negative,:two_negative)
            u=reactor_state(reactor)
            branch!==:positive && (u[1]=branch===:zero ? 0. : -1e-8)
            branch===:two_negative && (u[2]=-1e-8)
            u[end-1]+=1-sum(view(u,1:gas.n_species))
            original=copy(u)
            jac(J,u,nothing,0.)
            full=ForwardDiff.jacobian(u) do state
                rhs=_signed_reactor_rhs(reactor,eltype(state))
                out=zero(state);rhs(out,state,nothing,0.);out
            end
            row_scale=maximum(abs.(full);dims=2).+1e-50
            @test u==original
            @test all(isfinite,J)
            @test maximum(abs.(J.-full)./row_scale)<1e-11
            @test norm(sum(J[1:end-1,:];dims=1),Inf)<1e-11*max(norm(J,Inf),1.)
            @test norm(gas.ele_matrix*(J[1:end-1,:]./gas.MW),Inf)<1e-11*max(norm(J,Inf),1.)
        end
        # This native hydrogen fixture uses argon as its inert species.
        oxidizer=real ? RedlichKwongReactor(gas,model;temperature=300.,pressure=40one_atm,
            mole_fractions=Dict("O2"=>1.,"AR"=>3.76)) :
            IdealGasReactor(gas;temperature=300.,pressure=40one_atm,
                mole_fractions=Dict("O2"=>1.,"AR"=>3.76),constraint=:constant_volume)
        u=reactor_state(oxidizer)
        shocktube_ad_jacobian(oxidizer;trial_policy=policy)(J,u,nothing,0.)
        full=ForwardDiff.jacobian(u) do state
            rhs=_signed_reactor_rhs(oxidizer,eltype(state))
            out=zero(state);rhs(out,state,nothing,0.);out
        end
        absent=findall(gas.ele_matrix[findfirst(==("H"),gas.elements),:].>0)
        active=vcat(setdiff(1:gas.n_species,absent),gas.n_species+1)
        @test all(iszero,u[absent])
        @test all(iszero,full[absent,active])
        @test all(iszero,J[absent,active])
        @test J[absent,active]==full[absent,active]
        if real
            bad_model=deepcopy(model)
            bad_model.reference.nasa_high[1,1]+=1.
            bad=RedlichKwongReactor(gas,bad_model;temperature=1200.,pressure=40one_atm,
                mole_fractions=collect(1.:gas.n_species))
            @test_throws ArgumentError shocktube_ad_jacobian(bad;trial_policy=policy)
        end
    end
end

end
