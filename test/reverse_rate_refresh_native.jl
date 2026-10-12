module ReverseRateRefreshContracts

using Arrhenius, LinearAlgebra, Test

@testset "Internal reverse factors preserve public kinetics contracts" begin
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    reaction=deepcopy(gas.reaction)
    # Public queries must also provide Kc for irreversible reactions.
    reaction.is_reversible[1]=false
    plan=Arrhenius._ReversibleRatePlan(reaction)
    T=1100.;P=3one_atm
    X=mole_fractions(gas,collect(1.:gas.n_species))
    C=P/(R*T).*X
    h=cal_h_RT(gas,T,P,X).*(R*T)
    S=cal_s0_R(gas,T,P,X).*R
    workspace=KineticsWorkspace(reaction)
    partial=KineticsWorkspace(reaction)
    out=fill(-123.,gas.n_species)
    result=wdot!(out,reaction,T,C,S,h,workspace;get_rate_constants=true)
    @test keys(result)==(:forward,:reverse,:equilibrium)
    @test result.forward===workspace.kf
    @test result.reverse===workspace.kr
    @test result.equilibrium===workspace.equilibrium_constants
    @test all(==(-123.),out)
    expected=exp.((transpose(reaction.vk)*S)./R .-
        (transpose(reaction.vk)*h)./(R*T) .+ log(one_atm/(R*T)).*reaction.vk_sum)
    @test result.equilibrium≈expected rtol=3e-14
    @test all(isfinite,result.equilibrium)
    @test all(>(0),result.equilibrium)
    @test result.reverse[1]==0.
    reference=(copy(result.forward),copy(result.reverse),copy(result.equilibrium))

    fill!(partial.equilibrium_constants,NaN)
    @test Arrhenius._rate_factors!(reaction,T,C,S,h,partial,plan)===nothing
    @test partial.kf==reference[1]
    @test partial.kr==reference[2]
    @test partial.equilibrium_constants[reaction.is_reversible]==reference[3][reaction.is_reversible]
    @test isnan(partial.equilibrium_constants[1])
    # A public query always repairs every equilibrium entry in the same workspace.
    complete=wdot!(out,reaction,T,C,S,h,partial;get_rate_constants=true)
    @test complete.equilibrium==reference[3]

    cache=Arrhenius._KineticsTemperatureCache(reaction)
    @test Arrhenius._rate_factors!(reaction,T,C,S,h,partial,plan;temperature_cache=cache)===nothing
    @test cache.temperature==T
    @test cache.equilibrium==reference[3]
    @test partial.equilibrium_constants==reference[3]
    fill!(partial.equilibrium_constants,NaN)
    Arrhenius._rate_factors!(reaction,T,C,S,h,partial,plan)
    complete=wdot!(out,reaction,T,C,S,h,partial;temperature_cache=cache,get_rate_constants=true)
    @test complete.equilibrium==reference[3]
    changed_h=copy(h);changed_h[1]+=1e6
    cache.temperature=NaN
    Arrhenius._rate_factors!(reaction,T,C,S,changed_h,partial,plan;temperature_cache=cache)
    fresh=wdot!(out,reaction,T,C,S,changed_h,workspace;get_rate_constants=true)
    @test cache.temperature==T
    @test partial.kf==fresh.forward
    @test partial.kr==fresh.reverse
    @test partial.equilibrium_constants==fresh.equilibrium

    # Return aliases and flag precedence remain the original public behavior.
    both=wdot!(out,reaction,T,C,S,h,workspace;get_qdot=true,get_rate_constants=true)
    @test both.equilibrium===workspace.equilibrium_constants
    progress=wdot!(out,reaction,T,C,S,h,workspace;get_qdot=true)
    @test progress===workspace.rates_of_progress
    expected_source=reaction.vk*copy(progress)
    @test all(==(-123.),out)
    source=wdot!(out,reaction,T,C,S,h,workspace)
    @test source===out
    @test source==expected_source
end

@testset "Prepared reverse plans reject stale mechanisms without writing factors" begin
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    reaction=deepcopy(gas.reaction)
    reaction.is_reversible[1]=false
    plan=Arrhenius._ReversibleRatePlan(reaction)
    T=1000.;C=fill(.01,gas.n_species);S=collect(1.:gas.n_species);h=1e6S
    work=KineticsWorkspace(reaction)
    Arrhenius._rate_factors!(reaction,T,C,S,h,work,plan)
    snapshot=[copy(getfield(work,i)) for i in 1:fieldcount(typeof(work))]
    untouched()=all(getfield(work,i)==snapshot[i] for i in eachindex(snapshot))
    @test_throws ArgumentError Arrhenius._rate_factors!(deepcopy(reaction),T,C,S,h,work,plan)
    @test untouched()
    # Both adding and removing a reversible reaction invalidate a prepared plan.
    for index in (1,findfirst(identity,reaction.is_reversible))
        old=reaction.is_reversible[index]
        try
            reaction.is_reversible[index]=!old
            @test_throws ArgumentError Arrhenius._rate_factors!(reaction,T,C,S,h,work,plan)
            @test untouched()
            @test_throws ArgumentError Arrhenius._rate_factors!(reaction,T,C,S,h,work,plan;
                temperature_cache=Arrhenius._KineticsTemperatureCache(reaction))
            @test untouched()
            rebuilt=Arrhenius._ReversibleRatePlan(reaction)
            current=KineticsWorkspace(reaction)
            @test Arrhenius._rate_factors!(reaction,T,C,S,h,current,rebuilt)===nothing
            full=wdot!(zeros(gas.n_species),reaction,T,C,S,h,KineticsWorkspace(reaction);get_rate_constants=true)
            @test current.kf==full.forward
            @test current.kr==full.reverse
        finally
            reaction.is_reversible[index]=old
        end
    end
    @test Arrhenius._validate_reversible_plan(plan,reaction)===nothing
    push!(reaction.is_reversible,false)
    try
        @test_throws ArgumentError Arrhenius._rate_factors!(reaction,T,C,S,h,work,plan)
        @test untouched()
        @test_throws DimensionMismatch Arrhenius._ReversibleRatePlan(reaction)
    finally
        pop!(reaction.is_reversible)
    end
end

end
