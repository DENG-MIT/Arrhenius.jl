# Usage: julia --project=SOLVER_ENV real_gas_jacobians.jl REFERENCE_DIR
# Prepare data with real_gas_reactor_cases.py. Requires ForwardDiff.
# Both mechanisms include nonzero third-body/falloff rates; h2-plog also checks
# derivatives of true EOS pressure in PLOG. These are numerical test states.
using Arrhenius, ForwardDiff, LinearAlgebra, Test
include(joinpath(@__DIR__,"..","example","reactors","real_gas_ad_jacobian.jl"))

function uncached_reactor_derivative(reactor,u)
    T=eltype(u)
    work=reactor isa RedlichKwongReactor ? RealGasKineticsWorkspace(reactor.gas,reactor.model,T) :
        Arrhenius.ReactorWorkspace(reactor.gas,T)
    rhs=Arrhenius.ReactorRHS(reactor,work,copy(u),zero(u),zero(u),zero(u))
    du=zero(u)
    rhs(du,u,nothing,0.)
    return du
end

directory=only(ARGS)
@testset "RK/IG cached Jacobians including PLOG" begin
    for name in ("dodecane","h2-plog"), phase in ("RK","IG"), temperature in (760.,1000.,1600.)
        gas=CreateSolution(joinpath(directory,name*"_IG.yaml"))
        model=RedlichKwongThermo(joinpath(directory,name*"_RK.yaml"))
        composition=collect(1.:gas.n_species)
        reactor=phase=="RK" ? RedlichKwongReactor(gas,model;temperature,pressure=40one_atm,mole_fractions=composition) :
            IdealGasReactor(gas;temperature,pressure=40one_atm,mole_fractions=composition,constraint=:constant_volume)
        u=reactor_state(reactor)
        full=ForwardDiff.jacobian(x -> uncached_reactor_derivative(reactor,x),u)
        J=similar(full)
        shocktube_ad_jacobian(reactor)(J,u,nothing,0.)
        error=maximum(abs.(J-full)./max.(maximum(abs.(full),dims=2),1e-20))
        @test error<2e-12
        # Newton iterations can use slightly negative trial species. The
        # Jacobian must differentiate the same clipped-concentration branch
        # as the RHS, rather than inserting a positive-species derivative.
        negative=copy(u)
        negative[2:3:end-1].=-1e-20
        negative_full=ForwardDiff.jacobian(x -> uncached_reactor_derivative(reactor,x),negative)
        negative_J=similar(J)
        shocktube_ad_jacobian(reactor)(negative_J,negative,nothing,0.)
        @test maximum(abs.(negative_J-negative_full)./max.(maximum(abs.(negative_full),dims=2),1e-20))<2e-12
        direction=sin.(collect(1.:length(u))).*max.(abs.(u),1e-6)
        rhs=reactor_rhs(reactor)
        plus,minus=zero(u),zero(u)
        step=1e-5
        if temperature==1000.
            # A central difference crosses the NASA7 coefficient switch. Use
            # a second-order one-sided difference in the selected lower region.
            direction[end]=-abs(direction[end])
            base=zero(u); rhs(base,u,nothing,0.)
            rhs(plus,u.+step.*direction,nothing,0.)
            rhs(minus,u.+2step.*direction,nothing,0.)
            numerical=(-3base.+4plus.-minus)./(2step)
        else
            rhs(plus,u.+step.*direction,nothing,0.)
            rhs(minus,u.-step.*direction,nothing,0.)
            numerical=(plus.-minus)./(2step)
        end
        @test maximum(abs.(J*direction-numerical)./max.(abs.(numerical),1))<2e-5
        println(name," ",phase," ",temperature," K: row-relative AD difference ",error)
    end
end
