# Usage: julia --project=SOLVER_ENV real_gas_adapter_contracts.jl MECHANISM_DIR
# The caller environment needs SciMLBase, OrdinaryDiffEqSDIRK,
# OrdinaryDiffEqBDF and ForwardDiff. Prepare the dodecane files using
# real_gas_reactor_cases.py. This validation performs no Cantera calls.
using Arrhenius, LinearAlgebra, Test
include(joinpath(@__DIR__,"..","example","reactors","non_ideal_shock_tube.jl"))
include(joinpath(@__DIR__,"..","example","reactors","real_gas_qndf_solver.jl"))

function adapter_contracts(directory)
    path=joinpath(directory,"dodecane_IG.yaml")
    gas=CreateSolution(path)
    model=RedlichKwongThermo(joinpath(directory,"dodecane_RK.yaml"))
    policy=SignedIntegerShockTubeTrials(gas,path)
    @testset "Caller-owned shock-tube solver contracts" begin
        for real in (false,true)
            reactor=real ? RedlichKwongReactor(gas,model;temperature=1000.,pressure=40one_atm,
                mole_fractions=SHOCKTUBE_COMPOSITION) :
                IdealGasReactor(gas;temperature=1000.,pressure=40one_atm,
                    mole_fractions=SHOCKTUBE_COMPOSITION,constraint=:constant_volume)
            u=reactor_state(reactor)
            # Existing calls selecting only finite differences return an
            # untransformed integrator with the original physical state.
            legacy=shocktube_integrator(reactor;jacobian=:finite_difference)
            @test legacy.u==u
            @test legacy.sol.prob.tspan==(0.,Inf)
            @test legacy.opts.reltol==1e-13
            @test legacy.opts.abstol==1e-26

            integrator,scale=shocktube_qndf_integrator(reactor;trial_policy=policy)
            @test scale==vcat(gas.MW./reactor.density,1000.)
            @test integrator.u==u./scale
            @test integrator.sol.prob.tspan==(0.,Inf)
            @test integrator.opts.reltol==1e-9
            @test all(==(1e-19),integrator.opts.abstol[1:end-1])
            @test integrator.opts.abstol[end]==1e-6/1000.
            @test integrator.alg.nlsolve.κ==1//100
            @test _shocktube_physical_state(integrator.u,scale)≈u rtol=2e-16
            @test _shocktube_physical_state(u,nothing)==u
            @test _shocktube_physical_state(u,nothing)!==u
            domain=_ScaledShockTubeDomain(scale)
            z=copy(integrator.u)
            @test !domain(z,nothing,0.)
            z[1]=-2e-13/scale[1]
            @test domain(z,nothing,0.)
            z[1]=-5e-14/scale[1]
            @test !domain(z,nothing,0.)
            @test_throws ArgumentError shocktube_ignition(reactor;solver=:unknown)
            @test_throws ArgumentError shocktube_ignition(reactor;solver=:qndf)
            @test_throws ArgumentError shocktube_ignition(reactor;solver=:qndf,
                jacobian=:finite_difference,trial_policy=policy)
            @test_throws ArgumentError shocktube_ignition(reactor;solver=:sdirk,trial_policy=policy)
            @test_throws ArgumentError shocktube_qndf_integrator(reactor;trial_policy=policy,end_time=0.)
            # True zeros preserve the complete absent-H/C invariant block.
            air=real ? RedlichKwongReactor(gas,model;temperature=300.,pressure=40one_atm,
                mole_fractions=Dict("o2"=>1.,"n2"=>3.76)) :
                IdealGasReactor(gas;temperature=300.,pressure=40one_atm,
                    mole_fractions=Dict("o2"=>1.,"n2"=>3.76),constraint=:constant_volume)
            state=reactor_state(air)
            J=zeros(length(state),length(state))
            shocktube_ad_jacobian(air;trial_policy=policy)(J,state,nothing,0.)
            full=ForwardDiff.jacobian(state) do trial
                rhs=_signed_reactor_rhs(air,eltype(trial));out=zero(trial)
                rhs(out,trial,nothing,0.);out
            end
            elements=findall(e -> e in ("H","C"),gas.elements)
            absent=findall(vec(sum(gas.ele_matrix[elements,:];dims=1)).>0)
            active=vcat(setdiff(1:gas.n_species,absent),gas.n_species+1)
            @test length(elements)==2
            @test all(iszero,state[absent])
            @test all(iszero,full[absent,active])
            @test all(iszero,J[absent,active])

        end
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("supply the prepared MECHANISM_DIR")
    adapter_contracts(ARGS[1])
end
