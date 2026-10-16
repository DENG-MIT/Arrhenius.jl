# Usage: julia --project=SOLVER_ENV example/reactors/non_ideal_shock_tube.jl MECHANISM_DIR
#        [--smoke] [--finite-difference] [--qndf]
# The caller-owned solver environment needs SciMLBase, OrdinaryDiffEqSDIRK and ForwardDiff.
# The optional --qndf route additionally needs OrdinaryDiffEqBDF. From Julia,
# include real_gas_qndf_solver.jl after this file before selecting solver=:qndf.
# MECHANISM_DIR contains dodecane_IG.yaml, its native kinetic sidecar, and
# dodecane_RK.yaml. validation/real_gas_reactor_cases.py prepares the paired input
# files and independent Cantera references. No Cantera call occurs in this file.
# Source: https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/reactors/non_ideal_shock_tube.py
using Arrhenius, LinearAlgebra, Printf
include(joinpath(@__DIR__,"real_gas_ode_solver.jl"))
BLAS.set_num_threads(1)

const SHOCKTUBE_TEMPERATURES = [1250,1170,1120,1080,1040,1010,990,970,950,930,910,880,850,820,790,760]
const SHOCKTUBE_COMPOSITION = Dict("c12h26"=>1.,"o2"=>18.5,"n2"=>69.56)

"""
    shocktube_ignition(reactor; end_time=0.005, save_stride=20, jacobian=:ad,
                       solver=:sdirk, trial_policy=ClippedShockTubeTrials())

Integrate a constant-volume adiabatic reactor until its first accepted step at
or beyond `end_time`, retaining every `save_stride`th accepted state. This is
Cantera's published shock-tube sampling convention, including final overshoot.
Choose the cached AD Jacobian or `jacobian=:finite_difference` with SDIRK.
The optional `solver=:qndf` route requires its included adapter, AD and an explicit
mechanism-bound `SignedIntegerShockTubeTrials` policy.
Returns sampled state vectors `[Y; T]`, OH-peak delay, endpoint and solver counts.
"""
function shocktube_ignition(reactor;end_time=.005,save_stride=20,jacobian=:ad,
        solver=:sdirk,trial_policy::ShockTubeTrialPolicy=ClippedShockTubeTrials())
    save_stride isa Integer && save_stride>0 || throw(ArgumentError("positive integer save stride required"))
    oh = findfirst(==("oh"),reactor.gas.species_names)
    isnothing(oh) && throw(ArgumentError("the source mechanism must contain species named oh"))
    integrator,scale = _shocktube_initialize(Val(solver),reactor,end_time,jacobian,trial_policy)
    return _shocktube_collect(integrator,oh,reactor.density,end_time,save_stride,scale)
end

# Initialization can choose distinct solver/Jacobian types. Dispatch once here
# so the accepted-step loop specializes on the concrete integrator type.
_shocktube_physical_state(u,::Nothing)=copy(u)
_shocktube_physical_state(u,scale::AbstractVector)=u.*scale

function _shocktube_collect(integrator,oh,density,end_time,save_stride,scale=nothing)
    times,states = Float64[],Vector{Float64}[]
    counter = 0
    while integrator.t < end_time
        step!(integrator)
        code = SciMLBase.check_error(integrator)
        code in (ReturnCode.Success,ReturnCode.Default) || error("integration failed: $code")
        counter += 1
        if counter % save_stride == 0
            push!(times,integrator.t)
            push!(states,_shocktube_physical_state(integrator.u,scale))
        end
    end
    isempty(times) && error("too few accepted steps for the source sampling stride")
    delay = times[argmax([state[oh] for state in states])]
    return (;delay,times,states,steps=counter,
        final_time=integrator.t,final_state=_shocktube_physical_state(integrator.u,scale),density,
        rhs_evaluations=integrator.stats.nf,jacobian_evaluations=integrator.stats.njacs,
        accepted_steps=integrator.stats.naccept,rejected_steps=integrator.stats.nreject,
        linear_solves=integrator.stats.nsolve,matrix_updates=integrator.stats.nw,
        nonlinear_iterations=integrator.stats.nnonliniter,
        nonlinear_convergence_failures=integrator.stats.nnonlinconvfail)
end

"Run the published 1000 K pair, then the 16 RK cases followed by the 16 ideal-gas cases."
function shocktube_calculations(gas,model;smoke=false,jacobian=:ad,solver=:sdirk,
        trial_policy::ShockTubeTrialPolicy=ClippedShockTubeTrials())
    cases = [("RK",1000),("IG",1000)]
    if !smoke
        append!(cases,[(phase,T) for phase in ("RK","IG") for T in SHOCKTUBE_TEMPERATURES])
    end
    return map(cases) do (phase,temperature)
        reactor = phase == "RK" ?
            RedlichKwongReactor(gas,model;temperature,pressure=40one_atm,mole_fractions=SHOCKTUBE_COMPOSITION) :
            IdealGasReactor(gas;temperature,pressure=40one_atm,mole_fractions=SHOCKTUBE_COMPOSITION,constraint=:constant_volume)
        merge((;phase,temperature),shocktube_ignition(reactor;jacobian,solver,trial_policy))
    end
end

"Load the prepared mechanism and run all 34 published calculations (`sweep=false` selects the 1000 K pair)."
function non_ideal_shock_tube(directory;sweep=true,jacobian=:ad,solver=:sdirk)
    gas = CreateSolution(joinpath(directory,"dodecane_IG.yaml"))
    model = RedlichKwongThermo(joinpath(directory,"dodecane_RK.yaml"))
    trial_policy = solver===:qndf ? SignedIntegerShockTubeTrials(gas,joinpath(directory,"dodecane_IG.yaml")) : ClippedShockTubeTrials()
    return shocktube_calculations(gas,model;smoke=!sweep,jacobian,solver,trial_policy)
end

if abspath(PROGRAM_FILE) == @__FILE__
    !isempty(ARGS) && all(arg -> arg in ("--smoke","--sweep","--finite-difference","--qndf"),ARGS[2:end]) ||
        error("supply MECHANISM_DIR [--smoke] [--finite-difference] [--qndf]")
    "--qndf" in ARGS && include(joinpath(@__DIR__,"real_gas_qndf_solver.jl"))
    results = non_ideal_shock_tube(ARGS[1];sweep=!("--smoke" in ARGS),
        jacobian="--finite-difference" in ARGS ? :finite_difference : :ad,
        solver="--qndf" in ARGS ? :qndf : :sdirk)
    for result in results
        @printf("%s T=%4d K: ignition delay %.9g s\n",result.phase,result.temperature,result.delay)
    end
    @printf("1000 K ideal-gas error: %.4g%%\n",100*(results[2].delay/results[1].delay-1))
end
