using SciMLBase, OrdinaryDiffEqSDIRK
include(joinpath(@__DIR__,"real_gas_ad_jacobian.jl"))

_shocktube_outside_domain(u,p,t) = !all(isfinite,u) || u[end] <= 0 ||
    minimum(view(u,1:length(u)-1)) < -1e-13

"Native L-stable SDIRK adapter for the stiff shock-tube examples."
function shocktube_sdirk(problem;kwargs...)
    f = ODEFunction(problem.f;jac=problem.jac,tgrad=problem.tgrad)
    ode = ODEProblem(f,problem.u0,problem.tspan,problem.p)
    solution = solve(ode,KenCarp4();isoutofdomain=_shocktube_outside_domain,kwargs...)
    SciMLBase.successful_retcode(solution) || error("shock-tube integration failed: $(solution.retcode)")
    return solution
end

"Initialize source-style adaptive stepping without imposing a terminal time stop."
function shocktube_integrator(reactor;end_time=.005,jacobian=:ad)
    isfinite(end_time) && end_time>0 || throw(ArgumentError("positive finite end time required"))
    jacobian in (:ad,:finite_difference) || throw(ArgumentError("jacobian must be :ad or :finite_difference"))
    problem = reactor_problem(reactor,(0.,end_time))
    jac = jacobian === :ad ? shocktube_ad_jacobian(reactor) : problem.jac
    f = ODEFunction(problem.f;jac,tgrad=problem.tgrad)
    ode = ODEProblem(f,problem.u0,(0.,Inf),problem.p)
    return init(ode,KenCarp4();reltol=1e-13,abstol=1e-26,
        isoutofdomain=_shocktube_outside_domain,save_everystep=false,save_start=false,
        save_end=false,dense=false,maxiters=1_000_000)
end

function _shocktube_initialize(::Val{:sdirk},reactor,end_time,jacobian,trial_policy)
    trial_policy isa ClippedShockTubeTrials ||
        throw(ArgumentError("the existing SDIRK route requires ClippedShockTubeTrials"))
    shocktube_integrator(reactor;end_time,jacobian),nothing
end
function _shocktube_initialize(::Val{S},reactor,end_time,jacobian,trial_policy) where S
    S===:qndf && throw(ArgumentError("include real_gas_qndf_solver.jl to enable the optional QNDF adapter"))
    throw(ArgumentError("solver must be :sdirk or :qndf"))
end
