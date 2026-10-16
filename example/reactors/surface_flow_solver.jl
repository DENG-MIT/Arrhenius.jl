using SciMLBase
using OrdinaryDiffEqBDF

function native_surface_flow_bdf(problem;kwargs...)
    f = ODEFunction(problem.f;jac=problem.jac,tgrad=problem.tgrad,mass_matrix=problem.mass_matrix)
    ode = ODEProblem(f,problem.u0,problem.tspan,problem.p)
    result = solve(ode,QNDF();initializealg=SciMLBase.NoInit(),
                   isoutofdomain=problem.isoutofdomain,kwargs...)
    SciMLBase.successful_retcode(result) || error("surface flow integration failed: $(result.retcode)")
    return result
end
