using SciMLBase
using OrdinaryDiffEqBDF

function native_network_bdf(problem; kwargs...)
    f = ODEFunction(problem.f; jac=problem.jac, tgrad=problem.tgrad)
    ode = ODEProblem(f, problem.u0, problem.tspan, problem.p)
    result = solve(ode, QNDF(); isoutofdomain=problem.isoutofdomain, kwargs...)
    SciMLBase.successful_retcode(result) || error("network integration failed: $(result.retcode)")
    return result
end
