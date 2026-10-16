using SciMLBase
using OrdinaryDiffEqBDF

function native_surface_bdf(problem; kwargs...)
    f = ODEFunction(problem.f; jac=problem.jac, tgrad=problem.tgrad)
    ode = ODEProblem(f, problem.u0, problem.tspan, problem.p)
    # Adsorption often starts on a much shorter timescale than the final horizon.
    # Let the caller override this conservative initial step; QNDF adapts afterward.
    options = haskey(kwargs,:dt) ? kwargs : (;dt=min(1e-12,(problem.tspan[2]-problem.tspan[1])/100),kwargs...)
    result = solve(ode, QNDF(); isoutofdomain=problem.isoutofdomain, options...)
    SciMLBase.successful_retcode(result) || error("surface integration failed: $(result.retcode)")
    return result
end
