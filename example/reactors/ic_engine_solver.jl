using SciMLBase, OrdinaryDiffEqSDIRK

"Native L-stable SDIRK integration for the periodically switched engine network."
function native_engine_sdirk(problem;max_temperature_change=20.0,dtmax=1/(360*50.0),kwargs...)
    temperature_index = problem.f.model.network.nodes.cylinder.initial.gas.n_species+1
    previous_temperature = Ref(problem.u0[temperature_index])
    outside = (u,p,t)->problem.isoutofdomain(u,p,t) ||
        abs(u[temperature_index]-previous_temperature[]) > max_temperature_change
    function accepted!(integrator)
        previous_temperature[] = integrator.u[temperature_index]
        SciMLBase.u_modified!(integrator,false)
    end
    limiter = DiscreteCallback((u,t,integrator)->true,accepted!;save_positions=(false,false))
    ode = ODEProblem(ODEFunction(problem.f;jac=problem.jac,tgrad=problem.tgrad),
                     problem.u0,problem.tspan,problem.p)
    solution = solve(ode,KenCarp4();isoutofdomain=outside,callback=limiter,dtmax,kwargs...)
    SciMLBase.successful_retcode(solution) || error("engine integration failed: $(solution.retcode)")
    return solution
end
