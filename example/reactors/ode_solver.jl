using SciMLBase
using OrdinaryDiffEqBDF
using LinearAlgebra

# The optional environment pins the cache interfaces used by the structured
# solver. Other environments retain the ordinary dense QNDF integration.
if !isdefined(@__MODULE__, :StructuredReactorSolver) &&
        all(name -> Base.find_package(name) !== nothing,
            ("ForwardDiff", "LinearSolve", "KLU", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve"))
    import ForwardDiff, LinearSolve, KLU, OrdinaryDiffEqDifferentiation, OrdinaryDiffEqNonlinearSolve
    if Base.pkgversion(ForwardDiff) == v"1.4.6" && Base.pkgversion(SciMLBase) == v"2.155.2" &&
            Base.pkgversion(OrdinaryDiffEqBDF) == v"1.26.0" &&
            Base.pkgversion(OrdinaryDiffEqDifferentiation) == v"2.9.0" &&
            Base.pkgversion(OrdinaryDiffEqNonlinearSolve) == v"1.28.0" &&
            Base.pkgversion(LinearSolve) == v"3.87.0" && Base.pkgversion(KLU) == v"0.6.0"
        include(joinpath(@__DIR__, "structured", "StructuredReactorSolver.jl"))
    end
end

"""
Integrate an Arrhenius reactor problem with Julia's adaptive QNDF solver.

The optional reactor environment enables an analytic Jacobian and sparse
bordered solve for supported constant-pressure, adiabatic ideal-gas reactors.
Other reactors use the dense Jacobian supplied by `reactor_problem`.

Set `linear_solver=:dense` to use the guarded analytic Jacobian and signed RHS
with QNDF's ordinary dense linear solver for supported reactors. The default
`linear_solver=:auto` retains automatic structured selection.
"""
function native_bdf(problem; linear_solver=:auto, kwargs...)
    linear_solver in (:auto, :dense) ||
        throw(ArgumentError("linear_solver must be :auto or :dense"))
    if problem.f isa Arrhenius.PlasmaEnergyRHS
        f = SciMLBase.ODEFunction(problem.f; jac=problem.jac, tgrad=problem.tgrad,
            jac_prototype=zeros(length(problem.u0), length(problem.u0)))
        ode = SciMLBase.ODEProblem(f, problem.u0, problem.tspan, problem.p)
        # The state is mass, total enthalpy, then signed mass fractions.
        # Recover temperature from enthalpy; negative enthalpy is valid.
        outside_domain = function (u, p, t)
            (!all(isfinite, u) || u[1] <= 0) && return true
            try
                props = Arrhenius.reactor_properties(problem.f, u)
                return !(isfinite(props.T) && props.T > 0 &&
                         isfinite(props.rho) && props.rho > 0)
            catch error
                (error isa DomainError || error isa ErrorException) || rethrow()
                return true
            end
        end
        solution = SciMLBase.solve(ode, OrdinaryDiffEqBDF.QNDF();
            isoutofdomain=outside_domain, kwargs...)
        SciMLBase.successful_retcode(solution) ||
            error("plasma energy integration failed: $(solution.retcode)")
        return solution
    end
    if problem.f isa Arrhenius.PlasmaRHS
        f = SciMLBase.ODEFunction(problem.f; jac=problem.jac, tgrad=problem.tgrad,
            jac_prototype=zeros(length(problem.u0), length(problem.u0)))
        ode = SciMLBase.ODEProblem(f, problem.u0, problem.tspan, problem.p)
        outside_domain = (u, p, t) -> !all(isfinite, u) ||
            dot(problem.f.density_weights, u) <= 0
        solution = SciMLBase.solve(ode, OrdinaryDiffEqBDF.QNDF();
            isoutofdomain=outside_domain, kwargs...)
        SciMLBase.successful_retcode(solution) ||
            error("plasma integration failed: $(solution.retcode)")
        return solution
    end
    if isdefined(@__MODULE__, :StructuredReactorSolver)
        if linear_solver === :dense
            guard = StructuredReactorSolver.try_guarded_jacobian(problem)
            if guard !== nothing
                f = SciMLBase.ODEFunction(
                    guard.analytic.float_rhs; jac=guard, tgrad=problem.tgrad)
                ode = SciMLBase.ODEProblem(f, problem.u0, problem.tspan, problem.p)
                outside_domain = (u, p, t) -> !all(isfinite, u) || u[end] <= 0
                solution = SciMLBase.solve(
                    ode, OrdinaryDiffEqBDF.QNDF(); isoutofdomain=outside_domain, kwargs...)
                SciMLBase.successful_retcode(solution) ||
                    error("reactor integration failed: $(solution.retcode)")
                return solution
            end
        else
            adapter = StructuredReactorSolver.try_structured_adapter(problem)
            if adapter !== nothing
                integrator_ref = Ref{Any}(nothing)
                precs = StructuredReactorSolver.adapter_precs(adapter, integrator_ref)
                linsolve = StructuredReactorSolver.adapter_linsolve(adapter, integrator_ref)
                f = SciMLBase.ODEFunction(
                    adapter.guard.analytic.float_rhs; jac=adapter, tgrad=problem.tgrad)
                ode = SciMLBase.ODEProblem(f, problem.u0, problem.tspan, problem.p)
                outside_domain = (u, p, t) -> !all(isfinite, u) || u[end] <= 0
                integrator = SciMLBase.init(
                    ode, OrdinaryDiffEqBDF.QNDF(precs=precs, linsolve=linsolve);
                    isoutofdomain=outside_domain, kwargs...)
                integrator_ref[] = integrator
                solution = SciMLBase.solve!(integrator)
                SciMLBase.successful_retcode(solution) ||
                    error("reactor integration failed: $(solution.retcode)")
                StructuredReactorSolver._lifecycle_valid(adapter) ||
                    error("structured reactor solver cache lifecycle failed")
                return solution
            end
        end
    end
    f = SciMLBase.ODEFunction(problem.f; jac=problem.jac, tgrad=problem.tgrad)
    ode = SciMLBase.ODEProblem(f, problem.u0, problem.tspan, problem.p)
    outside_domain = (u, p, t) -> !all(isfinite, u) || u[end] <= 0 ||
        minimum(view(u, 1:length(u)-1)) < -1e-13
    solution = SciMLBase.solve(ode, OrdinaryDiffEqBDF.QNDF(); isoutofdomain=outside_domain, kwargs...)
    SciMLBase.successful_retcode(solution) || error("reactor integration failed: $(solution.retcode)")
    return solution
end
