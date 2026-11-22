"""
Analytic (and code-generated) structured Jacobian, sparse KLU/bordered
solve, and the `native_bdf`/`native_rosenbrock23` integrators for
constant-pressure, adiabatic ideal-gas reactors.

Loading all of `KLU`, `LinearSolve`, `OrdinaryDiffEqBDF`,
`OrdinaryDiffEqDifferentiation`, `OrdinaryDiffEqNonlinearSolve`,
`OrdinaryDiffEqRosenbrock`, and `SciMLBase` alongside `Arrhenius` activates
this extension. Deliberate behavior change from the old
`example/reactors/ode_solver.jl` script: `native_bdf`/`native_rosenbrock23`
now require the *entire* weakdep list to be loaded, including
`OrdinaryDiffEqRosenbrock`, even for the plain dense-Jacobian QNDF path
that previously needed only `SciMLBase`+`OrdinaryDiffEqBDF`. This trades a
small amount of optionality for a single, Julia-native activation
mechanism instead of a hand-rolled `Base.find_package`/`Base.pkgversion`
gate.
"""
module ArrheniusStructuredReactorExt

using Arrhenius
using ForwardDiff
using KLU
using LinearAlgebra
using LinearSolve
using OrdinaryDiffEqBDF
using OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqRosenbrock
using SciMLBase
using SparseArrays
const LS = LinearSolve

struct UnsupportedStructuredReactor <: Exception
    message::String
end
Base.showerror(io::IO, err::UnsupportedStructuredReactor) = print(io, err.message)

include("signed_rhs.jl")
include("signed_jacobian.jl")
include("signed_product_gradient.jl")
include("structured_jacobian.jl")
include("generated_jacobian.jl")
include("guarded_jacobian.jl")
include("klu_structured_solver.jl")
include("qndf_structured_adapter.jl")

"""
`builder` selects which `concentration_core`-filling implementation the
guarded Jacobian uses: `structured_trial_jacobian` (default, hand-written
runtime loop) or `generated_structured_jacobian` (`generated_jacobian.jl`.
 Everything downstream (the ambiguous-zero-crossing fallback here,
 and the KLU/QNDF wiring in `qndf_structured_adapter`)
 is unaffected by this choice.
"""
function try_guarded_jacobian(problem; builder=structured_trial_jacobian)
    problem.f isa Arrhenius.ReactorRHS || return nothing
    eltype(problem.u0) === Float64 || return nothing
    try
        return guarded_structured_jacobian(problem.f.reactor; builder)
    catch err
        err isa UnsupportedStructuredReactor && return nothing
        rethrow()
    end
end

function try_structured_adapter(problem; builder=structured_trial_jacobian)
    guard = try_guarded_jacobian(problem; builder)
    guard === nothing && return nothing
    try
        return qndf_structured_adapter(guard)
    catch err
        err isa UnsupportedStructuredReactor && return nothing
        rethrow()
    end
end

export UnsupportedStructuredReactor, try_structured_adapter, try_guarded_jacobian
export adapter_precs, adapter_linsolve, adapter_report
export signed_trial_rhs, signed_trial_jacobian!
export generated_structured_jacobian

function Arrhenius.native_bdf(problem; linear_solver=:auto, kwargs...)
    linear_solver in (:auto, :dense, :generated) ||
        throw(ArgumentError("linear_solver must be :auto, :dense, or :generated"))
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
    if linear_solver === :dense || linear_solver === :generated
        builder = linear_solver === :generated ?
            generated_structured_jacobian :
            structured_trial_jacobian
        guard = try_guarded_jacobian(problem; builder)
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
        adapter = try_structured_adapter(problem)
        if adapter !== nothing
            integrator_ref = Ref{Any}(nothing)
            precs = adapter_precs(adapter, integrator_ref)
            linsolve = adapter_linsolve(adapter, integrator_ref)
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
            _lifecycle_valid(adapter) ||
                error("structured reactor solver cache lifecycle failed")
            return solution
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

function Arrhenius.native_rosenbrock23(problem; linear_solver=:analytic, kwargs...)
    linear_solver in (:analytic, :generated) ||
        throw(ArgumentError("linear_solver must be :analytic or :generated"))
    builder = linear_solver === :generated ?
        generated_structured_jacobian :
        structured_trial_jacobian
    guard = try_guarded_jacobian(problem; builder)
    guard === nothing && throw(ArgumentError(
        "native_rosenbrock23 requires a reactor in the structured Jacobian's supported scope"))
    f = SciMLBase.ODEFunction(guard.analytic.float_rhs; jac=guard, tgrad=problem.tgrad)
    ode = SciMLBase.ODEProblem(f, problem.u0, problem.tspan, problem.p)
    outside_domain = (u, p, t) -> !all(isfinite, u) || u[end] <= 0
    solution = SciMLBase.solve(ode, OrdinaryDiffEqRosenbrock.Rosenbrock23();
        isoutofdomain=outside_domain, kwargs...)
    SciMLBase.successful_retcode(solution) ||
        error("reactor integration failed: $(solution.retcode)")
    return solution
end

end
