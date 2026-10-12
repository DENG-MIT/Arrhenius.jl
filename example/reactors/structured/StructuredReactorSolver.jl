module StructuredReactorSolver

using Arrhenius
using ForwardDiff
using KLU
using LinearAlgebra
using OrdinaryDiffEqBDF
using SciMLBase
using SparseArrays
import LinearSolve
const LS = LinearSolve

struct UnsupportedStructuredReactor <: Exception
    message::String
end
Base.showerror(io::IO, err::UnsupportedStructuredReactor) = print(io, err.message)

include("signed_rhs.jl")
include("signed_jacobian.jl")
include("signed_product_gradient.jl")
include("structured_jacobian.jl")
include("guarded_jacobian.jl")
include("klu_structured_solver.jl")
include("qndf_structured_adapter.jl")

function try_guarded_jacobian(problem)
    problem.f isa Arrhenius.ReactorRHS || return nothing
    eltype(problem.u0) === Float64 || return nothing
    try
        return guarded_structured_jacobian(problem.f.reactor)
    catch err
        err isa UnsupportedStructuredReactor && return nothing
        rethrow()
    end
end

function try_structured_adapter(problem)
    guard = try_guarded_jacobian(problem)
    guard === nothing && return nothing
    try
        return qndf_structured_adapter(guard)
    catch err
        err isa UnsupportedStructuredReactor && return nothing
        rethrow()
    end
end

export UnsupportedStructuredReactor, try_structured_adapter
export adapter_precs, adapter_linsolve, adapter_report
export signed_trial_rhs, signed_trial_jacobian!

end
