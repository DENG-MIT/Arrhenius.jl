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

end
