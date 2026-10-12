if !isdefined(@__MODULE__, :StructuredTrialJacobian)
    include(joinpath(@__DIR__, "structured_jacobian.jl"))
end
if !isdefined(@__MODULE__, :signed_trial_jacobian!)
    include(joinpath(@__DIR__, "signed_jacobian.jl"))
end

mutable struct GuardedStructuredJacobian{A}
    analytic::A
    state_copy::Vector{Float64}
    analytic_calls::Int
    fallback_calls::Int
    preparation_rhs_calls::Int
    generation::Int
    last_used_fallback::Bool
    last_boundary_reaction::Int
end
function guarded_structured_jacobian(reactor)
    jac = structured_trial_jacobian(reactor)
    GuardedStructuredJacobian(jac, zeros(reactor.gas.n_species+1), 0, 0, 0, 0, false, 0)
end
function first_anyn_boundary(forward, reverse, reversible, C)
    for i in eachindex(forward)
        if _ambiguous_anyn_zero(forward[i], C) ||
                (reversible[i] && _ambiguous_anyn_zero(reverse[i], C))
            return i
        end
    end
    return 0
end
function (guard::GuardedStructuredJacobian)(J, state, p, t)
    jac = guard.analytic
    n = length(guard.state_copy)
    length(state)==n && size(J)==(n,n) || throw(DimensionMismatch("Jacobian/state dimensions differ"))
    eltype(state)===Float64 && eltype(J)===Float64 || throw(ArgumentError("Float64 state and output required"))
    # Validate the frozen mechanism on both branches; a boundary must not hide mutation.
    _check_mechanism!(jac)
    copyto!(guard.state_copy, state)
    jac.float_rhs(jac.float_derivative, state, nothing, t)
    guard.preparation_rhs_calls += 1
    workspace = jac.float_rhs.workspace
    capacity = 0.0
    for i in 1:n-1
        capacity += state[i]*workspace.cp_R[i]/jac.reactor.gas.MW[i]
    end
    capacity *= Arrhenius.R
    isfinite(capacity) && capacity>0 || throw(DomainError(capacity,"positive heat capacity required"))
    boundary = first_anyn_boundary(jac.float_rhs.forward_plans, jac.float_rhs.reverse_plans,
        jac.reactor.gas.reaction.is_reversible, workspace.C)
    if boundary>0
        signed_trial_jacobian!(J, state, jac.float_rhs, t)
    else
        jac(J, state, p, t)
    end
    all(isfinite,J) || throw(DomainError(J,"Jacobian is not finite"))
    isequal(state,guard.state_copy) || error("guarded Jacobian mutated input state")
    guard.generation += 1
    guard.last_used_fallback = boundary>0
    guard.last_boundary_reaction = boundary
    if boundary>0
        guard.fallback_calls += 1
    else
        guard.analytic_calls += 1
    end
    return nothing
end
