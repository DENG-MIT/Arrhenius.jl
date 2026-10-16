# Finite-difference Jacobian for the solver-internal signed trial RHS.
function signed_trial_jacobian!(J, u, rhs::SignedTrialRHS, t=0)
    n = length(u)
    size(J) == (n, n) || throw(DimensionMismatch("Jacobian must have state dimensions"))
    copyto!(rhs.jac_state, u)
    rhs(rhs.jac_base, u, nothing, t)
    relative_step = cbrt(eps(eltype(rhs.jac_state)))
    @inbounds for j in 1:n
        if j == n && rhs.reactor.energy === :isothermal
            for i in 1:n
                J[i, j] = zero(eltype(J))
            end
            continue
        end
        scale = j == n ? one(u[j]) : oftype(u[j], 1e-6)
        step = relative_step * max(abs(u[j]), scale)
        rhs.jac_state[j] = u[j] + step
        step = rhs.jac_state[j] - u[j]
        if j < n && u[j] < 0
            rhs.jac_state[j] = u[j] - step
            rhs(rhs.jac_plus, rhs.jac_state, nothing, t)
            rhs.jac_state[j] = u[j] - 2 * step
            rhs(rhs.jac_minus, rhs.jac_state, nothing, t)
            for i in 1:n
                J[i, j] = (3 * rhs.jac_base[i] - 4 * rhs.jac_plus[i] +
                           rhs.jac_minus[i]) / (2 * step)
            end
            rhs.jac_state[j] = u[j]
            continue
        end
        rhs(rhs.jac_plus, rhs.jac_state, nothing, t)
        if u[j] >= step
            rhs.jac_state[j] = u[j] - step
            rhs(rhs.jac_minus, rhs.jac_state, nothing, t)
            for i in 1:n
                J[i, j] = (rhs.jac_plus[i] - rhs.jac_minus[i]) / (2step)
            end
        else
            rhs.jac_state[j] = u[j] + 2step
            rhs(rhs.jac_minus, rhs.jac_state, nothing, t)
            for i in 1:n
                J[i, j] = (-3rhs.jac_base[i] + 4rhs.jac_plus[i] - rhs.jac_minus[i]) / (2step)
            end
        end
        rhs.jac_state[j] = u[j]
    end
    return nothing
end
