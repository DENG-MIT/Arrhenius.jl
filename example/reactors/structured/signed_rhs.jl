using Arrhenius
using LinearAlgebra
using SparseArrays

include(joinpath(@__DIR__, "signed_products.jl"))

struct SignedTrialRHS{R,W,T,FP,RP}
    reactor::R
    workspace::W
    jac_state::Vector{T}
    jac_base::Vector{T}
    jac_plus::Vector{T}
    jac_minus::Vector{T}
    forward_plans::FP
    reverse_plans::RP
end

function signed_trial_rhs(
    reactor;
    scalar_type=eltype(reactor.mass_fractions),
)
    reactor.constraint === :constant_pressure ||
        throw(ArgumentError("signed trial RHS requires a constant-pressure reactor"))
    reactor.energy === :adiabatic ||
        throw(ArgumentError("signed trial RHS requires an adiabatic reactor"))
    gas = reactor.gas
    workspace = Arrhenius.ReactorWorkspace(gas, scalar_type)
    state_length = gas.n_species + 1
    jac_state = zeros(scalar_type, state_length)
    jac_base = zeros(scalar_type, state_length)
    jac_plus = zeros(scalar_type, state_length)
    jac_minus = zeros(scalar_type, state_length)
    forward_plans = signed_product_plans(
        gas.reaction.reactant_stoich_coeffs,
        gas.reaction.reactant_orders,
    )
    reverse_plans = signed_product_plans(
        gas.reaction.product_stoich_coeffs,
        gas.reaction.product_stoich_coeffs,
    )
    length(forward_plans) == gas.n_reactions ||
        error("forward signed-product plan count mismatch")
    length(reverse_plans) == gas.n_reactions ||
        error("reverse signed-product plan count mismatch")
    return SignedTrialRHS(
        reactor,
        workspace,
        jac_state,
        jac_base,
        jac_plus,
        jac_minus,
        forward_plans,
        reverse_plans,
    )
end

function (rhs::SignedTrialRHS)(du, u, p, t)
    reactor = rhs.reactor
    workspace = rhs.workspace
    gas = reactor.gas
    reaction = gas.reaction
    ns = gas.n_species
    length(u) == ns + 1 ||
        throw(DimensionMismatch("state must contain n_species + 1 entries"))
    length(du) == ns + 1 ||
        throw(DimensionMismatch("derivative must match reactor state"))

    T, P, density, inverse_mw = Arrhenius._reactor_tpρ(reactor, u)
    gas_constant = oftype(T, R)
    @inbounds for k in 1:ns
        workspace.X[k] = u[k] / (gas.MW[k] * inverse_mw)
        workspace.C[k] = u[k] * density / gas.MW[k]
    end
    cal_cp_R!(workspace.cp_R, gas, T, P, workspace.X)
    cal_h_RT!(workspace.h_mole, gas, T, P, workspace.X)
    cal_s0_R!(workspace.entropy, gas, T, P, workspace.X)
    @inbounds for k in 1:ns
        workspace.h_mole[k] *= gas_constant * T
        workspace.entropy[k] *= gas_constant
    end

    kinetics = workspace.kinetics
    Arrhenius._rate_factors!(
        reaction,
        T,
        workspace.C,
        workspace.entropy,
        workspace.h_mole,
        kinetics;
        rate_multipliers=reactor.rate_multipliers,
        pressure=P,
    )
    @inbounds for i in 1:reaction.n_reactions
        forward = signed_multiply(
            rhs.forward_plans[i], workspace.C, kinetics.kf[i],
        )
        reverse = reaction.is_reversible[i] ?
            signed_multiply(rhs.reverse_plans[i], workspace.C, kinetics.kr[i]) :
            zero(kinetics.kr[i])
        kinetics.kf[i] = forward
        kinetics.kr[i] = reverse
        kinetics.rates_of_progress[i] = forward - reverse
    end
    mul!(workspace.wdot, reaction.vk, kinetics.rates_of_progress)

    capacity = zero(T)
    energy_source = zero(T)
    @inbounds for k in 1:ns
        du[k] = workspace.wdot[k] * gas.MW[k] / density
        capacity += u[k] * workspace.cp_R[k] / gas.MW[k]
        energy_source += workspace.wdot[k] * workspace.h_mole[k]
    end
    capacity *= gas_constant
    isfinite(capacity) && capacity > zero(capacity) ||
        throw(DomainError(capacity, "positive heat capacity required"))
    du[end] = -energy_source / (density * capacity)
    return nothing
end

