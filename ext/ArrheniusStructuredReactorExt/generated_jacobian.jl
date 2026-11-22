using Arrhenius

if !isdefined(@__MODULE__, :StructuredTrialJacobian)
    include(joinpath(@__DIR__, "structured_jacobian.jl"))
end

# Mechanism-specialized replacement for `_fill_concentration_core!`
# (`structured_jacobian.jl`), ported from the LLMImprovments fork's
# `@generated`-function Jacobian generator (`src/Jacobian.jl`). Builds one
# unrolled, branch-free block of code per reaction instead of a runtime loop
# over `StructuredReactionPlan`/`SignedProductPlan` -- the per-reaction `kind`
# classification (1/2/3-way closed-form product vs. the general fractional-
# order path) is resolved once, here, at code-generation time, so the
# generated code for a "kind 1/2/3" reaction contains no runtime kind
# dispatch at all, and integer exponents become literal Julia integers (so
# `C[i]^2` compiles to `C[i]*C[i]` via `Base.literal_pow`, not a generic
# floating-point power call).
#
# Everything else about `StructuredTrialJacobian` -- its struct, its
# `(jac)(J,state,p,t)` callable, `_check_mechanism!`, and all of
# `guarded_jacobian.jl`/`klu_structured_solver.jl`/`qndf_structured_adapter.jl`
# -- is untouched; `generated_structured_jacobian(reactor)` builds an ordinary
# `structured_trial_jacobian(reactor)` and only replaces its `fill_core!`
# field.
#
# Scope: reactions whose forward/reverse `SignedProductPlan.kind` is 1, 2, or
# 3 (closed-form, Cantera-compatible-negative-concentration-safe product and
# gradient -- see `signed_products.jl`/`signed_product_gradient.jl`) are
# fully unrolled. `kind == 4` (fractional orders, or more than 3 total
# repeated factors) reactions fall back, per reaction, to calling the
# existing `signed_product_gradient!` directly from within the generated
# function -- still one unrolled block per reaction, just not a closed-form
# one for that reaction's product term. Collider/falloff-pressure
# derivatives (`_collider_derivative`, reusing `src/FlameDerivatives.jl`'s
# `_flame_base_rate`/`_flame_falloff_derivative`) are likewise called as-is
# rather than re-derived symbolically: they apply only to the much smaller
# three-body/falloff subset of reactions, so inlining them is not where the
# generated path's benefit comes from, and re-deriving Troe/Lindemann
# derivative formulas symbolically would be a correctness risk for no
# measured gain.

const GENERATED_STRUCTURED_CORE_BODIES = Dict{UInt64,Expr}()

struct GeneratedStructuredKernel{ID} end

"""
    product_expr(Csym, indices, orders)

`∏ C[idx]^order` over `(idx, order)` in `zip(indices, orders)`, as an `Expr`.
Integer-valued `order`s (guaranteed by the structured Jacobian's own scope
check, `_validate_structured_scope`) become literal `Int` exponents, so
`Base.literal_pow` turns `C[idx]^1`/`C[idx]^2`/`C[idx]^3` into plain
multiplication at compile time -- the same effect the LLMImprovments fork's
`k_Pi_C_rectantOrder!` gets by hand-coding `order==1 ? Ci : order==2 ? ...`.
"""
function product_expr(Csym, indices, orders)
    expr = :(one(eltype($Csym)))
    for (idx, order) in zip(indices, orders)
        expr = :($expr * $Csym[$idx]^$(Int(order)))
    end
    return expr
end

"""
    product_derivative_expr(Csym, indices, orders, derivative_position)

`∂/∂C[indices[derivative_position]]` of `product_expr(Csym, indices, orders)`,
via the power rule, as an `Expr`.
"""
function product_derivative_expr(Csym, indices, orders, derivative_position::Int)
    expr = :($(Int(orders[derivative_position])))
    for (k, (idx, order)) in enumerate(zip(indices, orders))
        exponent = k == derivative_position ? Int(order) - 1 : Int(order)
        expr = :($expr * $Csym[$idx]^$exponent)
    end
    return expr
end

"""
    kind_zero_guard_expr(Csym, repeated)

The Cantera-compatible "≥2 negative factors -> zero" guard from
`signed_multiply`/`signed_product_gradient!` (`kind` 2 or 3), as an `Expr`
over the plan's `repeated` (with-multiplicity) species list. `nothing` for
`kind == 1` (single factor, no zero-suppression -- a single negative factor
stays signed per Cantera's rule, exactly as `signed_multiply`'s `kind==1`
branch has no such check at all).
"""
function kind_zero_guard_expr(Csym, repeated)
    length(repeated) <= 1 && return nothing
    terms = [:($Csym[$r] < 0) for r in repeated]
    return reduce((a, b) -> :($a && $b), terms)
end

"""
    closed_form_contributions!(gradient_terms, plan, rate_sym, to_dependency)

For a `kind in (1,2,3)` `SignedProductPlan`, push one `Expr` per dependency
(keyed by its position in the owning reaction's unified `dependencies` list,
via `to_dependency`) into `gradient_terms[dependency_position]`, guarded by
`kind_zero_guard_expr` where applicable. `rate_sym` is the (already-signed)
rate expression -- `:kf_i` for the forward/reactant side, `:(-kr_i)` for the
reverse/product side.
"""
function closed_form_contributions!(gradient_terms, plan, rate_sym, to_dependency)
    guard = kind_zero_guard_expr(:C, plan.repeated)
    for (local_index, dependency_position) in enumerate(to_dependency)
        derivative = product_derivative_expr(:C, plan.indices, plan.orders, local_index)
        contribution = :($rate_sym * $derivative)
        guarded = guard === nothing ? contribution : :($guard ? zero($rate_sym) : $contribution)
        push!(gradient_terms[dependency_position], guarded)
    end
    return nothing
end

"""
    fallback_contributions!(statements, gradient_terms, plan, rate_sym, to_dependency, tmp)

For a `kind == 4` `SignedProductPlan` (fractional orders, or more than 3
total repeated factors), call the existing generic `signed_product_gradient!`
for this one reaction/side -- reusing its exact (already Cantera-compatible)
math rather than re-deriving it symbolically -- and push the resulting
per-dependency values into `gradient_terms`. Also reproduces the ambiguous
zero-crossing check `_ambiguous_anyn_zero` already performs in the generic
path, since that check only ever fires for `kind == 4`.
"""
function fallback_contributions!(statements, gradient_terms, plan, rate_sym, to_dependency, tmp, i::Int)
    push!(statements, :(_ambiguous_anyn_zero($plan, C) && throw(DomainError(
        $i, "ambiguous AnyN order-one derivative at a zero concentration"))))
    push!(statements, :(signed_product_gradient!($tmp, $plan, C, $rate_sym)))
    for (local_index, dependency_position) in enumerate(to_dependency)
        push!(gradient_terms[dependency_position], :($tmp[$local_index]))
    end
    return nothing
end

function reaction_block(i::Int, plan, forward_plan, reverse_plan, is_reversible::Bool)
    ndeps = length(plan.dependencies)
    gradient_terms = [Expr[] for _ in 1:ndeps]
    statements = Expr[]

    push!(statements, :(kf_i = kinetics.kf[$i]))
    if forward_plan.kind in (1, 2, 3)
        closed_form_contributions!(gradient_terms, forward_plan, :kf_i, plan.forward_to_dependency)
    else
        push!(statements, :(fg_tmp = jac.forward_gradients[$i]))
        fallback_contributions!(statements, gradient_terms, forward_plan, :kf_i, plan.forward_to_dependency, :fg_tmp, i)
    end

    if is_reversible
        push!(statements, :(kr_i = kinetics.kr[$i]))
        if reverse_plan.kind in (1, 2, 3)
            closed_form_contributions!(gradient_terms, reverse_plan, :(-kr_i), plan.reverse_to_dependency)
        else
            push!(statements, :(rg_tmp = jac.reverse_gradients[$i]))
            fallback_contributions!(statements, gradient_terms, reverse_plan, :(-kr_i), plan.reverse_to_dependency, :rg_tmp, i)
        end
    end

    if plan.collider_kind != 0
        push!(statements, :(collider_i = _collider_derivative(
            jac, $i, $plan, T, C, workspace.h_mole, kinetics)))
        for (local_index, dependency_position) in enumerate(plan.efficiency_to_dependency)
            deviation = plan.efficiency_deviations[local_index]
            push!(gradient_terms[dependency_position], :(collider_i * $deviation))
        end
    end

    noutputs = length(plan.output_rows)
    for d in 1:ndeps
        isempty(gradient_terms[d]) && continue
        total = length(gradient_terms[d]) == 1 ? gradient_terms[d][1] :
                Expr(:call, :+, gradient_terms[d]...)
        gradient_symbol = Symbol(:grad_, i, :_, d)
        push!(statements, :($gradient_symbol = $total))
        for output in 1:noutputs
            slot = plan.scatter_slots[(output - 1) * ndeps + d]
            coefficient = plan.output_coefficients[output]
            push!(statements, :(core[$slot] += $coefficient * $gradient_symbol))
        end
    end

    return Expr(:block, statements...)
end

function structured_core_body(jac)
    r = jac.reactor.gas.reaction
    forward_plans, reverse_plans = jac.float_rhs.forward_plans, jac.float_rhs.reverse_plans
    blocks = [reaction_block(i, jac.plans[i], forward_plans[i], reverse_plans[i], r.is_reversible[i])
              for i in 1:r.n_reactions]
    return quote
        workspace = jac.float_rhs.workspace
        r = jac.reactor.gas.reaction
        C = workspace.C
        kinetics = workspace.kinetics
        Arrhenius._rate_factors!(r, T, C, workspace.entropy, workspace.h_mole, kinetics;
            rate_multipliers=jac.reactor.rate_multipliers, pressure=P)
        core = jac.concentration_core.nzval
        fill!(core, 0.0)
        $(blocks...)
        return nothing
    end
end

function structured_kernel_id(jac)
    r = jac.reactor.gas.reaction
    return UInt64(hash((
        r.n_reactions,
        [(p.kind, p.indices, p.orders, p.repeated) for p in jac.float_rhs.forward_plans],
        [(p.kind, p.indices, p.orders, p.repeated) for p in jac.float_rhs.reverse_plans],
        [(p.dependencies, p.forward_to_dependency, p.reverse_to_dependency,
          p.efficiency_to_dependency, p.efficiency_deviations, p.scatter_slots,
          p.output_coefficients, p.collider_kind) for p in jac.plans],
    ), UInt(0)))
end

function GeneratedStructuredKernel(jac)
    id = structured_kernel_id(jac)
    if !haskey(GENERATED_STRUCTURED_CORE_BODIES, id)
        GENERATED_STRUCTURED_CORE_BODIES[id] = structured_core_body(jac)
    end
    return GeneratedStructuredKernel{id}()
end

@generated function _generated_fill_concentration_core!(jac, ::GeneratedStructuredKernel{ID}, T, P) where {ID}
    return GENERATED_STRUCTURED_CORE_BODIES[ID]
end

"""
    generated_structured_jacobian(reactor)

Like `structured_trial_jacobian(reactor)` (same scope, same guard, same
sparsity, same KLU/QNDF wiring downstream), except `concentration_core` is
filled by a mechanism-specialized compiled kernel instead of a runtime loop
over the reaction plans. Drop-in anywhere `structured_trial_jacobian`'s
result is used, including `guarded_structured_jacobian`/
`GuardedStructuredJacobian` and `qndf_structured_adapter`.
"""
function generated_structured_jacobian(reactor)
    jac = structured_trial_jacobian(reactor)
    kernel = GeneratedStructuredKernel(jac)
    jac.fill_core! = (j, T, P) -> _generated_fill_concentration_core!(j, kernel, T, P)
    return jac
end

export generated_structured_jacobian
