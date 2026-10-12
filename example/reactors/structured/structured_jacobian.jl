using Arrhenius
using LinearAlgebra
using SparseArrays
import ForwardDiff

# The signed RHS owns the SignedProductPlan type extended below.
if !isdefined(@__MODULE__, :SignedTrialRHS)
    include(joinpath(@__DIR__, "signed_rhs.jl"))
end
if !isdefined(@__MODULE__, :signed_product_gradient!)
    include(joinpath(@__DIR__, "signed_product_gradient.jl"))
end

struct StructuredJacobianTag end
const StructuredJacobianDual = ForwardDiff.Dual{StructuredJacobianTag,Float64,1}

struct StructuredReactionPlan
    output_rows::Vector{Int}
    output_coefficients::Vector{Float64}
    dependencies::Vector{Int}
    forward_to_dependency::Vector{Int}
    reverse_to_dependency::Vector{Int}
    efficiency_to_dependency::Vector{Int}
    efficiency_deviations::Vector{Float64}
    scatter_slots::Vector{Int}
    falloff_position::Int
    collider_kind::UInt8 # 0 none, 1 three body, 2 falloff
end

struct StructuredMechanismSnapshot{T}
    data::T
end

mutable struct StructuredTrialJacobian{R,RF,RD,P,S}
    reactor::R
    float_rhs::RF
    temperature_rhs::RD
    plans::P
    structure_snapshot::S
    concentration_core::SparseMatrixCSC{Float64,Int}
    mass_core::SparseMatrixCSC{Float64,Int}
    pattern_colptr::Vector{Int}
    pattern_rowval::Vector{Int}
    q::Vector{Float64}
    rank_vector::Vector{Float64}
    temperature_column::Vector{Float64}
    energy_row::Vector{Float64}
    concentration_product::Vector{Float64}
    rate_gradient::Vector{Float64}
    forward_gradients::Vector{Vector{Float64}}
    reverse_gradients::Vector{Vector{Float64}}
    float_derivative::Vector{Float64}
    dual_state::Vector{StructuredJacobianDual}
    dual_derivative::Vector{StructuredJacobianDual}
    input_state::Vector{Float64}
end

function _modal_efficiency(efficiencies, reaction_index)
    checkbounds(efficiencies, :, reaction_index)
    values_in_column = @view efficiencies[:, reaction_index]
    all(isfinite, values_in_column) || throw(DomainError(
        reaction_index, "collider efficiencies must be finite"))
    all(>=(zero(eltype(efficiencies))), values_in_column) || throw(ArgumentError(
        "collider efficiencies must be nonnegative"))
    counts = Dict{Float64,Int}()
    @inbounds for species in axes(efficiencies, 1)
        value = Float64(efficiencies[species, reaction_index])
        counts[value] = get(counts, value, 0) + 1
    end
    largest = maximum(values(counts))
    modes = [value for (value, count) in counts if count == largest]
    length(modes) == 1 || throw(UnsupportedStructuredReactor(
        "collider efficiency column $reaction_index has no unique mode"))
    mode = only(modes)
    mode in (0.0, 1.0) || throw(UnsupportedStructuredReactor(
        "unsupported collider efficiency mode $mode in reaction $reaction_index"))
    return mode
end

function _stored_slot(A::SparseMatrixCSC, row::Int, column::Int)
    first = A.colptr[column]
    last = A.colptr[column + 1] - 1
    @inbounds for slot in first:last
        A.rowval[slot] == row && return slot
        A.rowval[slot] > row && break
    end
    error("internal sparse pattern omitted ($row, $column)")
end

function _mechanism_snapshot(reactor)
    r = reactor.gas.reaction
    multipliers = isnothing(reactor.rate_multipliers) ? nothing : copy(reactor.rate_multipliers)
    StructuredMechanismSnapshot((
        MW=copy(reactor.gas.MW),
        thermo_low=copy(reactor.gas.thermo.nasa_low),
        thermo_high=copy(reactor.gas.thermo.nasa_high),
        thermo_range=copy(reactor.gas.thermo.Trange),
        thermo_common=reactor.gas.thermo.isTcommon,
        product=copy(r.product_stoich_coeffs),
        reactant=copy(r.reactant_stoich_coeffs), orders=copy(r.reactant_orders),
        reversible=copy(r.is_reversible), arrhenius=copy(r.Arrhenius_coeffs),
        low=copy(r.Arrhenius_0), troe=copy(r.Troe_),
        three_body=copy(r.index_three_body), falloff=copy(r.index_falloff),
        falloff_troe=copy(r.index_falloff_Troe), efficiencies=copy(r.efficiencies_coeffs),
        vk=copy(r.vk), vk_sum=copy(r.vk_sum),
        plog_reactions=copy(r.plog.reaction_indices),
        plog_colliders=copy(r.plog.collider_indices),
        plog_groups=copy(r.plog.group_offsets), plog_pressures=copy(r.plog.pressures),
        plog_rates=copy(r.plog.rate_offsets), plog_arrhenius=copy(r.plog.Arrhenius_coeffs),
        bm_reactions=copy(r.blowers_masel.reaction_indices),
        bm_coefficients=copy(r.blowers_masel.coefficients),
        multipliers=multipliers,
    ))
end

function _check_mechanism!(jac::StructuredTrialJacobian)
    reactor = jac.reactor
    r = reactor.gas.reaction
    current = (
        MW=reactor.gas.MW, thermo_low=reactor.gas.thermo.nasa_low,
        thermo_high=reactor.gas.thermo.nasa_high,
        thermo_range=reactor.gas.thermo.Trange,
        thermo_common=reactor.gas.thermo.isTcommon, product=r.product_stoich_coeffs,
        reactant=r.reactant_stoich_coeffs, orders=r.reactant_orders,
        reversible=r.is_reversible, arrhenius=r.Arrhenius_coeffs,
        low=r.Arrhenius_0, troe=r.Troe_, three_body=r.index_three_body,
        falloff=r.index_falloff, falloff_troe=r.index_falloff_Troe,
        efficiencies=r.efficiencies_coeffs, vk=r.vk, vk_sum=r.vk_sum,
        plog_reactions=r.plog.reaction_indices, plog_colliders=r.plog.collider_indices,
        plog_groups=r.plog.group_offsets, plog_pressures=r.plog.pressures,
        plog_rates=r.plog.rate_offsets, plog_arrhenius=r.plog.Arrhenius_coeffs,
        bm_reactions=r.blowers_masel.reaction_indices,
        bm_coefficients=r.blowers_masel.coefficients,
        multipliers=reactor.rate_multipliers,
    )
    current == jac.structure_snapshot.data || throw(ArgumentError(
        "mechanism or rate multipliers changed; rebuild the structured Jacobian"))
    jac.concentration_core.colptr == jac.pattern_colptr &&
        jac.concentration_core.rowval == jac.pattern_rowval &&
        jac.mass_core.colptr == jac.pattern_colptr &&
        jac.mass_core.rowval == jac.pattern_rowval ||
        error("structured Jacobian sparse pattern changed")
    nothing
end

function _validate_structured_scope(reactor)
    reactor.constraint === :constant_pressure || throw(UnsupportedStructuredReactor(
        "structured Jacobian requires a constant-pressure reactor"))
    reactor.energy === :adiabatic || throw(UnsupportedStructuredReactor(
        "structured Jacobian requires an adiabatic reactor"))
    gas = reactor.gas
    eltype(gas.MW) === Float64 || throw(UnsupportedStructuredReactor(
        "structured Jacobian currently requires Float64 mechanism data"))
    isnothing(gas.thermo.extra) || throw(UnsupportedStructuredReactor(
        "non-NASA7 thermo is outside this structured Jacobian scope"))
    r = gas.reaction
    ns, nr = gas.n_species, gas.n_reactions
    expected = (ns, nr)
    size(r.reactant_stoich_coeffs) == expected || throw(DimensionMismatch(
        "reactant stoichiometry dimensions differ from the mechanism"))
    size(r.product_stoich_coeffs) == expected || throw(DimensionMismatch(
        "product stoichiometry dimensions differ from the mechanism"))
    size(r.reactant_orders) == expected || throw(DimensionMismatch(
        "reaction-order dimensions differ from the mechanism"))
    size(r.efficiencies_coeffs) == expected || throw(DimensionMismatch(
        "collider-efficiency dimensions differ from the mechanism"))
    size(r.vk) == expected || throw(DimensionMismatch(
        "net-stoichiometry dimensions differ from the mechanism"))
    length(r.is_reversible) == nr && length(r.vk_sum) == nr &&
        size(r.Arrhenius_coeffs, 1) == nr || throw(DimensionMismatch(
        "reaction metadata dimensions differ from the mechanism"))
    length(r.index_falloff) == length(r.index_falloff_Troe) || throw(DimensionMismatch(
        "falloff/Troe index dimensions differ"))
    length(r.plog.reaction_indices) == length(r.plog.collider_indices) ||
        throw(DimensionMismatch("PLOG reaction/collider index dimensions differ"))
    length(r.blowers_masel.reaction_indices) == size(r.blowers_masel.coefficients, 1) ||
        throw(DimensionMismatch("Blowers-Masel reaction/coefficient dimensions differ"))

    reactant_values = nonzeros(r.reactant_stoich_coeffs)
    product_values = nonzeros(r.product_stoich_coeffs)
    order_values = nonzeros(r.reactant_orders)
    efficiency_values = nonzeros(r.efficiencies_coeffs)
    all(isfinite, reactant_values) && all(isfinite, product_values) ||
        throw(DomainError(nothing, "stoichiometry must be finite"))
    all(>=(zero(eltype(r.reactant_stoich_coeffs))), reactant_values) &&
        all(>=(zero(eltype(r.product_stoich_coeffs))), product_values) ||
        throw(ArgumentError("stoichiometry must be nonnegative"))
    all(isfinite, order_values) || throw(DomainError(nothing, "reaction orders must be finite"))
    all(isfinite, efficiency_values) || throw(DomainError(nothing, "collider efficiencies must be finite"))
    all(>=(zero(eltype(r.efficiencies_coeffs))), efficiency_values) ||
        throw(ArgumentError("collider efficiencies must be nonnegative"))

    families = vcat(r.index_three_body, r.index_falloff,
                    r.plog.reaction_indices, r.blowers_masel.reaction_indices)
    all(i -> 1 <= i <= nr, families) || throw(BoundsError(1:nr, families))
    all(i -> 0 <= i <= ns, r.plog.collider_indices) ||
        throw(BoundsError(0:ns, r.plog.collider_indices))
    all(i -> i == -1 || 1 <= i <= size(r.Troe_, 1), r.index_falloff_Troe) ||
        throw(BoundsError(axes(r.Troe_, 1), r.index_falloff_Troe))

    r.reactant_orders == r.reactant_stoich_coeffs || throw(UnsupportedStructuredReactor(
        "custom reaction orders are outside this structured Jacobian scope"))
    all(isinteger, reactant_values) && all(isinteger, product_values) ||
        throw(UnsupportedStructuredReactor("integer stoichiometry is required"))
    isempty(r.blowers_masel.reaction_indices) || throw(UnsupportedStructuredReactor(
        "Blowers-Masel reactions are outside this scope"))
    all(iszero, r.plog.collider_indices) || throw(UnsupportedStructuredReactor(
        "PLOG reactions with concentration colliders are unsupported"))
    length(unique(families)) == length(families) || throw(UnsupportedStructuredReactor(
        "overlapping special-rate reaction families are unsupported"))
    nothing
end

function structured_trial_jacobian(reactor)
    _validate_structured_scope(reactor)
    gas = reactor.gas
    r = gas.reaction
    ns, nr = gas.n_species, gas.n_reactions
    float_rhs = signed_trial_rhs(reactor; scalar_type=Float64)
    temperature_rhs = signed_trial_rhs(reactor; scalar_type=StructuredJacobianDual)
    forward_plans, reverse_plans = float_rhs.forward_plans, float_rhs.reverse_plans
    three_body = Set(r.index_three_body)
    falloff_positions = Dict(index => position for (position, index) in enumerate(r.index_falloff))
    pattern_rows, pattern_columns = collect(1:ns), collect(1:ns)
    partial = Vector{NamedTuple}(undef, nr)
    vk_rows, vk_values = rowvals(r.vk), nonzeros(r.vk)

    for i in 1:nr
        output_slots = collect(nzrange(r.vk, i))
        outputs = collect(vk_rows[output_slots])
        coefficients = Float64.(vk_values[output_slots])
        dependencies = copy(forward_plans[i].indices)
        r.is_reversible[i] && append!(dependencies, reverse_plans[i].indices)
        collider_kind, falloff_position = UInt8(0), 0
        efficiency_indices, deviations = Int[], Float64[]
        if i in three_body || haskey(falloff_positions, i)
            collider_kind = i in three_body ? UInt8(1) : UInt8(2)
            falloff_position = get(falloff_positions, i, 0)
            mode = _modal_efficiency(r.efficiencies_coeffs, i)
            @inbounds for species in 1:ns
                deviation = r.efficiencies_coeffs[species, i] - mode
                if !iszero(deviation)
                    push!(efficiency_indices, species)
                    push!(deviations, deviation)
                end
            end
            append!(dependencies, efficiency_indices)
        end
        sort!(unique!(dependencies))
        append!(pattern_rows, repeat(outputs, inner=length(dependencies)))
        append!(pattern_columns, repeat(dependencies, outer=length(outputs)))
        partial[i] = (outputs=outputs, coefficients=coefficients,
            dependencies=dependencies,
            forward_map=[searchsortedfirst(dependencies, k) for k in forward_plans[i].indices],
            reverse_map=[searchsortedfirst(dependencies, k) for k in reverse_plans[i].indices],
            efficiency_map=[searchsortedfirst(dependencies, k) for k in efficiency_indices],
            deviations=deviations, falloff_position=falloff_position,
            collider_kind=collider_kind)
    end

    concentration_core = sparse(pattern_rows, pattern_columns,
                                ones(Float64, length(pattern_rows)), ns, ns)
    fill!(nonzeros(concentration_core), 0.0)
    mass_core = copy(concentration_core)
    plans = Vector{StructuredReactionPlan}(undef, nr)
    for i in 1:nr
        x = partial[i]
        scatter = Int[]
        for row in x.outputs, column in x.dependencies
            push!(scatter, _stored_slot(concentration_core, row, column))
        end
        plans[i] = StructuredReactionPlan(x.outputs, x.coefficients, x.dependencies,
            x.forward_map, x.reverse_map, x.efficiency_map, x.deviations,
            scatter, x.falloff_position, x.collider_kind)
    end
    maximum_dependencies = maximum(length(plan.dependencies) for plan in plans; init=0)
    n = ns + 1
    dual_zero = StructuredJacobianDual(0.0, ForwardDiff.Partials((0.0,)))
    StructuredTrialJacobian(reactor, float_rhs, temperature_rhs, plans,
        _mechanism_snapshot(reactor), concentration_core, mass_core,
        copy(concentration_core.colptr), copy(concentration_core.rowval),
        zeros(ns), zeros(ns), zeros(n), zeros(ns), zeros(ns),
        zeros(maximum_dependencies),
        [zeros(length(plan.indices)) for plan in forward_plans],
        [zeros(length(plan.indices)) for plan in reverse_plans],
        zeros(n), fill(dual_zero, n), fill(dual_zero, n), zeros(n))
end

function _ambiguous_anyn_zero(plan::SignedProductPlan, C)
    plan.kind == 4 || return false
    zero_position = 0
    @inbounds for i in eachindex(plan.indices)
        plan.orders[i] == 0 && continue
        value = C[plan.indices[i]]
        value < 0 && return false
        if iszero(value)
            zero_position == 0 || return false
            zero_position = i
        end
    end
    zero_position != 0 && plan.orders[zero_position] == 1.0
end

function _collider_derivative(jac, i, plan, T, C, h, kinetics)
    plan.collider_kind == 0 && return 0.0
    r = jac.reactor.gas.reaction
    multiplier = isnothing(jac.reactor.rate_multipliers) ? 1.0 : jac.reactor.rate_multipliers[i]
    high = Arrhenius._flame_base_rate(r, i, T, h)
    factor_derivative = if plan.collider_kind == 1
        high
    else
        collider = dot(@view(r.efficiencies_coeffs[:, i]), C)
        Arrhenius._flame_falloff_derivative(
            r, plan.falloff_position, i, T, collider, high)
    end
    Pf = signed_multiply(jac.float_rhs.forward_plans[i], C, 1.0)
    Pr = r.is_reversible[i] ? signed_multiply(jac.float_rhs.reverse_plans[i], C, 1.0) : 0.0
    balance = Pf - (r.is_reversible[i] ?
        Pr / kinetics.equilibrium_constants[i] : 0.0)
    factor_derivative * multiplier * balance
end

function _fill_concentration_core!(jac, T, P)
    reactor, workspace = jac.reactor, jac.float_rhs.workspace
    r, C, kinetics = reactor.gas.reaction, workspace.C, workspace.kinetics
    Arrhenius._rate_factors!(r, T, C, workspace.entropy, workspace.h_mole, kinetics;
        rate_multipliers=reactor.rate_multipliers, pressure=P)
    fill!(nonzeros(jac.concentration_core), 0.0)
    @inbounds for i in 1:r.n_reactions
        fp, rp, plan = jac.float_rhs.forward_plans[i], jac.float_rhs.reverse_plans[i], jac.plans[i]
        (_ambiguous_anyn_zero(fp, C) ||
         (r.is_reversible[i] && _ambiguous_anyn_zero(rp, C))) && throw(DomainError(
            i, "ambiguous AnyN order-one derivative at a zero concentration"))
        ndeps = length(plan.dependencies)
        for k in 1:ndeps
            jac.rate_gradient[k] = 0.0
        end
        fg = jac.forward_gradients[i]
        signed_product_gradient!(fg, fp, C, kinetics.kf[i])
        for k in eachindex(fg)
            jac.rate_gradient[plan.forward_to_dependency[k]] += fg[k]
        end
        if r.is_reversible[i]
            rg = jac.reverse_gradients[i]
            signed_product_gradient!(rg, rp, C, -kinetics.kr[i])
            for k in eachindex(rg)
                jac.rate_gradient[plan.reverse_to_dependency[k]] += rg[k]
            end
        end
        collider = _collider_derivative(jac, i, plan, T, C, workspace.h_mole, kinetics)
        if !iszero(collider)
            for k in eachindex(plan.efficiency_deviations)
                jac.rate_gradient[plan.efficiency_to_dependency[k]] +=
                    collider * plan.efficiency_deviations[k]
            end
        end
        for output in eachindex(plan.output_rows)
            coefficient = plan.output_coefficients[output]
            offset = (output - 1) * ndeps
            for dependency in 1:ndeps
                jac.concentration_core.nzval[plan.scatter_slots[offset + dependency]] +=
                    coefficient * jac.rate_gradient[dependency]
            end
        end
    end
    nothing
end

function _temperature_direction!(jac, state, t)
    n = length(state)
    @inbounds for i in 1:n
        jac.dual_state[i] = StructuredJacobianDual(
            state[i], ForwardDiff.Partials((i == n ? 1.0 : 0.0,)))
    end
    jac.temperature_rhs(jac.dual_derivative, jac.dual_state, nothing, t)
    @inbounds for i in 1:n
        jac.temperature_column[i] = ForwardDiff.partials(jac.dual_derivative[i])[1]
    end
    nothing
end

function (jac::StructuredTrialJacobian)(J, state, p, t)
    ns, n = jac.reactor.gas.n_species, jac.reactor.gas.n_species + 1
    length(state) == n || throw(DimensionMismatch("state must contain n_species + 1 entries"))
    size(J) == (n, n) || throw(DimensionMismatch("Jacobian must have state dimensions"))
    eltype(state) === Float64 || throw(ArgumentError("structured callback requires Float64 state"))
    eltype(J) === Float64 || throw(ArgumentError("structured callback requires Float64 output"))
    _check_mechanism!(jac)
    copyto!(jac.input_state, state)
    T, P, density, inverse_mw = Arrhenius._reactor_tpρ(jac.reactor, state)
    jac.float_rhs(jac.float_derivative, state, nothing, t)
    _fill_concentration_core!(jac, T, P)
    workspace, MW = jac.float_rhs.workspace, jac.reactor.gas.MW
    core, mass_core = jac.concentration_core, jac.mass_core
    @inbounds for column in 1:ns
        for slot in nzrange(core, column)
            row = core.rowval[slot]
            mass_core.nzval[slot] = MW[row] * core.nzval[slot] / MW[column]
        end
        jac.q[column] = 1.0 / (inverse_mw * MW[column])
    end
    mul!(jac.concentration_product, core, workspace.C)
    @inbounds for row in 1:ns
        jac.rank_vector[row] = jac.float_derivative[row] -
            MW[row] * jac.concentration_product[row] / density
    end
    fill!(J, 0.0)
    @inbounds for column in 1:ns, row in 1:ns
        J[row, column] = jac.rank_vector[row] * jac.q[column]
    end
    @inbounds for column in 1:ns
        for slot in nzrange(mass_core, column)
            J[mass_core.rowval[slot], column] += mass_core.nzval[slot]
        end
    end
    _temperature_direction!(jac, state, t)
    @inbounds for row in 1:ns
        J[row, n] = jac.temperature_column[row]
    end
    capacity = 0.0
    @inbounds for species in 1:ns
        capacity += state[species] * workspace.cp_R[species] / MW[species]
    end
    capacity *= Arrhenius.R
    isfinite(capacity) && capacity > 0 || throw(DomainError(capacity, "positive heat capacity required"))
    fT = jac.float_derivative[n]
    @inbounds for column in 1:ns
        enthalpy_term = 0.0
        for species in 1:ns
            enthalpy_term += workspace.h_mole[species] / MW[species] * J[species, column]
        end
        cp_gradient = Arrhenius.R * workspace.cp_R[column] / MW[column]
        jac.energy_row[column] = -enthalpy_term / capacity - fT * cp_gradient / capacity
        J[n, column] = jac.energy_row[column]
    end
    J[n, n] = jac.temperature_column[n]
    state == jac.input_state || error("structured Jacobian mutated its input state")
    all(isfinite, J) || throw(DomainError(J, "structured Jacobian is not finite"))
    nothing
end