"""
    ReorderedReactions(reactions::Reactions)

Build a [`ReorderedReactions`](@ref) from `reactions`, bucketing reactions
into contiguous ranges by kind -- plain, then three-body, then falloff --
instead of `Reactions`' original mechanism-file order. `_rate_factors!`,
`_mass_action!` and `wdot!` are untyped/duck-typed over their `reaction`
argument, so every existing evaluator keeps working unchanged against a
`ReorderedReactions` value: it exposes the exact same field names as
`Reactions`, just column/row-permuted so that `index_three_body` and
`index_falloff` become plain contiguous ranges (`index_three_body_start:
index_falloff_start-1` and `index_falloff_start:n_reactions`) rather than
scattered index lists. That removes the gather/scatter indirection identified
as the scaling bottleneck in `validation/kinetics_scalar_scaling.md`, while
every value stays bit-identical to the unreordered evaluation.

`wdot!`'s net production-rate output does not depend on reaction order (it
sums each reaction's contribution independently), so calling it with a
`ReorderedReactions` and no `rate_multipliers`/`get_qdot`/`get_rate_constants`
gives exactly the same `wdot` as calling it with the original `Reactions` --
reordering is purely an internal performance detail there. Anything that
reports *per-reaction* values in mechanism order (`rate_multipliers`,
`reaction_rate_constants`, `wdot!(...;get_qdot=true)`, future per-reaction
sensitivities) must translate through `original_index`/`reordered_index`,
via the [`original_reaction_index`](@ref)/[`reordered_reaction_index`](@ref)
accessors below.
"""
struct ReorderedReactions{T<:AbstractFloat} <: AbstractReactions{T}
    product_stoich_coeffs::SparseMatrixCSC{T,Int64}
    reactant_stoich_coeffs::SparseMatrixCSC{T,Int64}
    reactant_orders::SparseMatrixCSC{T,Int64}
    is_reversible::Array{Bool,1}
    Arrhenius_coeffs::Array{T,2}
    Arrhenius_0::Array{T,2}
    Troe_::Array{T,2}
    index_three_body::Array{Int64,1}
    index_falloff::Array{Int64,1}
    index_falloff_Troe::Array{Int64,1}
    efficiencies_coeffs::SparseMatrixCSC{T,Int64}
    i_reactant::Array{Array{Int64,1},1}
    i_product::Array{Array{Int64,1},1}
    n_reactions::Int64
    vk::SparseMatrixCSC{T,Int64}
    vk_sum::Array{T,1}
    plog::PlogData{T}
    reactant_mass_action_terms::Vector{NTuple{4,Int64}}
    product_mass_action_terms::Vector{NTuple{4,Int64}}
    blowers_masel::BlowersMaselData{T}
    index_three_body_start::Int64
    index_falloff_start::Int64
    # original_index[new_i] = old_i (this reaction's position in `reactions`'
    # own mechanism-file order); reordered_index is its exact inverse.
    original_index::Vector{Int64}
    reordered_index::Vector{Int64}
end

function ReorderedReactions(reactions::Reactions{T}) where {T}
    n = reactions.n_reactions
    is_plain = trues(n)
    @inbounds for i in reactions.index_three_body
        is_plain[i] = false
    end
    @inbounds for i in reactions.index_falloff
        is_plain[i] = false
    end
    plain = findall(is_plain)

    # A reaction's Cantera `type:` selects exactly one of plain/three-body/
    # falloff at parse time (GasMechanismKineticsIO.jl's type dispatch), so
    # these three groups are disjoint and together cover every reaction.
    perm = vcat(plain, reactions.index_three_body, reactions.index_falloff)
    length(perm) == n && allunique(perm) ||
        throw(ArgumentError("index_three_body/index_falloff overlap or are out of range; cannot reorder"))

    original_index = perm
    reordered_index = Vector{Int64}(undef, n)
    @inbounds for (new_i, old_i) in enumerate(perm)
        reordered_index[old_i] = new_i
    end

    index_three_body_start = length(plain) + 1
    index_falloff_start = index_three_body_start + length(reactions.index_three_body)

    remap(indices) = [reordered_index[i] for i in indices]

    plog = PlogData(
        remap(reactions.plog.reaction_indices),
        reactions.plog.collider_indices,
        reactions.plog.group_offsets,
        reactions.plog.pressures,
        reactions.plog.rate_offsets,
        reactions.plog.Arrhenius_coeffs,
        reactions.plog.lmr_offsets,
        reactions.plog.lmr_species,
        reactions.plog.lmr_eps,
    )
    blowers_masel = BlowersMaselData(
        remap(reactions.blowers_masel.reaction_indices),
        reactions.blowers_masel.coefficients,
    )

    return ReorderedReactions{T}(
        reactions.product_stoich_coeffs[:, perm],
        reactions.reactant_stoich_coeffs[:, perm],
        reactions.reactant_orders[:, perm],
        reactions.is_reversible[perm],
        reactions.Arrhenius_coeffs[perm, :],
        reactions.Arrhenius_0,
        reactions.Troe_,
        collect(index_three_body_start:(index_falloff_start - 1)),
        collect(index_falloff_start:(index_falloff_start + length(reactions.index_falloff) - 1)),
        reactions.index_falloff_Troe,
        reactions.efficiencies_coeffs[:, perm],
        reactions.i_reactant[perm],
        reactions.i_product[perm],
        n,
        reactions.vk[:, perm],
        reactions.vk_sum[perm],
        plog,
        reactions.reactant_mass_action_terms[perm],
        reactions.product_mass_action_terms[perm],
        blowers_masel,
        index_three_body_start,
        index_falloff_start,
        original_index,
        reordered_index,
    )
end

"""
    original_reaction_index(reactions, i)

The index reaction `i` (in `reactions`' own storage order) occupies in the
mechanism's original YAML/Cantera reaction list. Identity for `Reactions`;
for `ReorderedReactions`, translates a bucketed position back to mechanism
order.
"""
original_reaction_index(::Reactions, i::Integer) = i
original_reaction_index(reactions::ReorderedReactions, i::Integer) = reactions.original_index[i]

"""
    reordered_reaction_index(reactions, original_i)

The storage position of the mechanism's `original_i`-th reaction within
`reactions`. Identity for `Reactions`; for `ReorderedReactions`, the inverse
of `original_reaction_index`.
"""
reordered_reaction_index(::Reactions, original_i::Integer) = original_i
reordered_reaction_index(reactions::ReorderedReactions, original_i::Integer) = reactions.reordered_index[original_i]

export ReorderedReactions, original_reaction_index, reordered_reaction_index
