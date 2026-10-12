"""Flattened pressure groups and Arrhenius expressions for PLOG reactions.

Linear-Burke (LMR-R) reactions whose additional colliders give only an
efficiency are stored as the PLOG rate of collider `M` evaluated at the
effective pressure ``P ε_mix``; `lmr_offsets` index each PLOG reaction's
colliders in `lmr_species` and their efficiency Arrhenius rows `lmr_eps`
(A, b, Ea in cal/mol, relative to `M`).
"""
struct PlogData{T<:AbstractFloat}
    reaction_indices::Vector{Int64}
    collider_indices::Vector{Int64}
    group_offsets::Vector{Int64}
    pressures::Vector{T}
    rate_offsets::Vector{Int64}
    Arrhenius_coeffs::Matrix{T}
    lmr_offsets::Vector{Int64}
    lmr_species::Vector{Int64}
    lmr_eps::Matrix{T}
end

function PlogData(reaction_indices, collider_indices, group_offsets, pressures, rate_offsets, coefficients)
    return PlogData(reaction_indices, collider_indices, group_offsets, pressures, rate_offsets,
        coefficients, ones(Int64, length(reaction_indices)+1), Int64[], zeros(eltype(pressures), 0, 3))
end

function PlogData(reaction_indices, group_offsets, pressures, rate_offsets, coefficients)
    return PlogData(
        reaction_indices,
        zeros(Int64, length(reaction_indices)),
        group_offsets,
        pressures,
        rate_offsets,
        coefficients,
    )
end

@inline _plog_linear_burke(plog::PlogData, p) = @inbounds plog.lmr_offsets[p+1] > plog.lmr_offsets[p]

struct Reaction{T<:AbstractFloat}
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
    blowers_masel::BlowersMaselData{T}
end

# Preserve the pre-Blowers–Masel constructor used by downstream mechanisms.
function Reaction(product,reactant,orders,reversible,arrhenius,low,troe,third,falloff,
        falloff_troe,efficiencies,reactant_indices,product_indices,n,vk,vk_sum,plog::PlogData{T}) where T
    return Reaction(product,reactant,orders,reversible,arrhenius,low,troe,third,falloff,
        falloff_troe,efficiencies,reactant_indices,product_indices,n,vk,vk_sum,plog,
        BlowersMaselData(Int64[],zeros(T,0,4)))
end

abstract type Thermo end

struct Transport{T<:AbstractFloat}
    poly_order::Int64
    species_viscosities_poly::Array{T,2}
    thermal_conductivity_poly::Array{T,2}
    binary_diff_coeffs_poly::Array{T,2}
    model::Symbol
end
Transport(order::Integer, viscosity::Matrix{T}, conductivity::Matrix{T}, binary::Matrix{T}) where {T<:AbstractFloat} =
    Transport(Int64(order),viscosity,conductivity,binary,:mixture_averaged)
Transport{T}(order::Integer, viscosity::Matrix{T}, conductivity::Matrix{T}, binary::Matrix{T}) where {T<:AbstractFloat} =
    Transport{T}(Int64(order),viscosity,conductivity,binary,:mixture_averaged)

struct Solution{T<:AbstractFloat,Th<:Thermo}
    n_species::Int64
    n_reactions::Int64
    MW::Vector{T}
    species_names::Vector{String}
    elements::Vector{String}
    ele_matrix::Array{T,2}
    thermo::Th
    trans::Transport{T}
    reaction::Reaction{T}
end
