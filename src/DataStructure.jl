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
    reaction::Reactions{T}
end
