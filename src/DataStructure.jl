struct Reaction
    product_stoich_coeffs::Array{Float64, 2}
    reactant_stoich_coeffs::Array{Float64, 2}
    reactant_orders::Array{Float64, 2}
    is_reversible::Array{Bool, 1}
    Arrhenius_coeffs::Array{Float64, 2}
    Arrhenius_0::Array{Float64, 2}
    Troe_::Array{Float64, 2}
    index_three_body::Array{Int64, 1}
    index_falloff::Array{Int64, 1}
    index_falloff_troe::Array{Int64, 1}
    efficiencies_coeffs::Array{Float64, 2}
    i_reactant::Array{Array{Int64, 1},1}
    i_product::Array{Array{Int64, 1},1}
    n_reactions::Int64
    vk::Array{Float64, 2}
    vk_sum::Array{Float64, 1}
end

struct Thermo
    nasa_low::Array{Float64, 2}
    nasa_high::Array{Float64, 2}
    Trange::Array{Float64, 2}
    isTcommon::Bool
end

struct Solution
    n_species::Int64
    n_reactions::Int64
    MW::Vector{Float64}
    species_names::Vector{String}
    elements::Vector{String}
    ele_matrix::Array{Float64, 2}
    thermo::Thermo
    reaction::Reaction
end
