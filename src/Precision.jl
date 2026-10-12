function _precision_array(values, ::Type{T}, label) where {T<:AbstractFloat}
    converted = T.(values)
    overflow_index = findfirst(
        index -> isfinite(values[index]) && !isfinite(converted[index]),
        eachindex(values),
    )
    if !isnothing(overflow_index)
        throw(ArgumentError(
            "$label contains a value ($(values[overflow_index])) outside the " *
            "finite range of $T",
        ))
    end
    underflow_index = findfirst(
        index -> !iszero(values[index]) && iszero(converted[index]),
        eachindex(values),
    )
    if !isnothing(underflow_index)
        throw(ArgumentError(
            "$label contains a nonzero value ($(values[underflow_index])) " *
            "that underflows to zero in $T",
        ))
    end
    return converted
end

function _precision_sparse(matrix, ::Type{T}, label) where {T<:AbstractFloat}
    return SparseMatrixCSC{T,Int64}(
        size(matrix, 1),
        size(matrix, 2),
        copy(matrix.colptr),
        copy(rowvals(matrix)),
        _precision_array(nonzeros(matrix), T, label),
    )
end

function _convert_precision(plog::PlogData, ::Type{T}) where {T<:AbstractFloat}
    return PlogData(
        copy(plog.reaction_indices),
        copy(plog.collider_indices),
        copy(plog.group_offsets),
        _precision_array(plog.pressures, T, "PLOG pressures"),
        copy(plog.rate_offsets),
        _precision_array(plog.Arrhenius_coeffs, T, "PLOG Arrhenius coefficients"),
    )
end

function _convert_precision(reaction::Reaction, ::Type{T}) where {T<:AbstractFloat}
    return Reaction(
        _precision_sparse(reaction.product_stoich_coeffs, T, "product stoichiometry"),
        _precision_sparse(reaction.reactant_stoich_coeffs, T, "reactant stoichiometry"),
        _precision_sparse(reaction.reactant_orders, T, "reactant orders"),
        copy(reaction.is_reversible),
        _precision_array(reaction.Arrhenius_coeffs, T, "Arrhenius coefficients"),
        _precision_array(reaction.Arrhenius_0, T, "low-pressure Arrhenius coefficients"),
        _precision_array(reaction.Troe_, T, "Troe parameters"),
        copy(reaction.index_three_body),
        copy(reaction.index_falloff),
        copy(reaction.index_falloff_Troe),
        _precision_sparse(reaction.efficiencies_coeffs, T, "third-body efficiencies"),
        [copy(indices) for indices in reaction.i_reactant],
        [copy(indices) for indices in reaction.i_product],
        reaction.n_reactions,
        _precision_sparse(reaction.vk, T, "net stoichiometry"),
        _precision_array(reaction.vk_sum, T, "net stoichiometric sums"),
        _convert_precision(reaction.plog, T),
        BlowersMaselData(copy(reaction.blowers_masel.reaction_indices),
            _precision_array(reaction.blowers_masel.coefficients,T,"Blowers–Masel coefficients")),
    )
end

function _convert_precision(
    thermo::IdealGasThermo,
    ::Type{T},
) where {T<:AbstractFloat}
    extra = isnothing(thermo.extra) ? nothing : SpeciesThermoData(
        copy(thermo.extra.models),
        [_precision_array(values,T,"species temperature ranges") for values in thermo.extra.ranges],
        [_precision_array(values,T,"extended species thermo coefficients") for values in thermo.extra.coefficients],
    )
    return IdealGasThermo(
        _precision_array(thermo.nasa_low, T, "low-temperature NASA coefficients"),
        _precision_array(thermo.nasa_high, T, "high-temperature NASA coefficients"),
        _precision_array(thermo.Trange, T, "thermodynamic temperature ranges"),
        thermo.isTcommon,
        extra,
    )
end

function _convert_precision(transport::Transport, ::Type{T}) where {T<:AbstractFloat}
    return Transport(
        transport.poly_order,
        _precision_array(transport.species_viscosities_poly, T, "viscosity coefficients"),
        _precision_array(transport.thermal_conductivity_poly, T, "thermal-conductivity coefficients"),
        _precision_array(transport.binary_diff_coeffs_poly, T, "binary-diffusion coefficients"),
        transport.model,
    )
end

"""
    convert_precision(gas::Solution, T)

Return an independent copy of `gas` whose floating-point mechanism,
thermodynamic, transport, and molecular data use type `T`. Index arrays and
species names retain their original types. The default YAML loader continues
to construct `Float64` solutions.
"""
function convert_precision(gas::Solution, ::Type{T}) where {T<:AbstractFloat}
    return Solution(
        gas.n_species,
        gas.n_reactions,
        _precision_array(gas.MW, T, "molecular weights"),
        copy(gas.species_names),
        copy(gas.elements),
        _precision_array(gas.ele_matrix, T, "element matrix"),
        _convert_precision(gas.thermo, T),
        _convert_precision(gas.trans, T),
        _convert_precision(gas.reaction, T),
    )
end
export convert_precision
