function read_species_basics(yaml)
    n_species = length(yaml["phases"][1]["species"])
    n_reactions = length(yaml["reactions"])
    species_names = yaml["phases"][1]["species"]
    elements = yaml["phases"][1]["elements"]
    n_elements = length(elements)

    ele_matrix = zeros(n_elements, n_species)

    _species_names =
        [yaml["species"][i]["name"] for i = 1:length(yaml["species"])]

    for (i, species) in enumerate(species_names)
        spec = yaml["species"][findfirst(x -> x == species, _species_names)]
        for j = 1:n_elements
            if haskey(spec["composition"], elements[j])
                ele_matrix[j, i] = spec["composition"][elements[j]]
            end
        end
    end

    return n_species, n_reactions, species_names, elements, n_elements, ele_matrix
end
"""
    CreateSolution(mech)
    
Reaction mechanism is interepreted here. Part of the infomation are read in
from the yaml file, pary of them are from the pre-processed .npz file from
ReacTorch and Cantera
test for math enviroment
"""
function CreateSolution(mech)
    yaml = YAML.load_file(mech)

    #### Basics
    n_species, n_reactions, species_names,
    elements, n_elements, ele_matrix = read_species_basics(yaml)

    #### Thermo
    if yaml["phases"][1]["thermo"] == "ideal-gas" # switch to work with Cantera standard
        thermo = IdealGasThermo(yaml)
    else
        constructorThermo = Symbol(yaml["phases"][1]["thermo"], :Thermo)
        thermo = @eval($constructorThermo)(yaml)
    end


    #### Kinetic data

    npz = npzread("$mech.npz")
    MW = npz["molecular_weights"]
    efficiencies_coeffs_full = npz["efficiencies_coeffs"]
    product_stoich_coeffs = sparse(npz["product_stoich_coeffs"])
    reactant_stoich_coeffs = sparse(npz["reactant_stoich_coeffs"])
    reactant_orders = sparse(npz["reactant_orders"])
    is_reversible = npz["is_reversible"]
    Arrhenius_coeffs = npz["Arrhenius_coeffs"]
    if haskey(npz, "Arrhenius_A0")
        Arrhenius_A0 = npz["Arrhenius_A0"]
        Arrhenius_b0 = npz["Arrhenius_b0"]
        Arrhenius_Ea0 = npz["Arrhenius_Ea0"]
    else
        Arrhenius_A0 = []
        Arrhenius_b0 = []
        Arrhenius_Ea0 = []
    end

    if haskey(npz, "Troe_A")
        Troe_A = npz["Troe_A"]
        Troe_T1 = npz["Troe_T1"]
        Troe_T2 = npz["Troe_T2"]
        Troe_T3 = npz["Troe_T3"]
    else
        Troe_A = []
        Troe_T1 = []
        Troe_T2 = []
        Troe_T3 = []
    end
    Arrhenius_0 = hcat(Arrhenius_A0, Arrhenius_b0, Arrhenius_Ea0)
    Troe_ = hcat(Troe_A, Troe_T1, Troe_T2, Troe_T3)

    index_three_body = []
    index_falloff = []
    index_falloff_Troe = []
    j = 1
    for i = 1:n_reactions
        reaction = yaml["reactions"][i]
        if haskey(reaction, "type")
            if reaction["type"] == "three-body"
                push!(index_three_body, i)
            end
            if reaction["type"] == "falloff"
                push!(index_falloff, i)
                if haskey(reaction, "Troe")
                    push!(index_falloff_Troe, j)
                    j = j + 1
                else
                    push!(index_falloff_Troe, -1)
                end
            end
        end
    end

    i_reactant = []
    i_product = []
    for i = 1:n_reactions
        push!(i_reactant, findall(reactant_orders[:, i] .> 0.01))
        push!(i_product, findall(product_stoich_coeffs[:, i] .> 0.01))
    end

    vk = sparse(product_stoich_coeffs - reactant_stoich_coeffs)
    vk_sum = sum(vk, dims=1)[1, :]

    for i in 1:n_reactions
        if !((i in index_three_body) | (i in index_falloff))
            efficiencies_coeffs_full[:, i] .= 0.0
        end
    end
    efficiencies_coeffs = sparse(efficiencies_coeffs_full)

    reaction = Reaction(
        product_stoich_coeffs,
        reactant_stoich_coeffs,
        reactant_orders,
        is_reversible,
        Arrhenius_coeffs,
        Arrhenius_0,
        Troe_,
        index_three_body,
        index_falloff,
        index_falloff_Troe,
        efficiencies_coeffs,
        i_reactant,
        i_product,
        n_reactions,
        vk,
        vk_sum,
    )


    #### Transport data

    if haskey(npz, "species_viscosities_poly")
        species_viscosities_poly = npz["species_viscosities_poly"]
        thermal_conductivity_poly = npz["thermal_conductivity_poly"]
        binary_diff_coeffs_poly = npz["binary_diff_coeffs_poly"]
        poly_order = size(species_viscosities_poly)[1]
    else
        species_viscosities_poly = zeros(2, 2)
        thermal_conductivity_poly = zeros(2, 2)
        binary_diff_coeffs_poly = zeros(2, 2)
        poly_order = 6
    end

    trans = Transport(poly_order,
                      species_viscosities_poly, 
                      thermal_conductivity_poly, 
                      binary_diff_coeffs_poly)

    gas = Solution(
        n_species,
        n_reactions,
        MW,
        species_names,
        elements,
        ele_matrix,
        thermo,
        trans,
        reaction,
    )
    return gas
end

export CreateSolution