"Reaction mechanism is interepreted here. Part of the infomation are read in
from the yaml file, pary of them are from the pre-processed .npz file from
ReacTorch and Cantera"
function CreateSolution(mech)
    yaml = YAML.load_file(mech)
    n_species = length(yaml["phases"][1]["species"])
    n_reactions = length(yaml["reactions"])
    species_names = yaml["phases"][1]["species"]
    elements = yaml["phases"][1]["elements"]
    n_elements = length(elements)

    ele_matrix = zeros(n_elements, n_species)

    nasa_low = zeros(n_species, 7)
    nasa_high = zeros(n_species, 7)

    _species_names = [yaml["species"][i]["name"] for i in 1:length(yaml["species"])]

    for (i, species) in enumerate(species_names)
        spec = yaml["species"][findfirst(x -> x == species, _species_names)]
        nasa_low[i, :] = spec["thermo"]["data"][1]
        nasa_high[i, :] = spec["thermo"]["data"][2]

        for j = 1:n_elements
            if haskey(spec["composition"], elements[j])
                ele_matrix[j, i] = spec["composition"][elements[j]]
            end
        end
    end

    thermo = Thermo(nasa_low, nasa_high)

    npz = npzread("$mech.npz")
    MW = npz["molecular_weights"]
    efficiencies_coeffs = npz["efficiencies_coeffs"]
    product_stoich_coeffs = npz["product_stoich_coeffs"]
    reactant_stoich_coeffs = npz["reactant_stoich_coeffs"]
    reactant_orders = npz["reactant_orders"]
    is_reversible = npz["is_reversible"]
    Arrhenius_coeffs = npz["Arrhenius_coeffs"]
    if haskey(npz, "Arrhenius_A0")
        Arrhenius_A0 = npz["Arrhenius_A0"]
        Arrhenius_b0 = npz["Arrhenius_b0"]
        Arrhenius_Ea0 = npz["Arrhenius_Ea0"]
        Troe_A = npz["Troe_A"]
        Troe_T1 = npz["Troe_T1"]
        Troe_T2 = npz["Troe_T2"]
        Troe_T3 = npz["Troe_T3"]
    else
        Arrhenius_A0 = []
        Arrhenius_b0 = []
        Arrhenius_Ea0 = []
        Troe_A = []
        Troe_T1 = []
        Troe_T2 = []
        Troe_T3 = []
    end

    index_three_body = []
    index_falloff = []
    index_falloff_troe = []
    for i = 1:n_reactions
        reaction = yaml["reactions"][i]
        if haskey(reaction, "type")
            if yaml["reactions"][i]["type"] == "three-body"
                push!(index_three_body, i)
            end
            if yaml["reactions"][i]["type"] == "falloff"
                push!(index_falloff, i)
                if haskey(yaml["reactions"][i], "Troe")
                    push!(index_falloff_troe, i)
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

    vk = product_stoich_coeffs - reactant_stoich_coeffs
    vk_sum = sum(vk, dims=1)[1, :]

    reaction = Reaction(
        product_stoich_coeffs,
        reactant_stoich_coeffs,
        reactant_orders,
        is_reversible,
        Arrhenius_coeffs,
        Arrhenius_A0,
        Arrhenius_b0,
        Arrhenius_Ea0,
        Troe_A,
        Troe_T1,
        Troe_T2,
        Troe_T3,
        index_three_body,
        index_falloff,
        index_falloff_troe,
        efficiencies_coeffs,
        i_reactant,
        i_product,
        n_reactions,
        vk,
        vk_sum,
    )

    gas = Solution(
        n_species,
        n_reactions,
        MW,
        species_names,
        elements,
        ele_matrix,
        thermo,
        reaction,
    )
    return gas
end

export CreateSolution
