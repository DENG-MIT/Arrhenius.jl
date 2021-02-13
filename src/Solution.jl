function set_states(gas::Solution, T, P, Y)
    mean_MW = 1.0 / dot(Y, 1 ./ gas.MW)
    ρ_mass = P / R / T * mean_MW
    X = Y2X(gas, Y, mean_MW)
    C = Y2C(gas, Y, ρ_mass)
    cp_mole = cp_mole_func(gas, T, X)
    cp_mass = cp_mass_func(gas, cp_mole, mean_MW)
    H_mole, h_mole = H_mole_func(gas, T, X)
    H_mass = H_mass_func(gas, h_mole, Y)
    S_mole, s_mole, S0 = S_mole_func(gas, T, P, X)
    S_mass = S_mass_func(gas, s_mole, Y)
    return wdot_func(gas.reaction, T, C, S0, h_mole)
end

function CreateSolution(mech)
    yaml = YAML.load_file(mech)
    n_species = length(yaml["species"])
    n_reactions = length(yaml["reactions"])
    species_names = yaml["phases"][1]["species"]

    nasa_low = zeros(n_species, 7)
    nasa_high = zeros(n_species, 7)

    for species in species_names
        j = findfirst(x -> x == species, yaml["phases"][1]["species"])
        nasa_low[j, :] = yaml["species"][j]["thermo"]["data"][1]
        nasa_high[j, :] = yaml["species"][j]["thermo"]["data"][2]
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
    for i = 1:n_reactions
        reaction = yaml["reactions"][i]
        if haskey(reaction, "type")
            if yaml["reactions"][i]["type"] == "three-body"
                push!(index_three_body, i)
            end
            if yaml["reactions"][i]["type"] == "falloff"
                push!(index_falloff, i)
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
    vk_sum = sum(vk, dims = 1)[1, :]

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
        efficiencies_coeffs,
        i_reactant,
        i_product,
        n_reactions,
        vk,
        vk_sum
    )

    gas = Solution(n_species, n_reactions, MW, species_names, thermo, reaction)
    return gas
end
