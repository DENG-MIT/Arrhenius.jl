using YAML, NPZ

struct Reaction
    product_stoich_coeffs
    reactant_stoich_coeffs
    Arrhenius_coeffs
    Arrhenius_A0
    Arrhenius_b0
    Arrhenius_Ea0
    Troe_A
    Troe_T1
    Troe_T2
    Troe_T3
    index_three_body
    index_falloff
    efficiencies_coeffs
end

struct Thermo
    nasa_low
    nasa_high
end

struct Solution
    n_species::Int
    n_reactions::Int
    molecular_weights
    species_names
    thermo::Thermo
    reaction::Reaction
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
    molecular_weights = npz["molecular_weights"]
    efficiencies_coeffs = npz["efficiencies_coeffs"]
    product_stoich_coeffs = npz["product_stoich_coeffs"]
    reactant_stoich_coeffs = npz["reactant_stoich_coeffs"]
    Arrhenius_coeffs = npz["Arrhenius_coeffs"]
    Arrhenius_A0 = npz["Arrhenius_A0"]
    Arrhenius_b0 = npz["Arrhenius_b0"]
    Arrhenius_Ea0 = npz["Arrhenius_Ea0"]
    Troe_A = npz["Troe_A"]
    Troe_T1 = npz["Troe_T1"]
    Troe_T2 = npz["Troe_T2"]
    Troe_T3 = npz["Troe_T3"]

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

    reaction = Reaction(product_stoich_coeffs, reactant_stoich_coeffs,
                        Arrhenius_coeffs, Arrhenius_A0, Arrhenius_b0, Arrhenius_Ea0,
                        Troe_A, Troe_T1, Troe_T2, Troe_T3, index_three_body, index_falloff, efficiencies_coeffs)

    gas = Solution(n_species, n_reactions, molecular_weights, species_names, thermo, reaction)
end
