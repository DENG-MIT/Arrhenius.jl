using YAML, NPZ
using LinearAlgebra

struct Reaction
    product_stoich_coeffs
    reactant_stoich_coeffs
    reactant_orders
    is_reversible
    list_type4_noTroe
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
    molecular_weights::Vector{Float64}
    species_names::Vector{String}
    thermo::Thermo
    reaction::Reaction
end

mutable struct MSolution
    Y::Vector{Float64}
    X::Vector{Float64}
    C::Vector{Float64}
    T::Float64
    P::Float64
    mean_molecular_weight::Float64
    ρ_mass::Float64
    cp_mole::Float64
    cp_mass::Float64
    H_mass::Float64
    H_mole::Float64
    h_mole::Vector{Float64}
    S_mass::Float64
    S_mole::Float64
    s_mole::Vector{Float64}
    S0::Vector{Float64}
    kf::Vector{Float64}
    kr::Vector{Float64}
    wdot::Vector{Float64}
    qdot::Vector{Float64}
end

function CreateMSolution(gas)
    ns = gas.n_species
    Y = zeros(ns)
    Y[1] = 1.
    X = zeros(ns)
    C = zeros(ns)
    T = 1200.
    P = one_atm
    mean_molecular_weight = 0.
    ρ_mass = 0.
    cp_mole = 0.
    cp_mass = 0.
    H_mass = 0.
    H_mole = 0.
    h_mole = zeros(ns)
    S_mass = 0.
    S_mole = 0.
    s_mole = zeros(ns)
    S0 = zeros(ns)
    kf = zeros(gas.n_reactions)
    kr = zeros(gas.n_reactions)
    wdot = zeros(ns)
    qdot = zeros(gas.n_reactions)
    mgas = MSolution(Y, X, C, T, P, mean_molecular_weight, ρ_mass, cp_mole, cp_mass,
                    H_mass, H_mole, h_mole, S_mass, S_mole, s_mole, S0, kf, kr, wdot, qdot)
end

function set_states(gas, mgas, T, P, Y)
    mgas.T = T
    mgas.P = P
    mgas.Y = Y
    mgas.mean_molecular_weight = 1. / dot(mgas.Y, 1 ./ gas.molecular_weights)
    mgas.ρ_mass = mgas.P / R / mgas.T * mgas.mean_molecular_weight
    Y2X(gas, mgas)
    Y2C(gas, mgas)
    cp_mole_func(gas, mgas)
    cp_mass_func(gas, mgas)
    H_mole_func(gas, mgas)
    H_mass_func(gas, mgas)
    S_mole_func(gas, mgas)
    S_mass_func(gas, mgas)
    wdot_func(gas, mgas)
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
    reactant_orders = npz["reactant_orders"]
    is_reversible = npz["is_reversible"]
    list_type4_noTroe = npz["list_type4_noTroe"]
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

    reaction = Reaction(product_stoich_coeffs, reactant_stoich_coeffs, reactant_orders, is_reversible, list_type4_noTroe,
                        Arrhenius_coeffs, Arrhenius_A0, Arrhenius_b0, Arrhenius_Ea0,
                        Troe_A, Troe_T1, Troe_T2, Troe_T3, index_three_body, index_falloff, efficiencies_coeffs)

    gas = Solution(n_species, n_reactions, molecular_weights, species_names, thermo, reaction)
end
