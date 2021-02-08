struct Reaction
    product_stoich_coeffs::Array{Float64}
    reactant_stoich_coeffs::Array{Float64}
    reactant_orders::Array{Float64, 2}
    is_reversible::Array{Bool, 1}
    list_type4_noTroe::Array{Int32, 1}
    Arrhenius_coeffs::Array{Float64, 2}
    Arrhenius_A0::Array{Float64, 1}
    Arrhenius_b0::Array{Float64, 1}
    Arrhenius_Ea0::Array{Float64, 1}
    Troe_A::Array{Float64, 1}
    Troe_T1::Array{Float64, 1}
    Troe_T2::Array{Float64, 1}
    Troe_T3::Array{Float64, 1}
    index_three_body::Array{Int32, 1}
    index_falloff::Array{Int32, 1}
    efficiencies_coeffs::Array{Float64, 2}
    i_reactant::Array{Array{Int32, 1},1}
    i_product::Array{Array{Int32, 1},1}
end

struct Thermo
    nasa_low::Array{Float64, 2}
    nasa_high::Array{Float64, 2}
end

struct Solution
    n_species::Int32
    n_reactions::Int32
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
