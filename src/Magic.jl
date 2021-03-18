"get mole fraction (X) from mass fraction (Y)"
function Y2X(gas, Y, mean_MW)
    return Y * mean_MW ./ gas.MW
end
export Y2X

"get concentration (C) from mass fraction (Y)"
function Y2C(gas, Y, ρ_mass)
    return Y * ρ_mass ./ gas.MW
end
export Y2C

"get mole fraction (X) from concentration (C)"
function C2X(gas, C)
    return C ./ sum(C)
end
export C2X

"get concentration (C) from mole fraction (X)"
function X2C(gas, X, ρ_mass)
    return X * ρ_mass ./ gas.MW
end
export X2C

"get mass fraction (Y) from mole fraction (X)"
function X2Y(gas, X, mean_MW)
    return X .* gas.MW / mean_MW
end
export X2Y

"""
Get species index of a species
# Example
```
    species_index(gas, "O2")
```
"""
function species_index(gas, species)
    return findfirst(gas.species_names .== species)
end
export species_index

"""
"mainly for testing code, will be removed in the future"
"please customize such functions following this example"
"""
function set_states(gas::Solution, T, P, Y)
    mean_MW = 1.0 / dot(Y, 1 ./ gas.MW)
    ρ_mass = P / R / T * mean_MW
    X = Y2X(gas, Y, mean_MW)
    C = Y2C(gas, Y, ρ_mass)
    cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
    h_mole = get_H(gas, T, Y, X)
    S0 = get_S(gas, T, P, X)
    wdot = wdot_func(gas.reaction, T, C, S0, h_mole)
    return wdot
end
export set_states