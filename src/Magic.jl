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
