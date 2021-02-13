function Y2X(gas, Y, mean_MW)
    return Y * mean_MW ./ gas.MW
end

function Y2C(gas, Y, ρ_mass)
    return Y * ρ_mass ./ gas.MW
end

function C2X(gas, C)
    return C ./ sum(C)
end

function X2C(gas, X, ρ_mass)
    return X * ρ_mass ./ gas.MW
end

function X2Y(gas, X, mean_MW)
    return X .* gas.MW / mean_MW
end

function species_index(gas, species)
    return findfirst(gas.species_names .== species)
end
