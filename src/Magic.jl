function Y2X(Y, mean_MW)
    return Y * mean_MW ./ gas.MW
end

function Y2C(Y, ρ_mass)
    return Y * ρ_mass ./ gas.MW
end

function C2X(C)
    return C ./ sum(C)
end

function X2C(X, ρ_mass)
    return X * ρ_mass ./ gas.MW
end

function X2Y(X, mean_MW)
    return X .* gas.MW / mean_MW
end

function species_index(species)
    return findfirst(gas.species_names .== species)
end
