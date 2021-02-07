function Y2X(gas, mgas)
    mgas.X = mgas.Y * mgas.mean_molecular_weight ./ gas.molecular_weights
end

function Y2C(gas, mgas)
    mgas.C = mgas.Y * mgas.ρ_mass ./ gas.molecular_weights
end

function C2X(gas, mgas)
    mgas.X = mgas.C ./ sum(mgas.C)
end

function X2C(gas, mgas)
    mgas.C = mgas.X * mgas.ρ_mass ./ gas.molecular_weights
end

function X2Y(gas, mgas)
    mgas.Y = mgas.X .* gas.molecular_weights / mgas.mean_molecular_weight
end
