function cp_mole_func(gas, mgas)
    T = mgas.T
    cp_T = [1., T, T^2, T^3, T^4]
    if T <= 1000.
        cp = @view(gas.thermo.nasa_low[:, 1:5]) * cp_T
    else
        cp = @view(gas.thermo.nasa_high[:, 1:5]) * cp_T
    end
    mgas.cp_mole = dot(cp, mgas.X) * R
end

function cp_mass_func(gas, mgas)
    mgas.cp_mass = mgas.cp_mole / mgas.mean_molecular_weight
end

function H_mole_func(gas, mgas)
    T = mgas.T
    H_T = [1., T / 2., T^2 / 3., T^3 / 4., T^4 / 5., 1. /T]
    if T <= 1000.
        mgas.h_mole = @view(gas.thermo.nasa_low[:, 1:6]) * H_T  * R * T
    else
        mgas.h_mole = @view(gas.thermo.nasa_high[:, 1:6]) * H_T  * R * T
    end
    mgas.H_mole = dot(mgas.h_mole, mgas.X)
end

function H_mass_func(gas, mgas)
    mgas.H_mass = dot(mgas.h_mole ./ gas.molecular_weights, mgas.Y)
end

function S_mole_func(gas, mgas)
    T = mgas.T
    S_T = [log(T), T, T^2 / 2., T^3 / 3., T^4 / 4., 1.]
    if T <= 1000.
        S0 = @view(gas.thermo.nasa_low[:, [1,2,3,4,5,7]]) * S_T  * R
    else
        S0 = @view(gas.thermo.nasa_high[:, [1,2,3,4,5,7]]) * S_T  * R
    end
    mgas.S0 = S0
    X = @. S0 - R * log(clamp(mgas.X, 1.e-30, Inf))
    mgas.s_mole = X .- R * (mgas.P / one_atm)
    mgas.S_mole = dot(mgas.s_mole, mgas.X)
end

function S_mass_func(gas, mgas)
    mgas.S_mass = dot(mgas.s_mole ./ gas.molecular_weights, mgas.Y)
end
