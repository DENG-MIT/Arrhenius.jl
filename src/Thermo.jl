function cp_mole_func(gas, T, X)
    if T <= 1000.0
        cp = @view(gas.thermo.nasa_low[:, 1:5]) * [1.0, T, T^2, T^3, T^4]
    else
        cp = @view(gas.thermo.nasa_high[:, 1:5]) * [1.0, T, T^2, T^3, T^4]
    end
    return dot(cp, X) * R
end

function cp_mass_func(gas, cp_mole, mean_MW)
    return cp_mole / mean_MW
end

function H_mole_func(gas, T, X)
    H_T = [1.0, T / 2.0, T^2 / 3.0, T^3 / 4.0, T^4 / 5.0, 1.0 / T]
    if T <= 1000.0
        h_mole = @view(gas.thermo.nasa_low[:, 1:6]) * H_T * R * T
    else
        h_mole = @view(gas.thermo.nasa_high[:, 1:6]) * H_T * R * T
    end
    H_mole = dot(h_mole, X)
    return H_mole, h_mole
end

function H_mass_func(gas, h_mole, Y)
    return dot(h_mole ./ gas.MW, Y)
end

function S_mole_func(gas, T, P, X)
    S_T = [log(T), T, T^2 / 2.0, T^3 / 3.0, T^4 / 4.0, 1.0]
    if T <= 1000.0
        S0 = @view(gas.thermo.nasa_low[:, [1, 2, 3, 4, 5, 7]]) * S_T * R
    else
        S0 = @view(gas.thermo.nasa_high[:, [1, 2, 3, 4, 5, 7]]) * S_T * R
    end
    _X = @. S0 - R * log(clamp(X, 1.e-30, Inf))
    s_mole = _X .- R * (P / one_atm)
    S_mole = dot(s_mole, X)
    return S_mole, s_mole, S0
end

function S_mass_func(gas, s_mole, Y)
    return dot(s_mole ./ gas.MW, Y)
end
