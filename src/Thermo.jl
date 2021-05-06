
"calculates the dimensionless partial enhtalpy of each species"
function cal_h_RT(gas, T, p, X)
    H_T = [1.0, T / 2.0, T^2 / 3.0, T^3 / 4.0, T^4 / 5.0, 1.0 / T]
    if T <= 1000.0
        h_mole = @view(gas.thermo.nasa_low[:, 1:6]) * H_T 
    else
        h_mole = @view(gas.thermo.nasa_high[:, 1:6]) * H_T 
    end
    if !gas.thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < gas.thermo.Trange[:, 2])
        h_mole[ind_correction] .=
            @view(gas.thermo.nasa_low[ind_correction, 1:6]) * H_T 
    end
    # H_mole = dot(h_mole, X)
    return h_mole
end

"calculates the dimensionless entropy of each species at the reference state"
function cal_s0_R(gas, T,p, X)
    S_T = [log(T), T, T^2 / 2.0, T^3 / 3.0, T^4 / 4.0, 1.0]
    if T <= 1000.0
        S0 = @view(gas.thermo.nasa_low[:, [1, 2, 3, 4, 5, 7]]) * S_T 
    else
        S0 = @view(gas.thermo.nasa_high[:, [1, 2, 3, 4, 5, 7]]) * S_T 
    end
    if !gas.thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < gas.thermo.Trange[:, 2])
        S0[ind_correction] .=
            @view(gas.thermo.nasa_low[ind_correction, [1, 2, 3, 4, 5, 7]]) *
            S_T 
    end
    return S0 
end
cal_s0_R(gas,T)=cal_s0_R(gas, T, 0, [])
export cal_s0_R

"calculates the dimensionless entropy of each species"
function cal_s_R(gas, T, p, X)
   return cal_s0_R(gas, T) - log.(max.(X,1e-30)) .- log(p/one_atm)
end

"calculates the dimensionless Gibbs energy of each species"
function cal_g_RT(gas, T, p, X)
    return cal_h_RT(gas, T, p, X) - cal_s_R(gas, T, p, X)
end

"calculates the dimensionless internal energy of each species"
function cal_u_RT(gas, T, p, X)
    return cal_h_RT(gas, T, p, X) .- 1
end

"calculates the dimensionless Helmholz free energy of each species"
function cal_a_RT(gas, T, p, X)
    return cal_u_RT(gas, T, p, X) - cal_s_R(gas, T, p, X)
end
"calculates the dimensionless heat capacity at constant pressure of each species"
function cal_cp_R(gas, T, p, X)
    cp_T = [1.0, T, T^2, T^3, T^4]
    if T <= 1000.0
        cp = @view(gas.thermo.nasa_low[:, 1:5]) * cp_T
    else
        cp = @view(gas.thermo.nasa_high[:, 1:5]) * cp_T
    end
    # TODO: not sure if inplace operation will be an issue for AD
    if !gas.thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < gas.thermo.Trange[:, 2])
        cp[ind_correction] .=
            @view(gas.thermo.nasa_low[ind_correction, 1:5]) * cp_T
    end
    return cp
end

"calculates the dimensionless heat capacity at constant volume of each species"
function cal_cv_R(gas, T, p, X)
   return cal_cp_R(gas, T, p, X) .- 1
end

# Metaprogramming loop to generate and export all mass and mean functions
for phi in [:cv, :cp, :s, :s0, :h, :a, :g, :u]
    if phi in (:cv, :cp, :s, :s0)
        dimmless = Symbol(:_R)
        RRT =:(R)
    else
        dimmless = Symbol(:_RT)
        RRT =:(R*T)
    end
    cal_phi_dimless = Symbol(:cal_,phi,dimmless)
    cal_phi = Symbol(:cal_,phi)
    cal_phimass =  Symbol(:cal_,phi,:mass)
    cal_phi_mean = Symbol(cal_phi,:_mean)
    cal_phimass_mean = Symbol(cal_phimass,:_mean)

    @eval begin 
        export $cal_phi_dimless
        "calculates the partial molar phi for each species"
        function $cal_phi(gas, T, p, X)
            return $cal_phi_dimless(gas, T, p, X) * $RRT
        end
        export $cal_phi
        "calculates the mean mole based phi of the mixture"
        function $cal_phi_mean(gas, T, p, X)
            return dot(X,$cal_phi_dimless(gas, T, p, X)) * $RRT
        end
        export $cal_phi_mean
        "calculates the partial mass based phi for each species"
        function $cal_phimass(gas, T, p, X)
            return $cal_phi(gas, T, p, X) ./gas.MW
        end
        export $cal_phimass
        "calculates the mean mass based phi of the mixture"
        function $cal_phimass_mean(gas, T, p, X)
            return $cal_phi_mean(gas, T, p, X) / dot(X,gas.MW)
        end
        export $cal_phimass_mean
    end
end


# ------ Deprecated ----------# 

"get specific of heat capacity"
function get_cp(gas, T, X, mean_MW)
    cp_T = [1.0, T, T^2, T^3, T^4]
    if T <= 1000.0
        cp = @view(gas.thermo.nasa_low[:, 1:5]) * cp_T
    else
        cp = @view(gas.thermo.nasa_high[:, 1:5]) * cp_T
    end
    # TODO: not sure if inplace operation will be an issue for AD
    if !gas.thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < gas.thermo.Trange[:, 2])
        cp[ind_correction] .=
            @view(gas.thermo.nasa_low[ind_correction, 1:5]) * cp_T
    end
    cp_mole = dot(cp, X) * R
    cp_mass = cp_mole / mean_MW
    return cp_mole, cp_mass
end
export get_cp


"get specific of heat capacity"
function get_cv(cp_mole, cp_mass, mean_MW)
    cv_mole = cp_mole - R
    cv_mass = cv_mole / mean_MW
    return cv_mole, cv_mass
end
export get_cv


"get enthaphy (H) per mole"
function get_H(gas, T, Y, X)
    H_T = [1.0, T / 2.0, T^2 / 3.0, T^3 / 4.0, T^4 / 5.0, 1.0 / T]
    if T <= 1000.0
        h_mole = @view(gas.thermo.nasa_low[:, 1:6]) * H_T * R * T
    else
        h_mole = @view(gas.thermo.nasa_high[:, 1:6]) * H_T * R * T
    end
    if !gas.thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < gas.thermo.Trange[:, 2])
        h_mole[ind_correction] .=
            @view(gas.thermo.nasa_low[ind_correction, 1:6]) * H_T * R * T
    end
    # H_mole = dot(h_mole, X)
    return h_mole
end
export get_H


"get enthaphy (H) per mass"
function H_mass_func(gas, h_mole, Y)
    return dot(h_mole ./ gas.MW, Y)
end
export H_mass_func


"get enthaphy (U) per mole"
function get_U(h_mole, T)
    u_mole = h_mole .- (R * T)
    return u_mole
end
export get_U


"get enthaphy (U) per mass"
function U_mass_func(gas, u_mole, Y)
    return dot(u_mole ./ gas.MW, Y)
end
export U_mass_func


"get entropy (S)"
function get_S(gas, T, P, X)
    S_T = [log(T), T, T^2 / 2.0, T^3 / 3.0, T^4 / 4.0, 1.0]
    if T <= 1000.0
        S0 = @view(gas.thermo.nasa_low[:, [1, 2, 3, 4, 5, 7]]) * S_T * R
    else
        S0 = @view(gas.thermo.nasa_high[:, [1, 2, 3, 4, 5, 7]]) * S_T * R
    end
    if !gas.thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < gas.thermo.Trange[:, 2])
        S0[ind_correction] .=
            @view(gas.thermo.nasa_low[ind_correction, [1, 2, 3, 4, 5, 7]]) *
            S_T * R
    end
    # _X = @. S0 - R * log(clamp(X, 1.e-30, Inf))
    # s_mole = _X .- R * (P / one_atm)
    # S_mole = dot(s_mole, X)
    return S0
end
export get_S

"get entropy (S) per unit mass"
function S_mass_func(gas, s_mole, Y)
    return dot(s_mole ./ gas.MW, Y)
end
export S_mass_func
