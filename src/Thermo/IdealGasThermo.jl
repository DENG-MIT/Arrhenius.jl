"""Struct for the ideal gas thermo.

nasa_low: Array with low temperature nasa coeff. for each species

nasa_high: Array with high temperature nasa coeff. for each species

Trange: Array with temperature ranges for each species

isTcommon: bool which indicates if both polynoms share same T at intersection

"""
struct IdealGasThermo <: Thermo 
    nasa_low::Array{Float64,2}
    nasa_high::Array{Float64,2}
    Trange::Array{Float64,2}
    isTcommon::Bool

end

"""
    cal_h_RT(gas, T, p, X)

calculates the dimensionless mole based enthalpy (h) for each species
"""
function cal_h_RT(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    H_T = [1.0, T / 2.0, T^2 / 3.0, T^3 / 4.0, T^4 / 5.0, 1.0 / T]
    if T <= 1000.0
        h_mole = @view(thermo.nasa_low[:, 1:6]) * H_T 
    else
        h_mole = @view(thermo.nasa_high[:, 1:6]) * H_T 
    end
    if !thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < thermo.Trange[:, 2])
        h_mole[ind_correction] .=
            @view(thermo.nasa_low[ind_correction, 1:6]) * H_T 
    end
    # H_mole = dot(h_mole, X)
    return h_mole
end

"""
    cal_s0_R(gas, T, p, X)

calculates the dimensionless mole based reference state entropy (s0) for each species
"""
function cal_s0_R(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    S_T = [log(T), T, T^2 / 2.0, T^3 / 3.0, T^4 / 4.0, 1.0]
    if T <= 1000.0
        S0 = @view(thermo.nasa_low[:, [1, 2, 3, 4, 5, 7]]) * S_T 
    else
        S0 = @view(thermo.nasa_high[:, [1, 2, 3, 4, 5, 7]]) * S_T 
    end
    if !thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < thermo.Trange[:, 2])
        S0[ind_correction] .=
            @view(thermo.nasa_low[ind_correction, [1, 2, 3, 4, 5, 7]]) *
            S_T 
    end
    return S0 
end
export cal_s0_R

"""
    cal_s_R(gas, T, p, X)

calculates the dimensionless mole based entropy (s) for each species
"""
function cal_s_R(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
   return cal_s0_R(gas,thermo, T,p,X) - log.(max.(X,1e-30)) .- log(p/one_atm)
end

"""
    cal_g_RT(gas, T, p, X)

calculates the dimensionless mole based free gibbs energy (g) for each species
"""
function cal_g_RT(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    return cal_h_RT(gas,thermo, T, p, X) - cal_s_R(gas,thermo, T, p, X)
end

"""
    cal_u_RT(gas, T, p, X)

calculates the dimensionless mole based internal energy (u) for each species
"""
function cal_u_RT(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    return cal_h_RT(gas,thermo, T, p, X) .- 1
end

"""
    cal_a_RT(gas, T, p, X)

calculates the dimensionless mole based helmholz free energy (a) for each species
"""
function cal_a_RT(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    return cal_u_RT(gas,thermo, T, p, X) - cal_s_R(gas,thermo, T, p, X)
end

"""
    cal_cp_R(gas, T, p, X)

calculates the dimensionless mole based heat capacity 
at constant pressure (cp) for each species
"""
function cal_cp_R(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    cp_T = [1.0, T, T^2, T^3, T^4]
    if T <= 1000.0
        cp = @view(thermo.nasa_low[:, 1:5]) * cp_T
    else
        cp = @view(thermo.nasa_high[:, 1:5]) * cp_T
    end
    # TODO: not sure if inplace operation will be an issue for AD
    if !thermo.isTcommon
        ind_correction = @. (T > 1000.0) & (T < thermo.Trange[:, 2])
        cp[ind_correction] .=
            @view(thermo.nasa_low[ind_correction, 1:5]) * cp_T
    end
    return cp
end

"""
    cal_cv_R(gas, T, p, X)

calculates the dimensionless mole based heat capacity 
at constant volume (cv) for each species
"""
function cal_cv_R(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
   return cal_cp_R(gas,thermo, T, p, X) .- 1
end

