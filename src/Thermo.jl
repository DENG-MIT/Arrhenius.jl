# include all different thermo files here 
#TODO Automate it to include all jl files in the Thermo folder
include("Thermo/IdealGasThermo.jl")

# --- Allocation-free mixture means for IdealGasThermo ---
# Which of the dimensionless cp_R, h_RT, s0_R each mixture property accumulates.
@inline _mean_needs(::Val{:cp}) = (true, false, false)
@inline _mean_needs(::Val{:cv}) = (true, false, false)
@inline _mean_needs(::Val{:h}) = (false, true, false)
@inline _mean_needs(::Val{:u}) = (false, true, false)
@inline _mean_needs(::Val{:s0}) = (false, false, true)
@inline _mean_needs(::Val{:s}) = (false, false, true)
@inline _mean_needs(::Val{:g}) = (false, true, true)
@inline _mean_needs(::Val{:a}) = (false, true, true)

# Per-species contribution to the dimensionless mixture mean, matching the
# per-species cal_* definitions term by term, including the log(max(X, tiny))
# entropy cutoff and the one_atm pressure reference.
@inline _mean_term(::Val{:cp}, Xi, cp, h, s, tiny, logp) = Xi * cp
@inline _mean_term(::Val{:cv}, Xi, cp, h, s, tiny, logp) = Xi * (cp - 1)
@inline _mean_term(::Val{:h}, Xi, cp, h, s, tiny, logp) = Xi * h
@inline _mean_term(::Val{:u}, Xi, cp, h, s, tiny, logp) = Xi * (h - 1)
@inline _mean_term(::Val{:s0}, Xi, cp, h, s, tiny, logp) = Xi * s
@inline _mean_term(::Val{:s}, Xi, cp, h, s, tiny, logp) = Xi * ((s - log(max(Xi, tiny))) - logp)
@inline _mean_term(::Val{:g}, Xi, cp, h, s, tiny, logp) = Xi * (h - ((s - log(max(Xi, tiny))) - logp))
@inline _mean_term(::Val{:a}, Xi, cp, h, s, tiny, logp) = Xi * ((h - 1) - ((s - log(max(Xi, tiny))) - logp))

function _nasa7_mean_dimless(::Val{phi}, thermo::IdealGasThermo,
                             T::Real, p::Real, X::AbstractArray) where {phi}
    need_cp, need_h, need_s = _mean_needs(Val(phi))
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    logT = need_s ? log(T) : zero(T)
    tiny = oftype(T, 1.0e-30)
    logp = phi in (:s,:g,:a) ? log(p / oftype(T, one_atm)) : zero(T)
    z = zero(promote_type(typeof(T), eltype(thermo.nasa_low)))
    acc = zero(promote_type(typeof(T), eltype(thermo.nasa_low), eltype(X), typeof(logp)))
    @inbounds for i in eachindex(X)
        Xi = X[i]
        # Exact zeros contribute nothing for ordinary floating compositions.
        # Dual-valued fractions retain their species term: their composition
        # derivative can be nonzero even when the primal fraction is zero.
        Xi isa AbstractFloat && iszero(Xi) && continue
        nasa = T <= thermo.Trange[i, 2] ? thermo.nasa_low : thermo.nasa_high
        cp = need_cp ? nasa[i, 1] + nasa[i, 2] * T + nasa[i, 3] * T2 +
             nasa[i, 4] * T3 + nasa[i, 5] * T4 : z
        h = need_h ? nasa[i, 1] + nasa[i, 2] * T / 2 + nasa[i, 3] * T2 / 3 +
            nasa[i, 4] * T3 / 4 + nasa[i, 5] * T4 / 5 + nasa[i, 6] / T : z
        s = need_s ? nasa[i, 1] * logT + nasa[i, 2] * T + nasa[i, 3] * T2 / 2 +
            nasa[i, 4] * T3 / 3 + nasa[i, 5] * T4 / 4 + nasa[i, 7] : z
        acc += _mean_term(Val(phi), Xi, cp, h, s, tiny, logp)
    end
    return acc
end

function _extended_mean_dimless(::Val{phi}, thermo::IdealGasThermo,
                                T::Real, p::Real, X::AbstractArray) where {phi}
    need_cp, need_h, need_s = _mean_needs(Val(phi))
    tiny = oftype(T, 1.0e-30)
    logp = phi in (:s,:g,:a) ? log(p / oftype(T, one_atm)) : zero(T)
    z = zero(promote_type(typeof(T), eltype(thermo.nasa_low)))
    acc = zero(promote_type(typeof(T), eltype(thermo.nasa_low), eltype(X), typeof(logp)))
    @inbounds for i in eachindex(X)
        Xi = X[i]
        Xi isa AbstractFloat && iszero(Xi) && continue
        species_cp, species_h, species_s = _extended_thermo(thermo, i, T)
        cp = need_cp ? species_cp : z
        h = need_h ? species_h : z
        s = need_s ? species_s : z
        acc += _mean_term(Val(phi), Xi, cp, h, s, tiny, logp)
    end
    return acc
end

"Dimensionless ideal-gas mixture mean `dot(X, cal_phi_dimless)` without species vectors."
function _mixture_mean_dimless(::Val{phi}, gas::Arrhenius.Solution, thermo::IdealGasThermo,
                               T::Real, p::Real, X::AbstractArray) where {phi}
    length(X) == gas.n_species || throw(DimensionMismatch(
        "expected $(gas.n_species) species fractions, got $(length(X))"))
    if isnothing(thermo.extra)
        return _nasa7_mean_dimless(Val(phi), thermo, T, p, X)
    end
    return _extended_mean_dimless(Val(phi), thermo, T, p, X)
end

# Metaprogramming loop to generate and export all mass and mean functions
property_names=((:cv,"Heat capacity at constant volume (cv)"),
                (:cp,"Heat capacity at constant pressure (cp)"), 
                (:s,"entropy (s)"), 
                (:s0,"reference entropy (s0)"), 
                (:h,"enthalpy (h)"),
                (:a,"helmholz free energy (a)"), 
                (:g,"gibbs free energy (g)"), 
                (:u,"internal energy"))
for (phi, doc_name) in property_names
    # decide if dimless factor is R or RT
    if phi in (:cv, :cp, :s, :s0)
        dimmless = Symbol(:_R)
        RRT = :(oftype(T, R))
    else
        dimmless = Symbol(:_RT)
        RRT = :(oftype(T, R) * T)
    end
    #define the 5 functions for each quantity
    cal_phi_dimless = Symbol(:cal_,phi,dimmless)
    cal_phi = Symbol(:cal_,phi)
    cal_phimass =  Symbol(:cal_,phi,:mass)
    cal_phi_mean = Symbol(cal_phi,:_mean)
    cal_phimass_mean = Symbol(cal_phimass,:_mean)

    @eval begin 
        # Dispatches the call from a solution object to it's thermo object
        $cal_phi_dimless(gas::Arrhenius.Solution, T, p, X)=$cal_phi_dimless(gas,gas.thermo, T, p, X)
        export $cal_phi_dimless 
        """
            $($cal_phi)(Solution, T, p, X)

        calculates the molar $($doc_name) for each species
        """
        function $cal_phi(gas::Arrhenius.Solution,thermo::Thermo,
                           T::Real, p::Real, X::AbstractArray)
            return $cal_phi_dimless(gas,thermo, T, p, X) * $RRT
        end
        $cal_phi(gas::Arrhenius.Solution, T, p, X)=$cal_phi(gas,gas.thermo, T, p, X)
        export $cal_phi 
        """
            $($cal_phi_mean)(Solution, T, p, X)
        
        calculates the mean mole based $($doc_name) of the mixture
        """
        function $cal_phi_mean(gas::Arrhenius.Solution,thermo::Thermo,
                                T::Real, p::Real, X::AbstractArray)
            return dot(X,$cal_phi_dimless(gas,thermo, T, p, X)) * $RRT
        end
        # Allocation-free specialization: accumulates the same mean directly.
        function $cal_phi_mean(gas::Arrhenius.Solution,thermo::IdealGasThermo,
                                T::Real, p::Real, X::AbstractArray)
            return _mixture_mean_dimless(Val{$(QuoteNode(phi))}(), gas, thermo, T, p, X) * $RRT
        end
        $cal_phi_mean(gas::Arrhenius.Solution, T, p, X)=$cal_phi_mean(gas,gas.thermo, T, p, X)
        export $cal_phi_mean 
        """
            $($cal_phimass)(Solution, T, p, X)

        calculates the partial mass based $($doc_name) for each species
        """
        function $cal_phimass(gas::Arrhenius.Solution,thermo::Thermo,
                              T::Real, p::Real, X::AbstractArray)
            return $cal_phi(gas,thermo, T, p, X) ./gas.MW
        end
        $cal_phimass(gas::Arrhenius.Solution, T, p, X)=$cal_phimass(gas,gas.thermo, T, p, X)
        export $cal_phimass 
        """
            $($cal_phimass_mean)(Solution, T, p, X)

        calculates the mean mass based $($doc_name) of the mixture
        """
        function $cal_phimass_mean(gas::Arrhenius.Solution,thermo::Thermo,
                                    T::Real, p::Real, X::AbstractArray)
            return $cal_phi_mean(gas,thermo, T, p, X) / dot(X,gas.MW)
        end        
        $cal_phimass_mean(gas::Arrhenius.Solution, T, p, X)=$cal_phimass_mean(gas,gas.thermo, T, p, X)
        export $cal_phimass_mean 
    end
end


# ------ Deprecated ----------# 

"get specific of heat capacity"
function get_cp(gas, T, X, mean_MW)
    cp = cal_cp_R(gas, T, one_atm, X)
    cp_mole = dot(cp, X) * oftype(T, R)
    cp_mass = cp_mole / mean_MW
    return cp_mole, cp_mass
end
export get_cp


"get specific of heat capacity"
function get_cv(cp_mole, cp_mass, mean_MW)
    cv_mole = cp_mole - oftype(cp_mole, R)
    cv_mass = cv_mole / mean_MW
    return cv_mole, cv_mass
end
export get_cv


"get enthaphy (H) per mole"
function get_H(gas, T, Y, X)
    return cal_h_RT(gas, T, one_atm, X) * oftype(T, R) * T
end
export get_H


"get enthaphy (H) per mass"
function H_mass_func(gas, h_mole, Y)
    return dot(h_mole ./ gas.MW, Y)
end
export H_mass_func


"get enthaphy (U) per mole"
function get_U(h_mole, T)
    u_mole = h_mole .- (oftype(T, R) * T)
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
    return cal_s0_R(gas, T, P, X) * oftype(T, R)
end
export get_S

"get entropy (S) per unit mass"
function S_mass_func(gas, s_mole, Y)
    return dot(s_mole ./ gas.MW, Y)
end
export S_mass_func
