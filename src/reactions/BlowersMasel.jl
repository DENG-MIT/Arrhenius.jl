"""
    BlowersMaselRate(A, b, Ea0, w)

Enthalpy-dependent modified Arrhenius rate. `Ea0` and the average bond energy
`w` use J/kmol, temperature uses K, and the units of `A` follow reaction order.
The intrinsic barrier must be nonnegative and `w > Ea0`.
"""
struct BlowersMaselRate{T<:AbstractFloat}
    A::T
    b::T
    Ea0::T
    w::T
    function BlowersMaselRate(A::T,b::T,Ea0::T,w::T) where {T<:AbstractFloat}
        all(isfinite,(A,b,Ea0,w)) && Ea0 >= 0 && w > Ea0 ||
            throw(ArgumentError("finite coefficients, nonnegative intrinsic barrier and w > Ea0 required"))
        new{T}(A,b,Ea0,w)
    end
end
BlowersMaselRate(A,b,Ea0,w) = BlowersMaselRate(promote(float(A),float(b),float(Ea0),float(w))...)

# Blowers and Masel, AIChE Journal 46 (2000), doi:10.1002/aic.690461015.
# https://cantera.org/dev/reference/kinetics/rate-constants.html#blowers-masel-reactions
function _blowers_masel_barrier(Ea0,w,delta_h)
    delta_h <= -4Ea0 && return zero(delta_h)
    delta_h > 4Ea0 && return delta_h
    vp = 2w*(w+Ea0)/(w-Ea0)
    return (w+delta_h/2)*(vp-2w+delta_h)^2/(vp^2-4w^2+delta_h^2)
end

"Effective activation energy [J/kmol] at the supplied reaction enthalpy [J/kmol]."
activation_energy(rate::BlowersMaselRate,delta_h) = _blowers_masel_barrier(rate.Ea0,rate.w,delta_h)

"Evaluate a Blowers–Masel forward rate constant from temperature and reaction enthalpy."
function rate_constant(rate::BlowersMaselRate,T,delta_h)
    isfinite(T) && T > 0 && isfinite(delta_h) || throw(ArgumentError("finite positive temperature and finite enthalpy required"))
    return rate.A*exp(rate.b*log(T)-activation_energy(rate,delta_h)/(R*T))
end

struct BlowersMaselData{T<:AbstractFloat}
    reaction_indices::Vector{Int64}
    coefficients::Matrix{T}
end

export BlowersMaselRate, activation_energy, rate_constant
