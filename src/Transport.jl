"""
    mix_trans(gas::A, P, T, X, mean_MW) where {A <: Arrhenius.Solution}

Compute the tranposrt properties of a mixture using mixture average formula

> Equations Ref. https://personal.ems.psu.edu/~radovic/ChemKin_Theory_PaSR.pdf
> Equations. 5-50/51/52 for viscosity and thermal conductivity

Pure species viscosities [Pa-s]

Thermal conductivity. [W/m/K].

> Equation 5-46 for diffusion

Mixture-averaged diffusion coefficients [m^2/s] relating the mass-averaged diffusive fluxes 
(with respect to the mass averaged velocity) to gradients in the species mole fractions.

Test this module in _transport_test.jl

See also implementations in ReacTorch
"""
function mix_trans(gas::A, P, T, X, mean_MW) where {A <: Arrhenius.Solution}
    gas.trans.model == :ionized_gas && throw(ArgumentError("use ionized_transport! for ionized-gas transport"))

    if gas.trans.poly_order == 5
        workspace = TransportWorkspace(gas)
        viscosity, conductivity = mixture_transport!(workspace, gas, P, T, X)
        return viscosity, conductivity, workspace.diffusion
    end
    size(gas.trans.species_viscosities_poly) == (7, gas.n_species) ||
        throw(ArgumentError("transport data missing; regenerate the mechanism sidecar"))

    logT = log(T)
    trans_T = [logT^6, logT^5, logT^4, logT^3, logT^2, logT, 1.0]
    
    ## species_viscosities_poly

    η = trans_T' * gas.trans.species_viscosities_poly

    Wk_over_Wj = gas.MW * (1.0 ./ gas.MW')

    Wj_over_Wk = 1 ./ Wk_over_Wj

    ηk_over_ηj = η' * (1.0 ./ η)

    Φ = @. 1.0 / sqrt(8.0) / sqrt(1.0 + Wk_over_Wj) * 
        (1.0 + sqrt(ηk_over_ηj) * (Wj_over_Wk)^(0.25))^2

    η_mix = sum((X .* η') ./ (Φ * X))

    ## thermal_conductivity_func
    λ = trans_T' * gas.trans.thermal_conductivity_poly

    λ_mix = (sum(X .* λ') + 1.0 / sum(X ./ λ')) / 2.0

    ## binary_diff_coeffs_func

    X_eps = clamp.(X, 1.e-12, Inf)
    D = reshape(trans_T' * gas.trans.binary_diff_coeffs_poly, 
                gas.n_species, gas.n_species)

    XjWj = X_eps' * gas.MW .- X_eps .* gas.MW

    XjDjk = sum(X_eps .* (1.0 ./ D'), dims=1) .- (X_eps ./ diag(D))'

    Dkm = XjWj ./ XjDjk' ./ mean_MW / P * one_atm

    return η_mix, λ_mix, Dkm
end
export mix_trans

"Reusable storage for mixture-averaged ideal-gas transport properties."
# `inverse_binary` holds reciprocal binary diffusion coefficients; the two pair-weight
# matrices are stored transposed ([j,k] for species k) for column access.
struct TransportWorkspace
    viscosity::Vector{Float64}
    conductivity::Vector{Float64}
    diffusion::Vector{Float64}
    binary::Matrix{Float64}
    weight_quarter_ratio::Matrix{Float64}
    inverse_phi_denominator::Matrix{Float64}
    molecular_weights::Vector{Float64}
    inverse_binary::Matrix{Float64}
    root_viscosity::Vector{Float64}
    inverse_root_viscosity::Vector{Float64}
    weighted_mole_fraction::Vector{Float64}
end
TransportWorkspace(gas::Solution) = TransportWorkspace(
    zeros(gas.n_species), zeros(gas.n_species), zeros(gas.n_species),
    zeros(gas.n_species, gas.n_species),
    [sqrt(sqrt(gas.MW[j]/gas.MW[k])) for j in 1:gas.n_species,k in 1:gas.n_species],
    [inv(sqrt(8*(1+gas.MW[k]/gas.MW[j]))) for j in 1:gas.n_species,k in 1:gas.n_species],
    copy(gas.MW), zeros(gas.n_species, gas.n_species),
    zeros(gas.n_species), zeros(gas.n_species), zeros(gas.n_species),
)

@inline function _transport_polynomial(coefficients, k, logT)
    return evalpoly(logT, (coefficients[1,k], coefficients[2,k],
        coefficients[3,k], coefficients[4,k], coefficients[5,k]))
end

"""
    mixture_transport!(workspace, gas, P, T, X)

Return viscosity [Pa s] and conductivity [W/(m K)], and fill
`workspace.diffusion` with mole-gradient mixture diffusion coefficients [m²/s].
`X` must contain normalized, nonnegative mole fractions. Requires transport
polynomials from `mechanism/export_sidecar.py`.
"""
function mixture_transport!(w::TransportWorkspace, gas::Solution, P, T, X; basis=:mole)
    gas.trans.model == :ionized_gas && throw(ArgumentError("use ionized_transport! for ionized-gas transport"))
    n = gas.n_species
    basis in (:mole,:mass) || throw(ArgumentError("diffusion gradient basis must be :mole or :mass"))
    length(X) == n || throw(DimensionMismatch("one mole fraction per species required"))
    isfinite(P) && P > 0 && isfinite(T) && T > 0 ||
        throw(ArgumentError("temperature and pressure must be finite and positive"))
    gas.trans.poly_order == 5 || throw(ArgumentError(
        "native transport fits missing; regenerate the mechanism sidecar"))
    length(w.viscosity) == length(w.conductivity) == length(w.diffusion) == n &&
        length(w.molecular_weights) == n && size(w.binary) == size(w.weight_quarter_ratio) ==
        size(w.inverse_phi_denominator) == (n,n) || throw(DimensionMismatch("transport workspace size mismatch"))
    size(gas.trans.species_viscosities_poly) == size(gas.trans.thermal_conductivity_poly) == (5,n) &&
        size(gas.trans.binary_diff_coeffs_poly) == (5,n*n) ||
        throw(DimensionMismatch("invalid transport polynomial dimensions"))
    all(isfinite,X) && all(>=(0),X) && isapprox(sum(X),1;atol=1e-10,rtol=1e-10) ||
        throw(ArgumentError("normalized nonnegative mole fractions required"))
    if w.molecular_weights != gas.MW
        for j in 1:n,k in 1:n
            w.weight_quarter_ratio[j,k] = sqrt(sqrt(gas.MW[j]/gas.MW[k]))
            w.inverse_phi_denominator[j,k] = inv(sqrt(8*(1+gas.MW[k]/gas.MW[j])))
        end
        copyto!(w.molecular_weights,gas.MW)
    end
    logT = log(T)
    sqrtT = sqrt(T)
    root,inverse_root,mass_fraction_weight = w.root_viscosity,w.inverse_root_viscosity,w.weighted_mole_fraction
    @inbounds for k in 1:n
        w.viscosity[k] = sqrtT * _transport_polynomial(gas.trans.species_viscosities_poly, k, logT)^2
        w.conductivity[k] = sqrtT * _transport_polynomial(gas.trans.thermal_conductivity_poly, k, logT)
        root[k] = sqrt(w.viscosity[k]); inverse_root[k] = inv(root[k])
        mass_fraction_weight[k] = X[k]*gas.MW[k]
    end
    # Binary diffusion fits are symmetric: evaluate each pair once and keep the
    # reciprocal coefficient, read column-wise below.
    scale = P/(T*sqrtT)
    @inbounds for j in 1:n, i in 1:j
        value = _transport_polynomial(gas.trans.binary_diff_coeffs_poly, i + (j-1)*n, logT)/scale
        w.binary[i,j] = value; w.binary[j,i] = value
        inverse = inv(value)
        w.inverse_binary[i,j] = inverse; w.inverse_binary[j,i] = inverse
    end
    viscosity = 0.0
    arithmetic = 0.0
    harmonic = 0.0
    mean_MW = dot(X, gas.MW)
    total_mass = sum(mass_fraction_weight)
    @inbounds for k in 1:n
        denominator = 0.0
        invdiff = 0.0
        mass_invdiff = 0.0
        rk = root[k]
        for j in 1:n
            phi = muladd(rk*inverse_root[j], w.weight_quarter_ratio[j,k], 1)^2 * w.inverse_phi_denominator[j,k]
            denominator = muladd(X[j], phi, denominator)
            invdiff = muladd(X[j], w.inverse_binary[j,k], invdiff)
            mass_invdiff = muladd(mass_fraction_weight[j], w.inverse_binary[j,k], mass_invdiff)
        end
        # Remove the diagonal pair term accumulated above.
        invdiff -= X[k]*w.inverse_binary[k,k]
        mass_invdiff -= mass_fraction_weight[k]*w.inverse_binary[k,k]
        othermass = total_mass - mass_fraction_weight[k]
        viscosity += X[k] * w.viscosity[k] / denominator
        arithmetic += X[k] * w.conductivity[k]
        harmonic += X[k] / w.conductivity[k]
        w.diffusion[k] = if n == 1
            w.binary[k,k]
        elseif basis == :mass
            othermass > 0 ? inv(invdiff+X[k]*mass_invdiff/othermass) : 0.0
        else
            invdiff > 0 ? othermass / (mean_MW * invdiff) : 0.0
        end
    end
    return viscosity, (arithmetic + 1/harmonic)/2
end
export TransportWorkspace, mixture_transport!
