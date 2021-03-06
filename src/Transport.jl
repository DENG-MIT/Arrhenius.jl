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