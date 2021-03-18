"compute reaction source term `dC/dt`"
function wdot_func(reaction, T, C, S0, h_mole; get_qdot = false)

    @inbounds _kf = @. @view(reaction.Arrhenius_coeffs[:, 1]) * exp(
        @view(reaction.Arrhenius_coeffs[:, 2]) * log(T) -
        @view(reaction.Arrhenius_coeffs[:, 3]) * (4184.0 / R / T),
    )

    for i in reaction.index_three_body
        @inbounds _kf[i] *= dot(@view(reaction.efficiencies_coeffs[:, i]), C)
    end

    for (j, i) in enumerate(reaction.index_falloff)
        @inbounds A0 = reaction.Arrhenius_A0[j]
        @inbounds b0 = reaction.Arrhenius_b0[j]
        @inbounds Ea0 = reaction.Arrhenius_Ea0[j]
        @inbounds k0 = A0 * exp(b0 * log(T) - Ea0 * 4184.0 / R / T)
        @inbounds Pr =
            k0 * dot(@view(reaction.efficiencies_coeffs[:, i]), C) / _kf[i]
        lPr = log10(Pr)
        _kf[i] *= (Pr / (1 + Pr))

        # reference:
        # http://web.mit.edu/2.62/cantera/doc/html/classCantera_1_1Troe4.html#a38aa787421d426dfd0a587fd6fc8108e
        if reaction.index_falloff_Troe[j] > 0
            k = reaction.index_falloff_Troe[j]
            @inbounds F_cent =
                (1 - reaction.Troe_A[k]) * exp(-T / reaction.Troe_T3[k]) +
                reaction.Troe_A[k] * exp(-T / reaction.Troe_T1[k]) +
                exp(-reaction.Troe_T2[k] / T)
            lF_cent = log10(F_cent)
            _C = -0.4 - 0.67 * lF_cent
            N = 0.75 - 1.27 * lF_cent
            @inbounds f1 = (lPr + _C) / (N - 0.14 * (lPr + _C))
            @inbounds _kf[i] *= exp(log(10.0) * lF_cent / (1 + f1^2))
        end
    end

    @inbounds ΔS_R = reaction.vk' * S0 / R
    @inbounds ΔH_RT = reaction.vk' * h_mole / (R * T)
    @inbounds Keq =
        @. exp(ΔS_R - ΔH_RT + log(one_atm / R / T) * reaction.vk_sum)
    @inbounds _kr = @. _kf / Keq * reaction.is_reversible

    for i = 1:reaction.n_reactions
        @inbounds for j in reaction.i_reactant[i]
            @inbounds _kf[i] *= C[j]^reaction.reactant_orders[j, i]
        end
        if reaction.is_reversible[i]
            @inbounds for j in reaction.i_product[i]
                @inbounds _kr[i] *= C[j]^reaction.product_stoich_coeffs[j, i]
            end
        end
    end

    if get_qdot
        return _kf - _kr
    end
    return reaction.vk * (_kf - _kr)
end
export wdot_func