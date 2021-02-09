function wdot_func(T, C, S0, h_mole)

    _kf = @. @view(reaction.Arrhenius_coeffs[:, 1]) * exp(
        @view(reaction.Arrhenius_coeffs[:, 2]) * log(T) -
        @view(reaction.Arrhenius_coeffs[:, 3]) * (4184.0 / R / T),
    )

    for i in reaction.index_three_body
        _kf[i] *= dot(@view(reaction.efficiencies_coeffs[:, i]), C)
    end

    jj = 1
    for (j, i) in enumerate(reaction.index_falloff)
        A0 = reaction.Arrhenius_A0[j]
        b0 = reaction.Arrhenius_b0[j]
        Ea0 = reaction.Arrhenius_Ea0[j]
        k0 = A0 * exp(b0 * log(T) - Ea0 * 4184.0 / R / T)
        Pr = k0 * dot(@view(reaction.efficiencies_coeffs[:, i]), C) / _kf[i]
        lPr = log10(Pr)

        _kf[i] *= (Pr / (1 + Pr))

        if i in reaction.list_type4_noTroe
            continue
        end

        F_cent =
            (1 - reaction.Troe_A[jj]) * exp(-T / reaction.Troe_T3[jj]) +
            reaction.Troe_A[jj] * exp(-T / reaction.Troe_T1[jj]) +
            exp(-reaction.Troe_T2[jj] / T)
        lF_cent = log10(F_cent)
        _C = -0.4 - 0.67 * lF_cent
        N = 0.75 - 1.27 * lF_cent
        f1 = (lPr + _C) / (N - 0.14 * (lPr + _C))
        _kf[i] *= exp(log(10.0) * lF_cent / (1 + f1^2))
        jj += 1
    end

    ΔS_R = vk' * S0 / R
    ΔH_RT = vk' * h_mole / (R * T)
    Keq = exp.(ΔS_R .- ΔH_RT .+ log(one_atm / R / T) .* sum(vk, dims = 1)[1, :])
    _kr = @. _kf / Keq * reaction.is_reversible

    _qdot = similar(_kf)
    for i = 1:gas.n_reactions
        rop_f = _kf[i]
        for j in i_reactant[i]
            rop_f *= C[j]^reactant_orders[j, i]
        end

        if reaction.is_reversible[i]
            rop_r = _kr[i]
            for j in i_product[i]
                rop_r *= C[j]^product_stoich_coeffs[j, i]
            end
            _qdot[i] = rop_f - rop_r
        else
            _qdot[i] = rop_f
        end
    end

    return vk * _qdot  #, _qdot, _kf, _kr
end
