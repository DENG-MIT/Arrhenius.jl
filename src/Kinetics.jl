function wdot_func(gas, mgas)
    T = mgas.T

    mgas.kf = @. @view(reaction.Arrhenius_coeffs[:, 1]) * exp(
        @view(reaction.Arrhenius_coeffs[:, 2]) * log(T) -
        @view(reaction.Arrhenius_coeffs[:, 3]) * (4184.0 / R / T),
    )

    for i in reaction.index_three_body
        mgas.kf[i] *= dot(@view(reaction.efficiencies_coeffs[:, i]), mgas.C)
    end

    jj = 1
    for (j, i) in enumerate(reaction.index_falloff)
        A0 = reaction.Arrhenius_A0[j]
        b0 = reaction.Arrhenius_b0[j]
        Ea0 = reaction.Arrhenius_Ea0[j]
        k0 = A0 * exp(b0 * log(T) - Ea0 * 4184.0 / R / T)
        Pr =
            k0 * dot(@view(reaction.efficiencies_coeffs[:, i]), mgas.C) /
            mgas.kf[i]
        lPr = log10(Pr)

        mgas.kf[i] *= (Pr / (1 + Pr))

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
        mgas.kf[i] *= exp(log(10.0) * lF_cent / (1 + f1^2))
        jj += 1
    end

    ΔS_R = vk' * mgas.S0 / R
    ΔH_RT = vk' * mgas.h_mole / (R * T)
    Keq = exp.(ΔS_R .- ΔH_RT .+ log(one_atm / R / T) .* sum(vk, dims = 1)[1, :])
    mgas.kr = @. mgas.kf / Keq * reaction.is_reversible

    for i = 1:gas.n_reactions
        rop_f = mgas.kf[i]
        for j in i_reactant[i]
            rop_f *= mgas.C[j]^reactant_orders[j, i]
        end

        rop_r = mgas.kr[i]
        for j in i_product[i]
            rop_r *= mgas.C[j]^product_stoich_coeffs[j, i]
        end
        mgas.qdot[i] = rop_f - rop_r
    end
    mgas.wdot = vk * mgas.qdot
end
