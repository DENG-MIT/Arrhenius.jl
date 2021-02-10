function wdot_func(T, C, S0, h_mole)
    _kf = kf_func.(1:gas.n_reactions; T, C)
    ΔS_R = vk' * S0 / R
    ΔH_RT = vk' * h_mole / (R * T)
    Keq = exp.(ΔS_R .- ΔH_RT .+ log(one_atm / R / T) .* sum(vk, dims = 1)[1, :])
    _kr = @. _kf / Keq * reaction.is_reversible
    _qdot = qdot_func.(1:gas.n_reactions; _kf, _kr, C)
    return vk * _qdot  #, _qdot, _kf, _kr
end

function kf_func(i; T, C)
    kfi = Arrhenius_coeffs[i, 1] * T^Arrhenius_coeffs[i, 2] *
            exp(Arrhenius_coeffs[i, 3] * (4184.0 / R / T))

    if i in reaction.index_three_body
        kfi *= dot(@view(reaction.efficiencies_coeffs[:, i]), C)
    end

    if i in reaction.index_falloff

        j = findfirst(reaction.index_falloff .== i)

        A0 = reaction.Arrhenius_A0[j]
        b0 = reaction.Arrhenius_b0[j]
        Ea0 = reaction.Arrhenius_Ea0[j]
        k0 = A0 * exp(b0 * log(T) - Ea0 * 4184.0 / R / T)
        Pr = k0 * dot(@view(reaction.efficiencies_coeffs[:, i]), C) / kfi
        lPr = log10(Pr)

        kfi *= (Pr / (1 + Pr))

        if !(i in reaction.list_type4_noTroe)
            jj = j - length(findall(reaction.list_type4_noTroe .< i)) + 1
            F_cent =
                (1 - reaction.Troe_A[jj]) * exp(-T / reaction.Troe_T3[jj]) +
                reaction.Troe_A[jj] * exp(-T / reaction.Troe_T1[jj]) +
                exp(-reaction.Troe_T2[jj] / T)
            lF_cent = log10(F_cent)
            _C = -0.4 - 0.67 * lF_cent
            N = 0.75 - 1.27 * lF_cent
            f1 = (lPr + _C) / (N - 0.14 * (lPr + _C))
            kfi *= exp(log(10.0) * lF_cent / (1 + f1^2))
        end
    end
    kfi
end

function qdot_func(i; _kf, _kr, C)
    rop_f = _kf[i]
    for j in i_reactant[i]
        rop_f *= C[j]^reactant_orders[j, i]
    end
    if reaction.is_reversible[i]
        rop_r = _kr[i]
        for j in i_product[i]
            rop_r *= C[j]^product_stoich_coeffs[j, i]
        end
        rop_f -= rop_r
    end
    rop_f
end

    # _qdot = _kf * 1.0
    # _qdot = []
    # for i = 1:gas.n_reactions
    #     rop_f = _kf[i]
    #     for j in i_reactant[i]
    #         rop_f *= C[j]^reactant_orders[j, i]
    #     end
    #
    #     if reaction.is_reversible[i]
    #         rop_r = _kr[i]
    #         for j in i_product[i]
    #             rop_r *= C[j]^product_stoich_coeffs[j, i]
    #         end
    #         # _qdot[i] = rop_f - rop_r
    #         push!(_qdot, rop_f - rop_r)
    #     else
    #         # _qdot[i] = rop_f
    #         push!(_qdot, rop_f)
    #     end
    # end
