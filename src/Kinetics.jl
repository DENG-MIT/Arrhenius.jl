function wdot_func(gas, mgas)
    T = mgas.T
    C = mgas.C
    reaction = gas.reaction

    kf = @. reaction.Arrhenius_coeffs[:, 1] * exp(reaction.Arrhenius_coeffs[:, 2] * log(T) -
                                       reaction.Arrhenius_coeffs[:, 3] * (4184.0 / R / T))

    for i in reaction.index_three_body
        kf[i] *= dot(reaction.efficiencies_coeffs[:, i], C)
    end

    j = 1
    for i in reaction.index_falloff
        kinf = kf[i]
        A0 = reaction.Arrhenius_A0[j]
        b0 = reaction.Arrhenius_b0[j]
        Ea0 = reaction.Arrhenius_Ea0[j]
        k0 = A0 * exp(b0 * log(T) - Ea0 * 4184.0 / R / T)
        Pr = k0 * dot(reaction.efficiencies_coeffs[:, i], C) / kinf
        lPr = log10(Pr)

        kf[i] *= (Pr / (1 + Pr))

        if i in reaction.list_type4_noTroe
            continue
        end

        A = reaction.Troe_A[j]
        T1 = reaction.Troe_T1[j]
        T2 = reaction.Troe_T2[j]
        T3 = reaction.Troe_T3[j]
        j = j + 1

        F_cent = (1 - A) * exp(-T / T3) + A * exp(-T / T1) + exp(-T2 / T)
        lF_cent = log10(F_cent)
        _C = -0.4 - 0.67 * lF_cent
        N = 0.75 - 1.27 * lF_cent
        f1 = (lPr + _C) / (N - 0.14 * (lPr + _C))
        kf[i] *= exp(log(10.0) * lF_cent / (1 + f1^2))
    end

    vk = reaction.product_stoich_coeffs - reaction.reactant_stoich_coeffs
    ΔS_R = vk' * mgas.S0 / R
    ΔH_RT = vk' * mgas.h_mole / (R * T)
    Keq = exp.(ΔS_R .- ΔH_RT .+ log(one_atm / R / T) .* sum(vk, dims=1)[1, :])

    kr = @. kf / Keq * gas.reaction.is_reversible

    for i in 1:gas.n_reactions
        i_reactant = findall(reaction.reactant_orders[:, i] .> 0.01)
        i_product = findall(reaction.product_stoich_coeffs[:, i] .> 0.01)
        rop_f = kf[i]
        for j in i_reactant
            rop_f *= C[j]^reaction.reactant_orders[j, i]
        end
        rop_r = kr[i]
        for j in i_product
            rop_r *= C[j]^reaction.product_stoich_coeffs[j, i]
        end
        mgas.qdot[i] = rop_f - rop_r
    end

    mgas.kf = kf
    mgas.kr = kr
    mgas.wdot = (reaction.product_stoich_coeffs - reaction.reactant_orders) * mgas.qdot
end
