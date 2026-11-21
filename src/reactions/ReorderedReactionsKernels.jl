"""
    _rate_factors!(reaction::ReorderedReactions, T, C, S0, h_mole, workspace, reverse_plan=nothing; kwargs...)

Specialized rate-factor kernel for [`ReorderedReactions`](@ref). Identical
to the generic `_rate_factors!` in every supported feature -- PLOG/
linear-Burke, Blowers-Masel, `log_rate_data`, `temperature_cache` and
`rate_multipliers` all work exactly the same way -- except the
delta_s/delta_h/equilibrium-constant computation, which walks the
contiguous `reaction.index_rev_start:reaction.n_reactions` range directly
instead of consulting `is_reversible`/a `_ReversibleRatePlan` index list
per reaction: by construction every reaction before `index_rev_start` is
plain and irreversible, so it needs no reverse-rate work at all, and
skipping it is a single range check rather than a per-reaction branch or
gather. This mirrors the LLMImprovments fork's `ReorderedKinetics`
(`vk_rev` sized to just the reversible reactions) and is the specific
piece of its `wdot_func!` design that benchmarks ~7-12% faster than the
generic evaluator here (see `validation/kinetics_scalar_scaling.md`).
`reverse_plan` is accepted for call-site compatibility with the generic
method but unused -- `index_rev_start` already carries that information,
structurally, for every call.
"""
function _rate_factors!(reaction::ReorderedReactions,T,C,S0,h_mole,workspace,reverse_plan=nothing;
        rate_multipliers=nothing,log_rate_data=nothing,temperature_cache=nothing,pressure=nothing)
    kf = workspace.kf
    kr = workspace.kr
    logT = log(T)
    gas_constant = oftype(T, R)
    one_atmosphere = oftype(T, one_atm)
    activation_scale = oftype(T, 4184.0 / R) / T
    cached = temperature_cache !== nothing && log_rate_data === nothing && T isa Float64 &&
        isempty(reaction.blowers_masel.reaction_indices)
    refresh_temperature = !cached || temperature_cache.temperature != T
    if !isnothing(log_rate_data)
        length(log_rate_data.base_log_a) == length(kf) ||
            throw(DimensionMismatch(
                "log_rate_data must contain one base rate per reaction",
            ))
        RateScalar = eltype(log_rate_data.base_log_a)
        rate_temperature = RateScalar(T)
        rate_logT = log(rate_temperature)
        rate_activation_scale = RateScalar(4184.0 / R) / rate_temperature
    end
    if !refresh_temperature
        copyto!(kf,temperature_cache.forward)
    else
    @inbounds for i in eachindex(kf)
        if isnothing(log_rate_data)
            kf[i] = reaction.Arrhenius_coeffs[i, 1] * exp(
                reaction.Arrhenius_coeffs[i, 2] * logT -
                reaction.Arrhenius_coeffs[i, 3] * activation_scale,
            )
        else
            kf[i] = _mixed_elementary_rate(
                log_rate_data.base_log_a[i],
                log_rate_data.base_b[i],
                log_rate_data.base_ea[i],
                rate_logT,
                rate_activation_scale,
                T,
            )
        end
    end
    if cached
        copyto!(temperature_cache.forward,kf)
        for j in eachindex(temperature_cache.low)
            temperature_cache.low[j] = reaction.Arrhenius_0[j,1]*exp(
                reaction.Arrhenius_0[j,2]*logT-reaction.Arrhenius_0[j,3]*activation_scale)
        end
        for k in eachindex(temperature_cache.log_fcent)
            # Cantera clamps F_cent at SmallNumber before taking its logarithm.
            temperature_cache.log_fcent[k] = log10(max(
                _troe_fcent(reaction.Troe_[k,1],reaction.Troe_[k,2],
                    reaction.Troe_[k,3],reaction.Troe_[k,4],T),1e-300))
        end
    end
    end

    if !isempty(reaction.plog.reaction_indices)
        P = isnothing(pressure) ? sum(C) * R * T : pressure
        if cached
            # PLOG rates depend only on T and P; reuse them at an unchanged state.
            reuse = !refresh_temperature && temperature_cache.pressure == P
            plog = reaction.plog
            @inbounds for (plog_index, reaction_index) in enumerate(plog.reaction_indices)
                rate = if _plog_linear_burke(plog, plog_index)
                    # Composition-dependent effective pressure: never cached.
                    _plog_rate(plog, plog_index, T, _plog_effective_pressure(plog, plog_index, T, P, C), logT)
                else
                    reuse ? temperature_cache.plog[plog_index] :
                        (temperature_cache.plog[plog_index] = _plog_rate(plog, plog_index, T, P, logT))
                end
                collider_index = plog.collider_indices[plog_index]
                kf[reaction_index] = collider_index > 0 ? rate * C[collider_index] : rate
            end
            temperature_cache.pressure = P
        else
            for (plog_index, reaction_index) in enumerate(reaction.plog.reaction_indices)
                @inbounds kf[reaction_index] = if isnothing(log_rate_data)
                    _plog_rate_with_collider(
                        reaction.plog,
                        plog_index,
                        T,
                        _plog_effective_pressure(reaction.plog, plog_index, T, P, C),
                        logT,
                        C,
                    )
                else
                    _plog_rate_with_collider(
                        log_rate_data,
                        reaction.plog,
                        plog_index,
                        T,
                        _plog_effective_pressure(reaction.plog, plog_index, T, P, C),
                        rate_logT,
                        rate_activation_scale,
                        C,
                    )
                end
            end
        end
    end

    for (j,i) in enumerate(reaction.blowers_masel.reaction_indices)
        delta_h = dot(@view(reaction.vk[:,i]),h_mole)
        p = reaction.blowers_masel.coefficients
        barrier = _blowers_masel_barrier(p[j,3],p[j,4],delta_h)
        kf[i] = p[j,1]*exp(p[j,2]*logT-barrier/(R*T))
    end

    for i in reaction.index_three_body
        @inbounds kf[i] *= dot(@view(reaction.efficiencies_coeffs[:, i]), C)
    end

    for (j, i) in enumerate(reaction.index_falloff)
        @inbounds A0 = reaction.Arrhenius_0[j, 1]
        @inbounds b0 = reaction.Arrhenius_0[j, 2]
        @inbounds Ea0 = reaction.Arrhenius_0[j, 3]
        k0 = if cached
            temperature_cache.low[j]
        elseif isnothing(log_rate_data)
            A0 * exp(b0 * logT - Ea0 * activation_scale)
        else
            _mixed_elementary_rate(
                log_rate_data.low_log_a[j],
                log_rate_data.low_b[j],
                log_rate_data.low_ea[j],
                rate_logT,
                rate_activation_scale,
                T,
            )
        end
        @inbounds collider = dot(@view(reaction.efficiencies_coeffs[:, i]), C)
        if collider <= zero(collider)
            kf[i] = zero(collider)
            continue
        end
        @inbounds Pr = k0 * collider / kf[i]
        if Pr <= zero(Pr)
            kf[i] = zero(Pr)
            continue
        end
        lPr = log10(Pr)
        @inbounds kf[i] *= Pr / (1 + Pr)

        if reaction.index_falloff_Troe[j] > 0
            k = reaction.index_falloff_Troe[j]
            lF_cent = if cached
                temperature_cache.log_fcent[k]
            else
            @inbounds F_cent = _troe_fcent(
                reaction.Troe_[k,1],reaction.Troe_[k,2],
                reaction.Troe_[k,3],reaction.Troe_[k,4],T)
            log10(max(F_cent, oftype(F_cent, 1e-300)))
            end
            C_troe = -oftype(T, 0.4) - oftype(T, 0.67) * lF_cent
            N = oftype(T, 0.75) - oftype(T, 1.27) * lF_cent
            f1 = (lPr + C_troe) /
                (N - oftype(T, 0.14) * (lPr + C_troe))
            @inbounds kf[i] *= exp(
                log(oftype(T, 10.0)) * lF_cent / (one(T) + f1^2),
            )
        end
    end

    if !isnothing(rate_multipliers)
        length(rate_multipliers) == length(kf) || throw(DimensionMismatch(
            "rate_multipliers must contain one value per reaction",
        ))
        @inbounds for i in eachindex(kf)
            kf[i] *= rate_multipliers[i]
        end
    end

    if refresh_temperature
        rows,coefficients=rowvals(reaction.vk),nonzeros(reaction.vk)
        zero_s,zero_h=zero(eltype(workspace.delta_s)),zero(eltype(workspace.delta_h))
        log_reference_concentration = log(one_atmosphere / gas_constant / T)
        # Everything before index_rev_start is plain and irreversible by
        # construction (see ReorderedReactions's bucketing): skip it in one
        # range check instead of an is_reversible/_ReversibleRatePlan lookup
        # per reaction. equilibrium_constants[i] for i < index_rev_start is
        # left stale/unset, which is harmless -- the final kr assignment
        # below only reads it where is_reversible[i] is true, never there.
        @inbounds for i in reaction.index_rev_start:reaction.n_reactions
            ds,dh=zero_s,zero_h
            for j in nzrange(reaction.vk,i)
                row=rows[j];coefficient=coefficients[j]
                ds=muladd(coefficient,S0[row],ds)
                dh=muladd(coefficient,h_mole[row],dh)
            end
            # Match SparseArrays' ordered accumulation and final zero addition.
            workspace.delta_s[i]=ds+zero_s
            workspace.delta_h[i]=dh+zero_h
            workspace.equilibrium_constants[i] = exp(
                workspace.delta_s[i] / gas_constant -
                workspace.delta_h[i] / (gas_constant * T) +
                log_reference_concentration * reaction.vk_sum[i],
            )
        end
        if cached
            copyto!(temperature_cache.equilibrium,workspace.equilibrium_constants)
            temperature_cache.temperature = T
        end
    else
        copyto!(workspace.equilibrium_constants,temperature_cache.equilibrium)
    end
    @inbounds for i in eachindex(kf)
        kr[i] = reaction.is_reversible[i] ?
            kf[i] / workspace.equilibrium_constants[i] : zero(T)
    end

    return nothing
end
