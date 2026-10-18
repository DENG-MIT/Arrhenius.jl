@inline function _plog_group_lograte(plog::PlogData, group, T, logT)
    first=plog.rate_offsets[group]
    last=plog.rate_offsets[group+1]-1
    scale=oftype(T,4184.0/R)/T
    if first==last
        A=plog.Arrhenius_coeffs[first,1]
        return log(A)+plog.Arrhenius_coeffs[first,2]*logT-
            plog.Arrhenius_coeffs[first,3]*scale
    end
    rate=_positive_floor(T)
    @inbounds for row in first:last
        A=plog.Arrhenius_coeffs[row,1]
        b=plog.Arrhenius_coeffs[row,2]
        Ea=plog.Arrhenius_coeffs[row,3]
        rate+=A*exp(b*logT-Ea*scale)
    end
    log(rate)
end

@inline function _plog_rate(plog::PlogData, plog_index, T, P, logT)
    first_group = plog.group_offsets[plog_index]
    last_group = plog.group_offsets[plog_index + 1] - 1

    if P <= plog.pressures[first_group]
        return exp(_plog_group_lograte(plog, first_group, T, logT))
    elseif P >= plog.pressures[last_group]
        return exp(_plog_group_lograte(plog, last_group, T, logT))
    end

    lower_group = first_group
    @inbounds while P > plog.pressures[lower_group + 1]
        lower_group += 1
    end
    upper_group = lower_group + 1
    lower_log_rate = _plog_group_lograte(plog, lower_group, T, logT)
    upper_log_rate = _plog_group_lograte(plog, upper_group, T, logT)
    @inbounds fraction =
        log(P / plog.pressures[lower_group]) /
        log(plog.pressures[upper_group] / plog.pressures[lower_group])
    return exp(lower_log_rate + fraction * (upper_log_rate - lower_log_rate))
end

@inline _positive_floor(value::T) where {T<:AbstractFloat} = max(oftype(value, 1.0e-300), nextfloat(zero(T)))

# Effective pressure P*eps_mix of a linear-Burke (LMR-R) reaction: listed
# colliders weigh their mole fraction by eps_i(T); all other species by one.
@inline function _plog_effective_pressure(plog::PlogData, plog_index, T, P, C)
    _plog_linear_burke(plog, plog_index) || return P
    total = zero(eltype(C))
    @inbounds for c in C
        total += c
    end
    listed = zero(total); weighted = zero(total)
    logT = log(T); scale = oftype(T, 4184.0 / R) / T
    @inbounds for row in plog.lmr_offsets[plog_index]:(plog.lmr_offsets[plog_index+1]-1)
        x = C[plog.lmr_species[row]]
        weighted += x * plog.lmr_eps[row,1] * exp(plog.lmr_eps[row,2] * logT - plog.lmr_eps[row,3] * scale)
        listed += x
    end
    return P * (weighted + total - listed) / total
end

# d log k / d log P of a PLOG rate (zero outside the tabulated pressures).
@inline function _plog_log_slope(plog::PlogData, plog_index, T, P, logT)
    first_group = plog.group_offsets[plog_index]
    last_group = plog.group_offsets[plog_index + 1] - 1
    (P <= plog.pressures[first_group] || P >= plog.pressures[last_group]) && return zero(T)
    lower_group = first_group
    @inbounds while P > plog.pressures[lower_group + 1]
        lower_group += 1
    end
    lower_log_rate = _plog_group_lograte(plog, lower_group, T, logT)
    upper_log_rate = _plog_group_lograte(plog, lower_group + 1, T, logT)
    return (upper_log_rate - lower_log_rate) /
        log(plog.pressures[lower_group + 1] / plog.pressures[lower_group])
end
@inline _positive_floor(value) = oftype(value, 1.0e-300)

@inline function _troe_fcent(A, T1, T2, T3, T)
    fcent=(one(T)-A)*exp(-T/T3)+A*exp(-T/T1)
    if T2 != 0
        fcent += exp(-T2/T)
    end
    fcent
end

"""Log-domain Arrhenius parameters for mixed-precision rate evaluation.

The outer thermochemical state and source accumulation retain their existing
precision. Only the stored Arrhenius parameters and elementary exponential use
`T`, which is typically `Float32`.
"""
struct LogRateData{T<:AbstractFloat}
    base_log_a::Vector{T}
    base_b::Vector{T}
    base_ea::Vector{T}
    low_log_a::Vector{T}
    low_b::Vector{T}
    low_ea::Vector{T}
    plog_log_a::Vector{T}
    plog_b::Vector{T}
    plog_ea::Vector{T}
end

function LogRateData(reaction::Reaction, ::Type{T}=Float32) where {T<:AbstractFloat}
    arrays = (
        reaction.Arrhenius_coeffs,
        reaction.Arrhenius_0,
        reaction.plog.Arrhenius_coeffs,
    )
    all(array -> all(>(0), @view(array[:, 1])), arrays) ||
        throw(ArgumentError("log-domain rates require positive pre-exponentials"))
    convert_column(array, column) = T.(array[:, column])
    log_column(array) = T.(log.(@view array[:, 1]))
    return LogRateData(
        log_column(arrays[1]),
        convert_column(arrays[1], 2),
        convert_column(arrays[1], 3),
        log_column(arrays[2]),
        convert_column(arrays[2], 2),
        convert_column(arrays[2], 3),
        log_column(arrays[3]),
        convert_column(arrays[3], 2),
        convert_column(arrays[3], 3),
    )
end
export LogRateData

@inline function _mixed_elementary_rate(log_a, b, ea, logT, activation_scale, outer)
    return oftype(outer, exp(log_a + b * logT - ea * activation_scale))
end

@inline function _plog_group_lograte(
    rate_data::LogRateData,
    plog::PlogData,
    group,
    outer_temperature,
    rate_logT,
    rate_activation_scale,
)
    first=plog.rate_offsets[group]
    last=plog.rate_offsets[group+1]-1
    if first==last
        return rate_data.plog_log_a[first] +
            rate_data.plog_b[first]*rate_logT -
            rate_data.plog_ea[first]*rate_activation_scale
    end
    rate=_positive_floor(outer_temperature)
    @inbounds for row in first:last
        rate+=_mixed_elementary_rate(
            rate_data.plog_log_a[row],
            rate_data.plog_b[row],
            rate_data.plog_ea[row],
            rate_logT,
            rate_activation_scale,
            outer_temperature,
        )
    end
    log(rate)
end

@inline function _plog_rate(
    rate_data::LogRateData,
    plog::PlogData,
    plog_index,
    T,
    P,
    rate_logT,
    rate_activation_scale,
)
    first_group = plog.group_offsets[plog_index]
    last_group = plog.group_offsets[plog_index + 1] - 1
    if P <= plog.pressures[first_group]
return exp(_plog_group_lograte(
            rate_data, plog, first_group, T, rate_logT, rate_activation_scale))
    elseif P >= plog.pressures[last_group]
return exp(_plog_group_lograte(
            rate_data, plog, last_group, T, rate_logT, rate_activation_scale))
    end
    lower_group = first_group
    @inbounds while P > plog.pressures[lower_group + 1]
        lower_group += 1
    end
    upper_group = lower_group + 1
    lower_log_rate = _plog_group_lograte(
        rate_data, plog, lower_group, T, rate_logT, rate_activation_scale)
    upper_log_rate = _plog_group_lograte(
        rate_data, plog, upper_group, T, rate_logT, rate_activation_scale)
    @inbounds fraction =
        log(P / plog.pressures[lower_group]) /
        log(plog.pressures[upper_group] / plog.pressures[lower_group])
    return exp(lower_log_rate + fraction * (upper_log_rate - lower_log_rate))
end

@inline function _plog_rate_with_collider(
    rate_data::LogRateData,
    plog::PlogData,
    plog_index,
    T,
    P,
    rate_logT,
    rate_activation_scale,
    C,
)
    rate = _plog_rate(
        rate_data,
        plog,
        plog_index,
        T,
        P,
        rate_logT,
        rate_activation_scale,
    )
    @inbounds collider_index = plog.collider_indices[plog_index]
    return collider_index > 0 ? rate * C[collider_index] : rate
end

@inline function _plog_rate_with_collider(
    plog::PlogData,
    plog_index,
    T,
    P,
    logT,
    C,
)
    rate = _plog_rate(plog, plog_index, T, P, logT)
    @inbounds collider_index = plog.collider_indices[plog_index]
    return collider_index > 0 ? rate * C[collider_index] : rate
end

"Reusable arrays for allocation-free kinetics evaluation."
mutable struct KineticsWorkspace{T}
    kf::Vector{T}
    kr::Vector{T}
    delta_s::Vector{T}
    delta_h::Vector{T}
    equilibrium_constants::Vector{T}
    rates_of_progress::Vector{T}
end

function KineticsWorkspace(reaction::Reaction, ::Type{T}=Float64) where {T}
    n = reaction.n_reactions
    return KineticsWorkspace{T}(
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
    )
end
export KineticsWorkspace

# Per-grid-point cache for flame Jacobians. Temperature-dependent factors and
# PLOG rates at the cached pressure are reused; collider concentrations and mass
# action are recomputed for every state. The owning flame workspace fixes the
# mechanism.
mutable struct _KineticsTemperatureCache
    temperature::Float64
    forward::Vector{Float64}
    low::Vector{Float64}
    log_fcent::Vector{Float64}
    equilibrium::Vector{Float64}
    pressure::Float64
    plog::Vector{Float64}
end
_KineticsTemperatureCache(reaction::Reaction) = _KineticsTemperatureCache(NaN,
    zeros(reaction.n_reactions),zeros(length(reaction.index_falloff)),
    zeros(size(reaction.Troe_,1)),zeros(reaction.n_reactions),
    NaN,zeros(length(reaction.plog.reaction_indices)))

@inline _concentration_power(c, order) = c^order
@inline function _concentration_power(c, order::AbstractFloat)
    value = convert(promote_type(typeof(c),typeof(order)),c)
    order == one(order) && return value
    order == oftype(order,2) && return value*value
    isinteger(order) && return value^Int(order)
    return value > 0 ? value^order : zero(value)
end

function _mass_action!(workspace, reaction::Reaction, C)
    reactants = reaction.reactant_orders
    products = reaction.product_stoich_coeffs
    ri, rv = rowvals(reactants), nonzeros(reactants)
    pi, pv = rowvals(products), nonzeros(products)
    @inbounds for i in 1:reaction.n_reactions
        forward,reverse = workspace.kf[i],workspace.kr[i]
        forward *= _mass_action_product(C, ri, rv, reaction.reactant_stoich_coeffs, i,
            nzrange(reactants, i), reaction.reactant_mass_action_terms[i])
        if reaction.is_reversible[i]
            reverse *= _mass_action_product(C, pi, pv, products, i, nzrange(products, i),
                reaction.product_mass_action_terms[i])
        end
        workspace.kf[i],workspace.kr[i] = forward,reverse
        workspace.rates_of_progress[i] = forward-reverse
    end
end

@inline function _mass_action_product(C, species, orders, stoich, column, entries,
        plan::NTuple{4,Int64})
    i1,i2,i3,n = plan
    n == 0 && return _mass_action_product(C, species, orders, stoich, column, entries)
    @inbounds c1 = C[i1]
    n == 1 && return c1
    @inbounds c2 = C[i2]
    if n == 2
        return c1 < 0 && c2 < 0 ? zero(c1) : c1*c2
    end
    @inbounds c3 = C[i3]
    negatives = Int(c1 < 0) + Int(c2 < 0) + Int(c3 < 0)
    return negatives >= 2 ? zero(c1) : (c1*c2)*c3
end

"""Cantera's optimized C1/C2/C3 and general stoichiometry concentration product."""
@inline function _mass_action_product(C, species, orders, stoich, column, entries)
    nterms = 0
    standard = length(entries) <= 3
    if standard
        @inbounds for j in entries
            coefficient = stoich[species[j], column]
            if orders[j] == coefficient && coefficient >= 0 && isinteger(coefficient)
                nterms += Int(coefficient)
            else
                standard = false
                break
            end
        end
    end
    standard &= 1 <= nterms <= 3

    multiplier = isempty(entries) ? one(eltype(C)) : one(C[species[first(entries)]])
    if standard
        negative_factors = 0
        @inbounds for j in entries
            c = C[species[j]]
            for _ in 1:Int(orders[j])
                negative_factors += c < 0
                multiplier *= c
            end
        end
        return negative_factors >= 2 ? zero(multiplier) : multiplier
    end

    # Cantera's arbitrary-order path returns zero if any participating
    # concentration is nonpositive. Species with zero order are skipped.
    @inbounds for j in entries
        order = orders[j]
        order == 0 && continue
        c = C[species[j]]
        c > 0 || return zero(multiplier)
        multiplier *= c^order
    end
    return multiplier
end

@inline function _mass_action_standard(orders, stoich, species, column, entries)
    length(entries) <= 3 || return false
    @inbounds for j in entries
        coefficient = stoich[species[j], column]
        (orders[j] == coefficient && coefficient >= 0 && isinteger(coefficient)) || return false
    end
    return true
end

"""Accumulate the derivative of Cantera's active mass-action product branch."""
function _mass_action_product_gradient!(gradient, C, species, orders, stoich,
        column, entries, scale)
    standard = _mass_action_standard(orders, stoich, species, column, entries)
    nterms = 0
    if standard
        @inbounds for j in entries
            nterms += Int(stoich[species[j], column])
        end
    end
    standard &= 1 <= nterms <= 3
    if standard
        negative_factors = 0
        @inbounds for j in entries
            c = C[species[j]]
            negative_factors += Int(c < 0) * Int(stoich[species[j], column])
        end
        negative_factors >= 2 && return gradient
    else
        @inbounds for j in entries
            orders[j] == 0 && continue
            C[species[j]] > 0 || return gradient
        end
    end

    @inbounds for j in entries
        k = species[j]
        order = standard ? Int(stoich[k, column]) : orders[j]
        order == 0 && continue
        derivative = zero(scale)
        if standard
            # Repeat the species factor `order` times, omitting each occurrence
            # once in turn; this handles zero concentrations without division.
            for _ in 1:order
                term = C[k]^(order - 1)
                for other in entries
                    other == j && continue
                    term *= C[species[other]]^Int(stoich[species[other], column])
                end
                derivative += term
            end
        else
            derivative = order * C[k]^(order - 1)
            for other in entries
                other == j && continue
                derivative *= C[species[other]]^orders[other]
            end
        end
        gradient[k] += scale * derivative
    end
    return gradient
end

# Prepared indices are valid only for this mechanism and reversibility mask.
# Reaction arrays are mutable, so check the snapshot before every internal use.
struct _ReversibleRatePlan{R,V}
    reaction::R
    flags::V
    reversible::Vector{Int}
end
function _ReversibleRatePlan(reaction::Reaction)
    flags=copy(reaction.is_reversible)
    length(flags)==reaction.n_reactions || throw(DimensionMismatch("reaction reversibility mask has the wrong length"))
    _ReversibleRatePlan(reaction,flags,findall(flags))
end
@inline _validate_reversible_plan(::Nothing,reaction)=nothing
@inline function _validate_reversible_plan(plan::_ReversibleRatePlan,reaction)
    plan.reaction === reaction || throw(ArgumentError("reverse-rate plan belongs to a different mechanism; rebuild the signed reactor RHS"))
    plan.flags == reaction.is_reversible || throw(ArgumentError("reaction reversibility changed; rebuild the signed reactor RHS"))
    nothing
end

# y = vk' * x accumulated column-wise; matches the term order of the generic
# transpose sparse matvec without its dispatch overhead.
@inline function _vk_transpose_mul!(y::Vector{T}, vk, x) where {T}
    y .= zero(T)
    rows = rowvals(vk)
    vals = nonzeros(vk)
    @inbounds for j in 1:size(vk, 2), slot in nzrange(vk, j)
        y[j] += vals[slot] * x[rows[slot]]
    end
    return y
end

# Internal factor-only entry point. A prepared plan selects reverse columns only
# for uncached ordinary rates. All cached, mixed-rate and BM calls stay full.
function _rate_factors!(reaction,T,C,S0,h_mole,workspace,reverse_plan=nothing;
        rate_multipliers=nothing,log_rate_data=nothing,temperature_cache=nothing,pressure=nothing)
    _validate_reversible_plan(reverse_plan,reaction)
    reversible_only = reverse_plan !== nothing && isnothing(temperature_cache) &&
        isnothing(log_rate_data) && isempty(reaction.blowers_masel.reaction_indices)
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
        if reversible_only
            rows,coefficients=rowvals(reaction.vk),nonzeros(reaction.vk)
            zero_s,zero_h=zero(eltype(workspace.delta_s)),zero(eltype(workspace.delta_h))
            @inbounds for i in reverse_plan.reversible
                ds,dh=zero_s,zero_h
                for j in nzrange(reaction.vk,i)
                    row=rows[j];coefficient=coefficients[j]
                    ds=muladd(coefficient,S0[row],ds)
                    dh=muladd(coefficient,h_mole[row],dh)
                end
                # Match SparseArrays' ordered accumulation and final zero addition.
                workspace.delta_s[i]=ds+zero_s
                workspace.delta_h[i]=dh+zero_h
            end
        else
            _vk_transpose_mul!(workspace.delta_s, reaction.vk, S0)
            _vk_transpose_mul!(workspace.delta_h, reaction.vk, h_mole)
        end
        equilibrium_indices = reversible_only ? reverse_plan.reversible : eachindex(kf)
        log_reference_concentration = log(one_atmosphere / gas_constant / T)
        @inbounds for i in equilibrium_indices
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

"Compute reaction source terms into preallocated storage."
function wdot!(
    wdot,
    reaction,
    T,
    C,
    S0,
    h_mole,
    workspace;
    get_qdot=false,
    rate_multipliers=nothing,
    log_rate_data=nothing,
    temperature_cache=nothing,
    pressure=nothing,
    activity_concentrations=nothing,
    get_rate_constants=false,
)
    _rate_factors!(reaction,T,C,S0,h_mole,workspace;
        rate_multipliers,log_rate_data,temperature_cache,pressure)
    kf,kr=workspace.kf,workspace.kr
    get_rate_constants && return (;forward=kf,reverse=kr,equilibrium=workspace.equilibrium_constants)
    _mass_action!(workspace,reaction,isnothing(activity_concentrations) ? C : activity_concentrations)

    if get_qdot
        return workspace.rates_of_progress
    end
    # Direct column scatter of the stoichiometric matrix; SparseArrays' generic
    # mul! dominated small-mechanism source evaluation.
    fill!(wdot, zero(eltype(wdot)))
    rows, coefficients = rowvals(reaction.vk), nonzeros(reaction.vk)
    rates = workspace.rates_of_progress
    @inbounds for i in eachindex(rates)
        q = rates[i]
        for slot in nzrange(reaction.vk, i)
            wdot[rows[slot]] = muladd(coefficients[slot], q, wdot[rows[slot]])
        end
    end
    return wdot
end
export wdot!

"Forward, reverse and concentration-equilibrium constants at an ideal-gas state."
function reaction_rate_constants(gas::Solution;T,P=one_atm,X)
    isfinite(T) && T > 0 && isfinite(P) && P > 0 || throw(ArgumentError("positive finite temperature and pressure required"))
    x = mole_fractions(gas,X)
    c = P/(R*T).*x
    h = cal_h_RT(gas,T,P,x).*(R*T)
    s = cal_s0_R(gas,T,P,x).*R
    workspace = KineticsWorkspace(gas.reaction,promote_type(typeof(T),eltype(c)))
    return wdot!(similar(c),gas.reaction,T,c,s,h,workspace;pressure=P,get_rate_constants=true)
end
export reaction_rate_constants

"compute reaction source term `dC/dt`"
function wdot_func(
    reaction,
    T,
    C,
    S0,
    h_mole;
    get_qdot=false,
    rate_multipliers=nothing,
    log_rate_data=nothing,
    pressure=nothing,
)
    multiplier_type = isnothing(rate_multipliers) ? Float64 : eltype(rate_multipliers)
    pressure_type = isnothing(pressure) ? Float64 : typeof(pressure)
    workspace_type = promote_type(
        typeof(T),
        eltype(C),
        eltype(S0),
        eltype(h_mole),
        multiplier_type,
        pressure_type,
    )
    workspace = KineticsWorkspace(reaction, workspace_type)
    wdot = zeros(workspace_type, size(reaction.vk, 1))
    return wdot!(
        wdot,
        reaction,
        T,
        C,
        S0,
        h_mole,
        workspace;
        get_qdot=get_qdot,
        rate_multipliers=rate_multipliers,
        log_rate_data=log_rate_data,
        pressure=pressure,
    )
end
export wdot_func
