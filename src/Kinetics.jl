@inline function _plog_group_rate(plog::PlogData, group, T, logT)
    rate = zero(T)
    activation_scale = oftype(T, 4184.0 / R) / T
    @inbounds for row in plog.rate_offsets[group]:(plog.rate_offsets[group + 1] - 1)
        A = plog.Arrhenius_coeffs[row, 1]
        b = plog.Arrhenius_coeffs[row, 2]
        Ea = plog.Arrhenius_coeffs[row, 3]
        rate += A * exp(b * logT - Ea * activation_scale)
    end
    return rate
end

@inline function _plog_rate(plog::PlogData, plog_index, T, P, logT)
    first_group = plog.group_offsets[plog_index]
    last_group = plog.group_offsets[plog_index + 1] - 1

    if P <= plog.pressures[first_group]
        return _plog_group_rate(plog, first_group, T, logT)
    elseif P >= plog.pressures[last_group]
        return _plog_group_rate(plog, last_group, T, logT)
    end

    lower_group = first_group
    @inbounds while P > plog.pressures[lower_group + 1]
        lower_group += 1
    end
    upper_group = lower_group + 1
    lower_rate = _plog_group_rate(plog, lower_group, T, logT)
    upper_rate = _plog_group_rate(plog, upper_group, T, logT)
    @inbounds fraction =
        log(P / plog.pressures[lower_group]) /
        log(plog.pressures[upper_group] / plog.pressures[lower_group])
    # Cantera sums duplicate expressions at one pressure before taking the log.
    tiny = _positive_floor(lower_rate)
    return exp(log(lower_rate + tiny) +
               fraction * (log(upper_rate + tiny) - log(lower_rate + tiny)))
end

@inline _positive_floor(value::T) where {T<:AbstractFloat} = nextfloat(zero(T))
@inline _positive_floor(value) = oftype(value, 1.0e-300)

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

@inline function _plog_group_rate(
    rate_data::LogRateData,
    plog::PlogData,
    group,
    outer_temperature,
    rate_logT,
    rate_activation_scale,
)
    rate = zero(outer_temperature)
    @inbounds for row in plog.rate_offsets[group]:(plog.rate_offsets[group + 1] - 1)
        rate += _mixed_elementary_rate(
            rate_data.plog_log_a[row],
            rate_data.plog_b[row],
            rate_data.plog_ea[row],
            rate_logT,
            rate_activation_scale,
            outer_temperature,
        )
    end
    return rate
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
        return _plog_group_rate(
            rate_data,
            plog,
            first_group,
            T,
            rate_logT,
            rate_activation_scale,
        )
    elseif P >= plog.pressures[last_group]
        return _plog_group_rate(
            rate_data,
            plog,
            last_group,
            T,
            rate_logT,
            rate_activation_scale,
        )
    end
    lower_group = first_group
    @inbounds while P > plog.pressures[lower_group + 1]
        lower_group += 1
    end
    upper_group = lower_group + 1
    lower_rate = _plog_group_rate(
        rate_data,
        plog,
        lower_group,
        T,
        rate_logT,
        rate_activation_scale,
    )
    upper_rate = _plog_group_rate(
        rate_data,
        plog,
        upper_group,
        T,
        rate_logT,
        rate_activation_scale,
    )
    @inbounds fraction =
        log(P / plog.pressures[lower_group]) /
        log(plog.pressures[upper_group] / plog.pressures[lower_group])
    tiny = _positive_floor(lower_rate)
    return exp(
        log(lower_rate + tiny) +
        fraction * (log(upper_rate + tiny) - log(lower_rate + tiny)),
    )
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

# Per-grid-point cache for flame Jacobians. Only temperature-dependent factors
# are cached; collider concentrations, pressure-dependent rates and mass action
# are recomputed for every state. The owning flame workspace fixes the mechanism.
mutable struct _KineticsTemperatureCache
    temperature::Float64
    forward::Vector{Float64}
    low::Vector{Float64}
    log_fcent::Vector{Float64}
    equilibrium::Vector{Float64}
end
_KineticsTemperatureCache(reaction::Reaction) = _KineticsTemperatureCache(NaN,
    zeros(reaction.n_reactions),zeros(length(reaction.index_falloff)),
    zeros(size(reaction.Troe_,1)),zeros(reaction.n_reactions))

@inline _concentration_power(c, order) = c^order
@inline function _concentration_power(c, order::AbstractFloat)
    value = convert(promote_type(typeof(c),typeof(order)),c)
    order == one(order) && return value
    order == oftype(order,2) && return value*value
    return c^order
end

function _mass_action!(workspace, reaction::Reaction, C)
    reactants = reaction.reactant_orders
    products = reaction.product_stoich_coeffs
    ri, rv = rowvals(reactants), nonzeros(reactants)
    pi, pv = rowvals(products), nonzeros(products)
    @inbounds for i in 1:reaction.n_reactions
        forward,reverse = workspace.kf[i],workspace.kr[i]
        for j in nzrange(reactants,i)
            forward *= _concentration_power(C[ri[j]],rv[j])
        end
        if reaction.is_reversible[i]
            for j in nzrange(products,i)
                reverse *= _concentration_power(C[pi[j]],pv[j])
            end
        end
        workspace.kf[i],workspace.kr[i] = forward,reverse
        workspace.rates_of_progress[i] = forward-reverse
    end
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
            temperature_cache.log_fcent[k] = log10(
                (1-reaction.Troe_[k,1])*exp(-T/reaction.Troe_[k,4])+
                reaction.Troe_[k,1]*exp(-T/reaction.Troe_[k,2])+
                exp(-reaction.Troe_[k,3]/T))
        end
    end
    end

    if !isempty(reaction.plog.reaction_indices)
        P = isnothing(pressure) ? sum(C) * R * T : pressure
        for (plog_index, reaction_index) in enumerate(reaction.plog.reaction_indices)
            @inbounds kf[reaction_index] = if isnothing(log_rate_data)
                _plog_rate_with_collider(
                    reaction.plog,
                    plog_index,
                    T,
                    P,
                    logT,
                    C,
                )
            else
                _plog_rate_with_collider(
                    log_rate_data,
                    reaction.plog,
                    plog_index,
                    T,
                    P,
                    rate_logT,
                    rate_activation_scale,
                    C,
                )
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
            @inbounds F_cent =
                (one(T) - reaction.Troe_[k, 1]) * exp(-T / reaction.Troe_[k, 4]) +
                reaction.Troe_[k, 1] * exp(-T / reaction.Troe_[k, 2]) +
                exp(-reaction.Troe_[k, 3] / T)
                log10(F_cent)
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
            mul!(workspace.delta_s, transpose(reaction.vk), S0)
            mul!(workspace.delta_h, transpose(reaction.vk), h_mole)
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
    mul!(wdot, reaction.vk, workspace.rates_of_progress)
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
    return wdot!(similar(c),gas.reaction,T,c,s,h,workspace;get_rate_constants=true)
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
)
    multiplier_type = isnothing(rate_multipliers) ? Float64 : eltype(rate_multipliers)
    workspace_type = promote_type(
        typeof(T),
        eltype(C),
        eltype(S0),
        eltype(h_mole),
        multiplier_type,
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
    )
end
export wdot_func
