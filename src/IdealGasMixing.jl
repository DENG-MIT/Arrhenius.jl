"""
    mix_constant_pressure(gas::Solution, streams)

Steady adiabatic mixing of ideal-gas streams at a common pressure. `streams`
is a tuple or vector of named tuples, each with `T` (K), `P` (Pa), `X`
(composition) and `moles` (kmol). Every stream shares `gas` — one species
ordering and ideal-gas thermodynamics — and one pressure, matched to a
relative tolerance of 1e-10. Compositions accept any form `mole_fractions`
supports. Stream amounts must be finite and nonnegative with a positive
total; streams carrying zero amount do not participate in the mixture.
Inputs are never modified.

The mixed composition and totals follow from species and mass conservation.
The mixed temperature conserves the total enthalpy,
`sum(moles * cal_h_mean(gas, T, P, X))` in J, with molar enthalpies in
J/kmol from the species thermodynamic models; their reference states fix the
energy convention. A safeguarded Newton solve with derivative `cal_cp_mean`
brackets the root between the inlet temperature extremes; equal inlet
temperatures return that exact temperature. Temperatures and amounts are
Float64, matching the `IdealGasStates` APIs; the root solve itself does not
support generic automatic differentiation.

Returns `(T, P, X, Y, moles, mass, enthalpy)`: K, Pa, mole and mass
fractions, total kmol, total mass in kg and the conserved total enthalpy in
J. The state stays chemically frozen — call `equilibrate` separately to
equilibrate it.
"""
function mix_constant_pressure(gas::Solution, streams)
    gas.thermo isa IdealGasThermo || throw(ArgumentError("ideal-gas thermo required"))
    entries = collect(streams)
    isempty(entries) && throw(ArgumentError("at least one stream is required"))
    temperatures = Float64[]
    amounts = Float64[]
    fractions = Vector{Float64}[]
    pressure = 0.0
    for stream in entries
        stream isa NamedTuple ||
            throw(ArgumentError("each stream must be a named tuple with T, P, X and moles"))
        all(field -> haskey(stream, field), (:T, :P, :X, :moles)) ||
            throw(ArgumentError("each stream requires T, P, X and moles"))
        (stream.T isa Real && stream.P isa Real && stream.moles isa Real) ||
            throw(ArgumentError("stream T, P and moles must be real numbers"))
        T, P, moles = Float64(stream.T), Float64(stream.P), Float64(stream.moles)
        (isfinite(T) && T > 0) ||
            throw(ArgumentError("stream temperatures must be finite and positive"))
        (isfinite(P) && P > 0) ||
            throw(ArgumentError("stream pressures must be finite and positive"))
        (isfinite(moles) && moles >= 0) ||
            throw(ArgumentError("stream amounts must be finite and nonnegative"))
        if iszero(pressure)
            pressure = P
        else
            abs(P - pressure) <= 1e-10 * pressure ||
                throw(ArgumentError("all streams must share a common pressure"))
        end
        push!(temperatures, T)
        push!(amounts, moles)
        push!(fractions, mole_fractions(gas, stream.X))
    end
    active = findall(>(0), amounts)
    isempty(active) && throw(ArgumentError("the streams must carry a positive total amount"))
    Ts = temperatures[active]
    ns = amounts[active]
    Xs = fractions[active]
    total = sum(ns)
    isfinite(total) || throw(ArgumentError("total stream amount must be finite"))
    X = zeros(gas.n_species)
    mass = 0.0
    enthalpy = 0.0
    for (n, x, T) in zip(ns, Xs, Ts)
        X .+= n .* x
        mass += n * dot(gas.MW, x)
        enthalpy += n * cal_h_mean(gas, T, pressure, x)
    end
    X ./= total
    mw = dot(gas.MW, X)
    isfinite(mass) && isfinite(enthalpy) && isfinite(mw) && mw > 0 ||
        throw(ArgumentError("mixture mass and enthalpy must be finite"))
    if all(==(Ts[1]), Ts)
        temperature = Ts[1]
    else
        temperature = _mixing_temperature(gas, Ts, ns, total, pressure, X, enthalpy)
    end
    cp = cal_cp_mean(gas, temperature, pressure, X)
    isfinite(cp) && cp > 0 || throw(ArgumentError("positive finite mixture heat capacity required"))
    return (T=temperature, P=pressure, X=X, Y=X .* gas.MW ./ mw,
            moles=total, mass=mass, enthalpy=enthalpy)
end

# Safeguarded Newton solve of total * h_mean(T) = enthalpy, bracketed by the
# inlet temperature extremes with a bisection fallback. The residual is
# monotone while cp stays positive, so the bracket always contains the root.
function _mixing_temperature(gas, Ts, ns, total, P, X, enthalpy)
    target = enthalpy / total
    residual(t) = cal_h_mean(gas, t, P, X) - target
    scale = max(abs(target), R * minimum(Ts), one(target))
    tolerance = 1e-12 * scale
    lower, upper = minimum(Ts), maximum(Ts)
    flow, fup = residual(lower), residual(upper)
    abs(flow) <= tolerance && return lower
    abs(fup) <= tolerance && return upper
    (flow < 0 && fup > 0) ||
        throw(ArgumentError("mixture enthalpy is not bracketed by the inlet temperatures"))
    temperature = clamp(dot(ns, Ts) / total, lower, upper)
    for _ in 1:50
        f = residual(temperature)
        abs(f) <= tolerance && return temperature
        f > 0 ? (upper = temperature) : (lower = temperature)
        cp = cal_cp_mean(gas, temperature, P, X)
        (isfinite(cp) && cp > 0) ||
            throw(ArgumentError("mixture heat capacity must stay positive during the temperature solve"))
        next = temperature - f / cp
        temperature = lower < next < upper ? next : (lower + upper) / 2
    end
    error("constant-pressure mixing temperature solve did not converge")
end

export mix_constant_pressure
