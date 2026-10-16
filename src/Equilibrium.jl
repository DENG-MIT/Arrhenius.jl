"""
    mole_fractions(gas, composition; basis=:mole)

Return normalized mole fractions. Compositions may be vectors, species-name
dictionaries, `"H2:2, O2:1"`, or a single species name. `basis=:mass` interprets
the input amounts as masses. The input is never modified.
"""
function mole_fractions(gas::Solution, composition; basis=:mole)
    basis = Symbol(basis)
    basis in (:mole, :mass) || throw(ArgumentError("basis must be :mole or :mass"))
    if composition isa AbstractString
        entries = Dict{String,Float64}()
        for item in split(composition, ',')
            fields = strip.(split(strip(item), ':'; limit=2))
            isempty(fields[1]) && throw(ArgumentError("empty species name"))
            value = length(fields) == 1 ? 1.0 : tryparse(Float64, fields[2])
            isnothing(value) && throw(ArgumentError("invalid composition amount: $item"))
            isfinite(value) && value >= 0 || throw(ArgumentError("composition amounts must be finite and nonnegative"))
            entries[fields[1]] = get(entries, fields[1], 0.0) + value
        end
        composition = entries
    end
    if composition isa AbstractDict
        X = zeros(gas.n_species)
        for (name, value) in composition
            k = findfirst(==(String(name)), gas.species_names)
            isnothing(k) && throw(ArgumentError("unknown species: $name"))
            X[k] = value
        end
    elseif composition isa AbstractVector
        X = Float64.(composition)
    else
        throw(ArgumentError("composition must be a vector, dictionary, or string"))
    end
    length(X) == gas.n_species || throw(DimensionMismatch("one fraction per species required"))
    all(isfinite, X) && all(>=(0), X) && maximum(X) > 0 ||
        throw(ArgumentError("composition must be finite, nonnegative, and nonzero"))
    X ./= maximum(X) # avoid overflow when normalizing large finite amounts
    basis == :mass && (X ./= gas.MW)
    return X ./ sum(X)
end

"Return normalized mass fractions, interpreting input on the specified basis."
function mass_fractions(gas::Solution, composition; basis=:mole)
    X = mole_fractions(gas, composition; basis)
    return X .* gas.MW ./ dot(X, gas.MW)
end

function _oxygen_demand(gas)
    demand = zeros(gas.n_species)
    for (element, coefficient) in (("C", 1.0), ("H", 0.25), ("S", 1.0), ("O", -0.5))
        i = findfirst(==(element), gas.elements)
        isnothing(i) || (demand .+= coefficient .* gas.ele_matrix[i, :])
    end
    return demand
end

"""
    set_mixture_fraction(gas, Z; fuel, oxidizer, basis=:mole)

Return mole fractions with fuel mass fraction `Z` in the combined streams.
`basis` describes the two input stream compositions, not `Z`.
"""
function set_mixture_fraction(gas::Solution, Z; fuel, oxidizer, basis=:mole)
    isfinite(Z) && 0 <= Z <= 1 || throw(ArgumentError("mixture fraction must lie in [0, 1]"))
    fuelY = mass_fractions(gas, fuel; basis)
    oxidizerY = mass_fractions(gas, oxidizer; basis)
    return mole_fractions(gas, Z .* fuelY .+ (1-Z) .* oxidizerY; basis=:mass)
end

"""
    set_equivalence_ratio(gas, phi; fuel, oxidizer, basis=:mole,
                          diluent=nothing, fraction=nothing)

Return mole fractions at the requested equivalence ratio. Complete oxidation
forms CO2, H2O and SO2; other elements are inert in the mixture definition.
Optional dilution specifies exactly one `:diluent`, `:fuel`, or `:oxidizer`
fraction in a dictionary or named tuple, using the supplied `basis`.
"""
function set_equivalence_ratio(gas::Solution, phi; fuel, oxidizer, basis=:mole,
                               diluent=nothing, fraction=nothing)
    isfinite(phi) && phi >= 0 || throw(ArgumentError("equivalence ratio must be finite and nonnegative"))
    fuelX = mole_fractions(gas, fuel; basis)
    oxidizerX = mole_fractions(gas, oxidizer; basis)
    demand = _oxygen_demand(gas)
    df, dox = dot(demand, fuelX), dot(demand, oxidizerX)
    df > 0 && dox < 0 || throw(ArgumentError("fuel must require oxygen and oxidizer must supply oxygen"))
    # Scale before combining so very large finite phi cannot overflow.
    ratio = -dox / df
    fuelpart = phi == 0 ? 0.0 : 1 / (1 + (1/phi)/ratio)
    X = fuelpart .* fuelX .+ (1-fuelpart) .* oxidizerX
    if isnothing(diluent)
        isnothing(fraction) || throw(ArgumentError("fraction requires a diluent"))
        return X
    end
    (fraction isa AbstractDict || fraction isa NamedTuple) && length(fraction) == 1 ||
        throw(ArgumentError("dilution requires one fraction: diluent, fuel, or oxidizer"))
    key, value = first(pairs(fraction))
    key = Symbol(key)
    key in (:diluent, :fuel, :oxidizer) || throw(ArgumentError("unknown dilution fraction: $key"))
    isfinite(value) && 0 <= value <= 1 || throw(ArgumentError("dilution fraction must lie in [0, 1]"))
    diluentX = mole_fractions(gas, diluent; basis)
    if Symbol(basis) == :mass
        fuelpart = fuelpart * dot(fuelX, gas.MW) / dot(X, gas.MW)
        X = mass_fractions(gas, X)
        diluentX = mass_fractions(gas, diluentX)
    end
    if key == :diluent
        retained = 1-value
    else
        stream = key == :fuel ? fuelpart : 1-fuelpart
        stream > 0 || throw(ArgumentError("requested stream is absent from the undiluted mixture"))
        retained = value / stream
        retained <= 1+1e-14 || throw(ArgumentError("requested fraction cannot be obtained by dilution"))
        retained = min(retained,1.0)
    end
    return mole_fractions(gas, retained .* X .+ (1-retained) .* diluentX; basis)
end

"""
    mixture_fraction(gas, X; fuel, oxidizer, basis=:mole, element="Bilger")

Return fuel mass fraction inferred from conserved elements. `X` is a mole
composition; `basis` applies only to fuel and oxidizer. The default Bilger
definition combines C, H, S and O; an individual element may also be selected.
"""
function mixture_fraction(gas::Solution, X; fuel, oxidizer, basis=:mole, element="Bilger")
    mole = mole_fractions(gas, X)
    fuelX = mole_fractions(gas, fuel; basis)
    oxidizerX = mole_fractions(gas, oxidizer; basis)
    if String(element) == "Bilger"
        weights = _oxygen_demand(gas)
    else
        i = findfirst(==(String(element)), gas.elements)
        isnothing(i) && throw(ArgumentError("unknown element: $element"))
        weights = gas.ele_matrix[i, :]
    end
    beta(x) = dot(weights, x) / dot(gas.MW, x)
    bf, bo = beta(fuelX), beta(oxidizerX)
    abs(bf-bo) > eps(Float64)*max(abs(bf),abs(bo)) ||
        throw(ArgumentError("fuel and oxidizer have indistinguishable element fractions"))
    return clamp((beta(mole)-bo)/(bf-bo), 0.0, 1.0)
end

"""
    equivalence_ratio(gas, X; fuel=nothing, oxidizer=nothing, basis=:mole,
                      include_species=nothing)

Return equivalence ratio from a mole composition. With explicit streams, use
the conserved Bilger mixture fraction. Without streams, assume all C/H/S
originate in fuel and all O in oxidizer. Optionally restrict to named species.
"""
function equivalence_ratio(gas::Solution, X; fuel=nothing, oxidizer=nothing,
                           basis=:mole, include_species=nothing)
    mole = mole_fractions(gas, X)
    if !isnothing(include_species)
        selected = zeros(gas.n_species)
        for name in include_species
            k = findfirst(==(String(name)), gas.species_names)
            isnothing(k) && throw(ArgumentError("unknown species: $name"))
            selected[k] = mole[k]
        end
        mole = mole_fractions(gas, selected)
    end
    if isnothing(fuel) && isnothing(oxidizer)
        oxygen = findfirst(==("O"), gas.elements)
        supplied = isnothing(oxygen) ? 0.0 : dot(gas.ele_matrix[oxygen, :], mole)/2
        needed = dot(_oxygen_demand(gas), mole) + supplied
        return supplied == 0 ? (needed == 0 ? NaN : Inf) : needed/supplied
    end
    !isnothing(fuel) && !isnothing(oxidizer) || throw(ArgumentError("provide both fuel and oxidizer"))
    stoich = set_equivalence_ratio(gas, 1.0; fuel, oxidizer, basis)
    Z = mixture_fraction(gas, mole; fuel, oxidizer, basis)
    Zst = mixture_fraction(gas, stoich; fuel, oxidizer, basis)
    return Z == 1 ? Inf : Z/(1-Z) * (1-Zst)/Zst
end

function _equilibrium_system(gas, X)
    any(<(0), gas.ele_matrix) && return _signed_equilibrium_system(gas,X)
    b_all = gas.ele_matrix * X
    elements = findall(>(0), b_all)
    isempty(elements) && throw(ArgumentError("composition must contain a conserved element"))
    absent = findall(==(0), b_all)
    species = [k for k in 1:gas.n_species if all(gas.ele_matrix[e,k] == 0 for e in absent)]
    A = Float64.(gas.ele_matrix[elements, species])
    all(>(0), vec(sum(A; dims=1))) || throw(ArgumentError("species without conserved elements are unsupported"))
    # Select independent original rows, preserving nonnegative atom counts for
    # the logarithmic balance equations. An SVD basis could introduce negatives.
    factor = qr(transpose(A), ColumnNorm())
    diagonal = abs.(diag(factor.R))
    rankA = count(>(maximum(diagonal)*max(size(A)...)*eps(Float64)), diagonal)
    independent = factor.p[1:rankA]
    A = A[independent, :]
    ne,ns = size(A)
    return (gas=gas, X=X, species=species, A=A, At=Matrix(transpose(A)),
            b=Float64.(b_all[elements[independent]]), state=zeros(rankA+1),
            g=zeros(ns), initialized=Ref(false), logb=log.(b_all[elements[independent]]),
            mole=zeros(ns), potentials=zeros(ns), abar=zeros(ne), residual=zeros(ne+1),
            jacobian=zeros(ne+1,ne+1), factor=zeros(ne+1,ne+1),
            trial=zeros(ne+1),step=zeros(ne+1))
end

function _equilibrium_evaluate(system, state, logpressure, constant_volume; jacobian=false)
    hasproperty(system,:signed_balance) &&
        return _signed_equilibrium_evaluate(system,state,logpressure,constant_volume;jacobian)
    A, At, g, b = system.A, system.At, system.g, system.b
    ne = size(A,1)
    v,mole,abar,f = system.potentials,system.mole,system.abar,system.residual
    mul!(v,At,view(state,1:ne))
    offset = logpressure + (constant_volume ? state[end] : 0.)
    @. v = v-g-offset
    vmax = maximum(v)
    @. mole = exp(v-vmax)
    total = sum(mole)
    logsum = vmax + log(total)
    mole ./= total
    mul!(abar,A,mole)
    @inbounds for i in 1:ne
        f[i] = log(abar[i])+state[end]-system.logb[i]
    end
    f[end] = logsum
    jacobian || return f, mole
    J = system.jacobian
    @inbounds for j in 1:ne, i in 1:ne
        moment = 0.
        for k in eachindex(mole)
            moment += A[i,k]*A[j,k]*mole[k]
        end
        J[i,j] = moment/abar[i]-abar[j]
    end
    @inbounds for i in 1:ne
        J[i,end] = 1
        J[end,i] = abar[i]
    end
    J[end,end] = constant_volume ? -1 : 0
    return f, mole, J
end

function _equilibrium_step!(system,f,J)
    # Most element systems are small and nonsingular. Use a pivoted solve;
    # retain the rank-truncated SVD for cold, nearly degenerate equilibria.
    factor = lu!(copyto!(system.factor,J);check=false)
    largest,smallest = 0.,Inf
    @inbounds for i in axes(J,1)
        pivot = abs(factor.factors[i,i])
        largest,smallest = max(largest,pivot),min(smallest,pivot)
    end
    step = system.step
    if issuccess(factor) && smallest > 1e-12*largest
        @. step = -f
        ldiv!(factor,step)
    else
        mul!(step,pinv(J;rtol=1e-14),f,-1.,0.)
    end
    return step
end

function _equilibrium_newton!(system, logpressure, constant_volume)
    state = system.state
    for iteration in 1:120
        f, _, J = _equilibrium_evaluate(system,state,logpressure,constant_volume; jacobian=true)
        all(isfinite,f) && all(isfinite,J) || return false
        residual_max,residual_norm = norm(f,Inf),norm(f)
        residual_max < 5e-14 && return true
        step = _equilibrium_step!(system,f,J)
        all(isfinite,step) || return false
        alpha = min(1.0,20/max(norm(step,Inf),1e-30))
        accepted = false
        for backtrack in 1:35
            trial = system.trial
            @. trial = state + alpha*step
            ft, _ = _equilibrium_evaluate(system,trial,logpressure,constant_volume)
            if all(isfinite,ft) && norm(ft) < residual_norm
                state .= trial
                accepted = true
                break
            end
            alpha /= 2
        end
        accepted || return residual_max < 1e-9
    end
    return false
end

function _equilibrium_at!(system, T, P; constant_volume=false, initial_temperature=T)
    gas, X, species = system.gas, system.X, system.species
    logpressure = log(P/one_atm) + (constant_volume ? log(T/initial_temperature) : 0.0)
    gtarget = Float64.((cal_h_RT(gas,T,P,X)-cal_s0_R(gas,T,P,X))[species])
    system.g .= gtarget
    converged = system.initialized[] && _equilibrium_newton!(system,logpressure,constant_volume)
    if !converged && !system.initialized[]
        # Previously equilibrated compositions provide element-potential
        # information directly. Use resolved species only; fresh reactants or
        # rank-deficient cold products retain the homotopy initialization.
        resolved = findall(k -> X[species[k]] > 1e-20,eachindex(species))
        ne = size(system.A,1)
        if length(resolved) > ne
            B = system.At[resolved,:]
            factor = qr(B,ColumnNorm())
            diagonal = abs.(diag(factor.R))
            if minimum(diagonal) > 1e-12*maximum(diagonal)
                rhs = gtarget[resolved] .+ log.(X[species[resolved]]) .+ logpressure
                seed = factor\rhs
                if norm(B*seed-rhs,Inf) < 5.
                    system.state[1:ne] .= seed
                    system.state[end] = 0.
                    converged = _equilibrium_newton!(system,logpressure,constant_volume)
                end
            end
        end
    end
    if !converged
        # Homotopy in standard chemical potentials avoids cold-start collapse.
        # Clamp the starting temperature to the shared NASA validity interval.
        Tmin = maximum(gas.thermo.Trange[species,1])
        Tmax = minimum(gas.thermo.Trange[species,3])
        startT = Tmin <= Tmax ? clamp(3500.0,Tmin,Tmax) : T
        gstart = Float64.((cal_h_RT(gas,startT,P,X)-cal_s0_R(gas,startT,P,X))[species])
        xseed = max.(X[species],1e-8)
        system.state[1:end-1] .= system.At \ (gstart + log.(xseed) .+ logpressure)
        system.state[end] = 0
        seed = copy(system.state)
        # A larger chemical-potential increment avoids redundant Newton solves.
        # Retain the original smaller increments as a convergence fallback.
        for increment in (96.,12.,3.)
            system.state .= seed
            stages = max(2,ceil(Int,maximum(abs.(gtarget-gstart))/increment)+1)
            converged = true
            for fraction in range(0,1; length=stages)
                system.g .= (1-fraction).*gstart .+ fraction.*gtarget
                if !_equilibrium_newton!(system,logpressure,constant_volume)
                    converged = false
                    break
                end
            end
            converged && break
        end
        converged || error("ideal-gas equilibrium did not converge at $T K and $P Pa")
    end
    system.initialized[] = true
    _, mole = _equilibrium_evaluate(system,system.state,logpressure,constant_volume)
    result = zeros(gas.n_species)
    result[species] = mole
    # Check every element, including dependent rows removed from Newton.
    before = gas.ele_matrix*X ./ dot(gas.MW,X)
    after = gas.ele_matrix*result ./ dot(gas.MW,result)
    norm(after-before,Inf) <= 2e-9*max(norm(before,Inf),1e-30) ||
        error("equilibrium failed elemental conservation")
    finalP = constant_volume ? P*T/initial_temperature*dot(gas.MW,X)/dot(gas.MW,result) : P
    return (T=Float64(T), P=Float64(finalP), X=result,
            Y=result .* gas.MW ./ dot(gas.MW,result))
end

"""
    equilibrate(gas; T, P=one_atm, X, mode=:TP, temperature_bounds=(200,6000),
                property_rtol=1e-10)

Native single-phase ideal-gas equilibrium. Supported conserved pairs are `:TP`,
`:TV`, `:HP`, `:UV`, `:SP` and `:SV` (symbols or strings). The supplied state
defines the conserved specific enthalpy, internal energy, entropy, or volume.
Returns `(T, P, X, Y)`; no input is modified. Absent ordinary elements stay
absent; signed electron counts conserve the supplied net charge.

The temperature search uses the explicit `temperature_bounds`. Like the
underlying NASA property functions, evaluation outside species fit ranges
extrapolates the polynomials. Supply bounds within the shared fit interval when
extrapolation is unsuitable. This solver does not support nonideal, surface,
or multiphase equilibrium.
`property_rtol` controls the outer conserved-property temperature solve for HP,
UV, SP and SV states, with an enthalpy/energy scale floor of 1e6 J/kg or an
entropy scale floor of 1e3 J/(kg K). Element-balance tolerances are unchanged.
"""
function equilibrate(gas::Solution; T, P=one_atm, X, mode=:TP, temperature_bounds=(200.0,6000.0),
                     property_rtol=1e-10)
    gas.thermo isa IdealGasThermo || throw(ArgumentError("ideal-gas thermo required"))
    isfinite(T) && T > 0 && isfinite(P) && P > 0 || throw(ArgumentError("positive finite T and P required"))
    mode = Symbol(mode)
    mode in (:TP,:TV,:HP,:UV,:SP,:SV) || throw(ArgumentError("unsupported equilibrium mode: $mode"))
    isfinite(property_rtol) && 0 < property_rtol < 1 || throw(ArgumentError("property_rtol must lie between zero and one"))
    length(temperature_bounds) == 2 || throw(ArgumentError("provide lower and upper temperature bounds"))
    lo, hi = Float64.(temperature_bounds)
    isfinite(lo) && isfinite(hi) && 0 < lo < hi || throw(ArgumentError("temperature bounds must be finite, positive, and ordered"))
    initial = mole_fractions(gas,X)
    system = _equilibrium_system(gas,initial)
    return _equilibrate(system,initial,Float64(T),Float64(P),mode,lo,hi,property_rtol)
end

# A caller following the same conserved elements may reuse the element
# potentials across successive thermodynamic constraints. The target property
# still comes from the explicitly supplied initial state.
function _equilibrate(system,initial,T,P,mode,lo,hi,property_rtol)
    gas = system.gas
    constant_volume = mode in (:TV,:UV,:SV)
    state(temperature) = _equilibrium_at!(system,temperature,Float64(P); constant_volume, initial_temperature=Float64(T))
    mode in (:TP,:TV) && return state(Float64(T))
    property = mode == :HP ? cal_hmass_mean : mode == :UV ? cal_umass_mean : cal_smass_mean
    target = property(gas,T,P,initial)
    scale = mode in (:SP,:SV) ? max(abs(target),1e3) : max(abs(target),1e6)
    # Enthalpy/energy constraints can start from a hot equilibrium, where the
    # element system is well conditioned; a cold inlet TP equilibrium is not
    # needed to determine its adiabatic equilibrium temperature. Entropy
    # constraints retain the nearby input state for pressure perturbations.
    trialT = clamp(mode in (:HP,:UV) ? max(Float64(T),3500.) : Float64(T),lo,hi)
    trial = state(trialT)
    residual = property(gas,trial.T,trial.P,trial.X)-target
    abs(residual) < property_rtol*scale && return trial
    if residual < 0
        lower, fl = trialT, residual
        upper, fu = trialT, residual
        while fu < 0 && upper < hi
            upper = min(hi,max(upper+100,1.4*upper))
            trial = state(upper)
            fu = property(gas,trial.T,trial.P,trial.X)-target
        end
    else
        upper, fu = trialT, residual
        lower, fl = trialT, residual
        while fl > 0 && lower > lo
            lower = max(lo,min(lower-100,lower/1.4))
            trial = state(lower)
            fl = property(gas,trial.T,trial.P,trial.X)-target
        end
    end
    fl <= 0 <= fu || throw(ArgumentError("$mode equilibrium is not bracketed in $(lo)–$(hi) K"))
    for iteration in 1:80
        # Safeguarded secant gives fast convergence without a numerical derivative.
        width = upper-lower
        trialT = clamp((lower*fu-upper*fl)/(fu-fl),lower+0.05*width,upper-0.05*width)
        trial = state(trialT)
        residual = property(gas,trial.T,trial.P,trial.X)-target
        (abs(residual) <= property_rtol*scale || width < min(1e-7,property_rtol*trialT)) && return trial
        if residual > 0
            upper, fu = trialT, residual
        else
            lower, fl = trialT, residual
        end
    end
    error("$mode equilibrium temperature search did not converge")
end

_equilibrium_tp(gas,T,P,X) = equilibrate(gas; T,P,X,mode=:TP).X

export mole_fractions, mass_fractions, set_equivalence_ratio, equivalence_ratio,
       set_mixture_fraction, mixture_fraction, equilibrate
