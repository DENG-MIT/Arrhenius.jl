"Variable-region species models used only by phases containing non-NASA7 data."
struct SpeciesThermoData{T<:AbstractFloat}
    models::Vector{UInt8} # 0 NASA7, 1 NASA9, 2 Shomate, 3 constant-cp
    ranges::Vector{Vector{T}}
    coefficients::Vector{Matrix{T}}
end

"""Struct for the ideal gas thermo.

nasa_low: Array with low temperature nasa coeff. for each species

nasa_high: Array with high temperature nasa coeff. for each species

Trange: Array with temperature ranges for each species

isTcommon: bool which indicates if both polynoms share same T at intersection

extra: optional NASA9/Shomate/constant-cp regions. The NASA7 coefficient rows
are zero for these species; Trange retains their minimum, first split, maximum.

"""
struct IdealGasThermo{T<:AbstractFloat} <: Thermo
    nasa_low::Array{T,2}
    nasa_high::Array{T,2}
    Trange::Array{T,2}
    isTcommon::Bool
    extra::Union{Nothing,SpeciesThermoData{T}}
end

# Preserve the original constructor and storage for the NASA7 fast path.
IdealGasThermo(low::Matrix{T}, high::Matrix{T}, ranges::Matrix{T}, common::Bool) where {T<:AbstractFloat} =
    IdealGasThermo(low,high,ranges,common,nothing)
IdealGasThermo{T}(low,high,ranges,common::Bool) where {T<:AbstractFloat} =
    IdealGasThermo{T}(low,high,ranges,common,nothing)

function _thermo_quantity(value, dimension, units)
    pressure_factors = Dict("Pa"=>1.0,"kPa"=>1e3,"MPa"=>1e6,"GPa"=>1e9,
        "bar"=>1e5,"mbar"=>1e2,"atm"=>101325.0,"torr"=>101325/760,"dyn/cm^2"=>0.1)
    energy_factors = Dict("J"=>1.0,"kJ"=>1e3,"MJ"=>1e6,"cal"=>4.184,"kcal"=>4184.0,"erg"=>1e-7)
    quantity_factors = Dict("kmol"=>1.0,"mol"=>1e-3,"mole"=>1e-3,
        "molecule"=>1/6.02214076e26,"molecules"=>1/6.02214076e26)
    if value isa Real
        number, unit = Float64(value), ""
    elseif value isa AbstractString
        matched = match(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(.*?)\s*$",value)
        isnothing(matched) && throw(ArgumentError("invalid thermodynamic quantity: $value"))
        number, unit = parse(Float64,matched[1]), String(matched[2])
    else
        throw(ArgumentError("invalid thermodynamic quantity: $value"))
    end
    if dimension == :temperature
        isempty(unit) && (unit = String(get(units,"temperature","K")))
        unit == "K" || throw(ArgumentError("thermodynamic temperatures must use kelvin"))
        return number
    elseif dimension == :pressure
        isempty(unit) && (unit = String(get(units,"pressure","Pa")))
        haskey(pressure_factors,unit) || throw(ArgumentError("unsupported thermodynamic pressure unit: $unit"))
        return number*pressure_factors[unit]
    end
    if isempty(unit)
        energy = String(get(units,"energy","J"))
        quantity = String(get(units,"quantity","kmol"))
    else
        # Common Cantera spellings: J/kmol/K, cal/mol/K, J/(mol*K).
        parts = split(replace(unit," "=>"","("=>"",")"=>"","*"=>"/"),'/')
        expected = dimension == :molar_energy ? 2 : 3
        length(parts) == expected && (expected == 2 || parts[3] == "K") ||
            throw(ArgumentError("unsupported thermodynamic molar unit: $unit"))
        energy, quantity = parts[1:2]
    end
    haskey(energy_factors,energy) && haskey(quantity_factors,quantity) ||
        throw(ArgumentError("unsupported thermodynamic energy/amount unit: $energy/$quantity"))
    return number*energy_factors[energy]/quantity_factors[quantity]
end

@inline function _extended_thermo(thermo::IdealGasThermo, i, T)
    extra = thermo.extra
    model = extra.models[i]
    if model == 0
        a = _nasa7_coefficients(thermo,i,T)
        cp = a[i,1]+T*(a[i,2]+T*(a[i,3]+T*(a[i,4]+T*a[i,5])))
        h = a[i,1]+T*(a[i,2]/2+T*(a[i,3]/3+T*(a[i,4]/4+T*a[i,5]/5)))+a[i,6]/T
        s = a[i,1]*log(T)+T*(a[i,2]+T*(a[i,3]/2+T*(a[i,4]/3+T*a[i,5]/4)))+a[i,7]
        return cp,h,s
    end
    bounds = extra.ranges[i]
    # NASA9 selects the upper region exactly at an interior breakpoint;
    # NASA7 and Shomate select the lower one (Cantera source behavior).
    region = model == 1 ? searchsortedlast(bounds,T) : searchsortedfirst(bounds,T)-1
    region = clamp(region,1,size(extra.coefficients[i],2))
    a = view(extra.coefficients[i],:,region)
    if model == 1
        invT = inv(T)
        cp = a[1]*invT^2+a[2]*invT+a[3]+T*(a[4]+T*(a[5]+T*(a[6]+T*a[7])))
        h = -a[1]*invT^2+a[2]*log(T)*invT+a[3]+T*(a[4]/2+T*(a[5]/3+T*(a[6]/4+T*a[7]/5)))+a[8]*invT
        s = -a[1]*invT^2/2-a[2]*invT+a[3]*log(T)+T*(a[4]+T*(a[5]/2+T*(a[6]/3+T*a[7]/4)))+a[9]
    elseif model == 2
        t = T/1000
        cp = a[1]+t*(a[2]+t*(a[3]+t*a[4]))+a[5]/t^2
        h = a[1]+t*(a[2]/2+t*(a[3]/3+t*a[4]/4))-a[5]/t^2+a[6]/t
        s = a[1]*log(t)+t*(a[2]+t*(a[3]/2+t*a[4]/3))-a[5]/(2t^2)+a[7]
    else # coefficients are T0, h0/R, s0/R, cp0/R
        cp = a[4]+zero(T)
        h = (a[2]+a[4]*(T-a[1]))/T
        s = a[3]+a[4]*log(T/a[1])
    end
    return cp,h,s
end

function _extended_thermo!(output,thermo,T,property)
    for i in eachindex(output)
        output[i] = _extended_thermo(thermo,i,T)[property]
    end
    return output
end

@inline function _nasa7_coefficients(thermo::IdealGasThermo, species, T)
    return T <= thermo.Trange[species, 2] ? thermo.nasa_low : thermo.nasa_high
end

function cal_h_RT!(output, gas::Solution, thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    isnothing(thermo.extra) || return _extended_thermo!(output,thermo,T,2)
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    @inbounds for i in eachindex(output)
        nasa = _nasa7_coefficients(thermo, i, T)
        output[i] = nasa[i, 1] + nasa[i, 2] * T / 2 + nasa[i, 3] * T2 / 3 +
            nasa[i, 4] * T3 / 4 + nasa[i, 5] * T4 / 5 + nasa[i, 6] / T
    end
    return output
end
cal_h_RT!(output, gas::Solution, T, p, X) =
    cal_h_RT!(output, gas, gas.thermo, T, p, X)
export cal_h_RT!

function cal_s0_R!(output, gas::Solution, thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    isnothing(thermo.extra) || return _extended_thermo!(output,thermo,T,3)
    logT = log(T)
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    @inbounds for i in eachindex(output)
        nasa = _nasa7_coefficients(thermo, i, T)
        output[i] = nasa[i, 1] * logT + nasa[i, 2] * T + nasa[i, 3] * T2 / 2 +
            nasa[i, 4] * T3 / 3 + nasa[i, 5] * T4 / 4 + nasa[i, 7]
    end
    return output
end
cal_s0_R!(output, gas::Solution, T, p, X) =
    cal_s0_R!(output, gas, gas.thermo, T, p, X)
export cal_s0_R!

function cal_cp_R!(output, gas::Solution, thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    isnothing(thermo.extra) || return _extended_thermo!(output,thermo,T,1)
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    @inbounds for i in eachindex(output)
        nasa = _nasa7_coefficients(thermo, i, T)
        output[i] = nasa[i, 1] + nasa[i, 2] * T + nasa[i, 3] * T2 +
            nasa[i, 4] * T3 + nasa[i, 5] * T4
    end
    return output
end
cal_cp_R!(output, gas::Solution, T, p, X) =
    cal_cp_R!(output, gas, gas.thermo, T, p, X)
export cal_cp_R!
"""
Constructor for the idealGasThermo:

yaml:: Dict of the input yaml file
"""
function IdealGasThermo(yaml::AbstractDict)
    species_names = yaml["phases"][1]["species"]
    all(name -> name isa AbstractString,species_names) ||
        throw(ArgumentError("species thermo requires an explicit list of species names"))
    n_species = length(species_names)
    n_species > 0 || throw(ArgumentError("at least one species is required"))

    nasa_low = zeros(n_species, 7)
    nasa_high = zeros(n_species, 7)
    Trange = zeros(n_species, 3)
    models = zeros(UInt8,n_species)
    extended_ranges = [Float64[] for _ in 1:n_species]
    extended_coefficients = [zeros(0,0) for _ in 1:n_species]

    _species_names =
        [yaml["species"][i]["name"] for i = 1:length(yaml["species"])]

    for (i, species) in enumerate(species_names)
        spec = yaml["species"][findfirst(x -> x == species, _species_names)]
        thermo = spec["thermo"]
        model = get(thermo,"model","NASA7")
        model in ("NASA7","NASA9","Shomate","constant-cp") ||
            throw(ArgumentError("unsupported species thermodynamic model $model for $species"))
        units = merge(get(yaml,"units",Dict()),get(spec,"units",Dict()),get(thermo,"units",Dict()))
        pref = haskey(thermo,"reference-pressure") ?
            _thermo_quantity(thermo["reference-pressure"],:pressure,units) : one_atm
        isfinite(pref) && pref > 0 || throw(ArgumentError("reference pressure must be finite and positive"))
        # All downstream phase formulas use one_atm. Shift entropy once during
        # parsing so species with 1 bar reference data retain their true mu(T,P).
        entropy_shift = log(pref/one_atm)
        if model == "constant-cp"
            T0 = haskey(thermo,"T0") ? _thermo_quantity(thermo["T0"],:temperature,units) : 298.15
            h0 = _thermo_quantity(get(thermo,"h0",0.0),:molar_energy,units)/R
            s0 = _thermo_quantity(get(thermo,"s0",0.0),:molar_heat,units)/R+entropy_shift
            cp0 = _thermo_quantity(get(thermo,"cp0",0.0),:molar_heat,units)/R
            Tmin = _thermo_quantity(get(thermo,"T-min",0.0),:temperature,units)
            Tmax = _thermo_quantity(get(thermo,"T-max",Inf),:temperature,units)
            isfinite(T0) && T0 > 0 && all(isfinite,(h0,s0,cp0)) ||
                throw(ArgumentError("invalid constant-cp parameters for $species"))
            isfinite(Tmin) && 0 <= Tmin < Tmax || throw(ArgumentError("invalid constant-cp temperature bounds for $species"))
            models[i] = 3
            extended_ranges[i] = [Tmin,Tmax]
            extended_coefficients[i] = reshape([T0,h0,s0,cp0],4,1)
            Trange[i,:] .= (Tmin,clamp(T0,Tmin,Tmax),Tmax)
            continue
        end
        data = thermo["data"]
        ranges = [_thermo_quantity(value,:temperature,units) for value in thermo["temperature-ranges"]]
        all(isfinite,ranges) && all(>(0),ranges) && all(>(0),diff(ranges)) ||
            throw(ArgumentError("temperature ranges must be finite, positive, and increasing for $species"))
        length(ranges) == length(data)+1 && !isempty(data) ||
            throw(ArgumentError("one coefficient set per temperature interval required for $species"))
        count = model == "NASA9" ? 9 : 7
        all(row -> length(row) == count && all(isfinite,row),data) ||
            throw(ArgumentError("$model requires $count finite coefficients per region for $species"))
        model == "NASA9" || length(data) <= 2 ||
            throw(ArgumentError("$model supports only one or two regions for $species"))
        if model != "NASA7"
            coefficients = Float64.(hcat(data...))
            if model == "NASA9"
                models[i] = 1
                coefficients[9,:] .+= entropy_shift
            else
                models[i] = 2
                coefficients .*= 1000/R # Shomate is J/mol/K and kJ/mol
                coefficients[7,:] .+= entropy_shift
            end
            extended_ranges[i] = ranges
            extended_coefficients[i] = coefficients
            Trange[i,:] .= (first(ranges),ranges[2],last(ranges))
        elseif length(data) == 1 && length(ranges) == 2
            # A single NASA7 region is valid Cantera YAML. Duplicating the
            # polynomial preserves its value on both sides of an arbitrary
            # internal split while retaining the existing two-region layout.
            nasa_low[i, :] = data[1]
            nasa_high[i, :] = data[1]
            midpoint = clamp(1000.0, ranges[1], ranges[2])
            Trange[i, :] .= (ranges[1], midpoint, ranges[2])
        elseif length(data) == 2 && length(ranges) == 3
            nasa_low[i, :] = data[1]
            nasa_high[i, :] = data[2]
            Trange[i, :] .= ranges
        else
            throw(ArgumentError(
                "species $species must define one or two NASA7 regions; " *
                "found $(length(data)) data regions and $(length(ranges)) bounds",
            ))
        end
        if model == "NASA7"
            nasa_low[i,7] += entropy_shift
            nasa_high[i,7] += entropy_shift
        end

    end
    isTcommon = (maximum(Trange[:, 2]) - minimum(Trange[:, 2])) < 0.01
    extra = any(!=(0),models) ? SpeciesThermoData(models,extended_ranges,extended_coefficients) : nothing
    return IdealGasThermo(nasa_low,nasa_high,Trange,isTcommon,extra)
end

"""
    species_thermo(thermo::IdealGasThermo, T; P=one_atm)

Evaluate dimensionless ideal-gas standard-state `(cp_R, h_RT, s_R)` for every
species, in the phase's declared order. `T` is kelvin and `P` is pascals.
NASA7, multi-region NASA9, Shomate, and constant-cp models may coexist.
This standalone evaluation needs no reaction or transport sidecar. For example,
`species_thermo(IdealGasThermo(YAML.load_file("airNASA9.yaml")), 8000.0)`.
"""
function species_thermo(thermo::IdealGasThermo,T;P=one_atm)
    isfinite(T) && T > 0 && isfinite(P) && P > 0 ||
        throw(ArgumentError("positive finite temperature and pressure required"))
    scalar = promote_type(typeof(T),typeof(P),eltype(thermo.nasa_low))
    cp_R, h_RT, s_R = (Vector{scalar}(undef,size(thermo.nasa_low,1)) for _ in 1:3)
    for i in eachindex(cp_R)
        if isnothing(thermo.extra)
            a = _nasa7_coefficients(thermo,i,T)
            cp_R[i] = a[i,1]+T*(a[i,2]+T*(a[i,3]+T*(a[i,4]+T*a[i,5])))
            h_RT[i] = a[i,1]+T*(a[i,2]/2+T*(a[i,3]/3+T*(a[i,4]/4+T*a[i,5]/5)))+a[i,6]/T
            s_R[i] = a[i,1]*log(T)+T*(a[i,2]+T*(a[i,3]/2+T*(a[i,4]/3+T*a[i,5]/4)))+a[i,7]
        else
            cp_R[i],h_RT[i],s_R[i] = _extended_thermo(thermo,i,T)
        end
        s_R[i] -= log(P/one_atm)
    end
    return (;cp_R,h_RT,s_R)
end
export IdealGasThermo, species_thermo
"""
    cal_h_RT(gas, T, p, X)

calculates the dimensionless mole based enthalpy (h) for each species
"""
function cal_h_RT(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    isnothing(thermo.extra) || return _extended_thermo!(Vector{promote_type(typeof(T),eltype(thermo.nasa_low))}(undef,gas.n_species),thermo,T,2)
    H_T = [1.0, T / 2.0, T^2 / 3.0, T^3 / 4.0, T^4 / 5.0, 1.0 / T]
    if thermo.isTcommon && T <= thermo.Trange[1, 2]
        h_mole = @view(thermo.nasa_low[:, 1:6]) * H_T 
    elseif thermo.isTcommon
        h_mole = @view(thermo.nasa_high[:, 1:6]) * H_T 
    else
        h_mole = @view(thermo.nasa_high[:, 1:6]) * H_T
        use_low = T .<= @view(thermo.Trange[:, 2])
        h_mole[use_low] .= @view(thermo.nasa_low[use_low, 1:6]) * H_T
    end
    # H_mole = dot(h_mole, X)
    return h_mole
end

"""
    cal_s0_R(gas, T, p, X)

calculates the dimensionless mole based reference state entropy (s0) for each species
"""
function cal_s0_R(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    isnothing(thermo.extra) || return _extended_thermo!(Vector{promote_type(typeof(T),eltype(thermo.nasa_low))}(undef,gas.n_species),thermo,T,3)
    S_T = [log(T), T, T^2 / 2.0, T^3 / 3.0, T^4 / 4.0, 1.0]
    if thermo.isTcommon && T <= thermo.Trange[1, 2]
        S0 = @view(thermo.nasa_low[:, [1, 2, 3, 4, 5, 7]]) * S_T 
    elseif thermo.isTcommon
        S0 = @view(thermo.nasa_high[:, [1, 2, 3, 4, 5, 7]]) * S_T 
    else
        S0 = @view(thermo.nasa_high[:, [1, 2, 3, 4, 5, 7]]) * S_T
        use_low = T .<= @view(thermo.Trange[:, 2])
        S0[use_low] .=
            @view(thermo.nasa_low[use_low, [1, 2, 3, 4, 5, 7]]) * S_T
    end
    return S0 
end
export cal_s0_R

"""
    cal_s_R(gas, T, p, X)

calculates the dimensionless mole based entropy (s) for each species
"""
function cal_s_R(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
   tiny = oftype(T, 1.0e-30)
   return cal_s0_R(gas,thermo, T,p,X) - log.(max.(X, tiny)) .-
       log(p / oftype(T, one_atm))
end

"""
    cal_g_RT(gas, T, p, X)

calculates the dimensionless mole based free gibbs energy (g) for each species
"""
function cal_g_RT(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    return cal_h_RT(gas,thermo, T, p, X) - cal_s_R(gas,thermo, T, p, X)
end

"""
    cal_u_RT(gas, T, p, X)

calculates the dimensionless mole based internal energy (u) for each species
"""
function cal_u_RT(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    return cal_h_RT(gas,thermo, T, p, X) .- 1
end

"""
    cal_a_RT(gas, T, p, X)

calculates the dimensionless mole based helmholz free energy (a) for each species
"""
function cal_a_RT(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    return cal_u_RT(gas,thermo, T, p, X) - cal_s_R(gas,thermo, T, p, X)
end

"""
    cal_cp_R(gas, T, p, X)

calculates the dimensionless mole based heat capacity 
at constant pressure (cp) for each species
"""
function cal_cp_R(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
    isnothing(thermo.extra) || return _extended_thermo!(Vector{promote_type(typeof(T),eltype(thermo.nasa_low))}(undef,gas.n_species),thermo,T,1)
    cp_T = [1.0, T, T^2, T^3, T^4]
    if thermo.isTcommon && T <= thermo.Trange[1, 2]
        cp = @view(thermo.nasa_low[:, 1:5]) * cp_T
    elseif thermo.isTcommon
        cp = @view(thermo.nasa_high[:, 1:5]) * cp_T
    else
        cp = @view(thermo.nasa_high[:, 1:5]) * cp_T
        use_low = T .<= @view(thermo.Trange[:, 2])
        cp[use_low] .= @view(thermo.nasa_low[use_low, 1:5]) * cp_T
    end
    return cp
end

"""
    cal_cv_R(gas, T, p, X)

calculates the dimensionless mole based heat capacity 
at constant volume (cv) for each species
"""
function cal_cv_R(gas::Solution,thermo::IdealGasThermo, T::Real, p::Real, X::AbstractArray)
   return cal_cp_R(gas,thermo, T, p, X) .- 1
end

