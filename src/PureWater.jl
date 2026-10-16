# Native Julia adaptation of Cantera's TPX water model (Water.cpp and the
# saturation/density algorithms in Sub.cpp), source revision
# 726522be4e2a13454d8415b7ef799d621f665cf3. The equation of state is the
# Reynolds (1979) water model used by cantera.Water(backend="Reynolds").
# https://github.com/Cantera/cantera/tree/726522be4e2a13454d8415b7ef799d621f665cf3/src/tpx
#
# Copyright (c) 2001-2009, California Institute of Technology
# All rights reserved.
# Copyright (c) 2009 Sandia Corporation. Under the terms of Contract
# AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain
# rights in this software.
# Copyright (c) 2011-2026, Cantera Developers. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# - Neither the name of the California Institute of Technology, Sandia
#   Corporation nor the names of other contributors may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

const WATER_TMIN = 273.16
const WATER_TMAX = 1600.0
const WATER_TC = 647.286
const WATER_PC = 22.089e6
const WATER_RHOC = 317.0
const WATER_MW = 18.016
const _WATER_R = 461.51
const _WATER_A = [
    2.9492937e-2 -5.1985860e-3 6.8335354e-3 -1.5641040e-4 -6.3972405e-3 -3.9661401e-3 -6.9048554e-4;
    -1.3213917e-4 7.7779182e-6 -2.6149751e-5 -7.2546108e-7 2.6409282e-5 1.5453061e-5 2.7407416e-6;
    2.7464632e-7 -3.3301902e-8 6.5326396e-8 -9.2734289e-9 -4.7740374e-8 -2.9142470e-8 -5.1028070e-9;
    -3.6093828e-10 -1.6254622e-11 -2.6181978e-11 4.3125840e-12 5.6323130e-11 2.9568796e-11 3.9636085e-12;
    3.4218431e-13 -1.7731074e-13 0 0 0 0 0;
    -2.4450042e-16 1.2748742e-16 0 0 0 0 0;
    1.5518535e-19 1.3746153e-19 0 0 0 0 0;
    5.9728487e-24 1.5597836e-22 0 0 0 0 0;
    -4.1030848e-1 3.3731180e-1 -1.3746678e-1 6.7874983e-3 1.3687317e-1 7.984797e-2 1.3041253e-2;
    -4.1605860e-4 -2.0988866e-4 -7.3396848e-4 1.0401717e-5 6.4581880e-4 3.9917570e-4 7.1531353e-5
]
const _WATER_F = (-7.4192420,.29721,-.1155286,.008685635,.001094098,-.00439993,.002520658,-.0005218684)
const _WATER_D = (3.6711257,-28.512396,222.65240,-882.43852,2000.2765,-2612.2557,1829.7674,-533.50520)
const _WATER_G = (4.6e4,1011.249,.83893,-2.19989e-4,2.466619e-7,-9.704700e-11)

"""
    PureWater()

Reynolds liquid/vapor water equation of state, evaluated entirely in Julia.
Temperatures span 273.16–1600 K; saturation spans 273.16–647.286 K. This model
does not describe ice. Thermochemical reference offsets match Cantera's water
phase. A model holds a saturation cache; use one model per concurrent task.
"""
mutable struct PureWater
    energy_offset::Float64
    entropy_offset::Float64
    sat_T::Float64
    sat_P::Float64
    liquid_density::Float64
    vapor_density::Float64
end

"Water state in SI units: T [K], P [Pa], rho [kg/m³], h/u [J/kg], s [J/kg/K]."
struct WaterState
    model::PureWater
    T::Float64
    P::Float64
    rho::Float64
    Q::Float64
    h::Float64
    u::Float64
    s::Float64
end

function Base.getproperty(state::WaterState,name::Symbol)
    name == :v && return 1/getfield(state,:rho)
    name == :density && return getfield(state,:rho)
    name == :cp && return water_cp(state)
    name == :cv && return water_cv(state)
    name == :g && return getfield(state,:h)-getfield(state,:T)*getfield(state,:s)
    name == :a && return getfield(state,:u)-getfield(state,:T)*getfield(state,:s)
    return getfield(state,name)
end
Base.propertynames(::WaterState,private::Bool=false) =
    (:T,:P,:rho,:density,:v,:Q,:h,:u,:s,:cp,:cv,:g,:a,:model)

@inline function _water_C(T,j)
    tau = 1000/T
    return j == 0 ? _WATER_R*T : _WATER_R*T*(tau-1.544912)*(tau-2.5)^(j-1)
end
@inline function _water_Cprime(T,j)
    tau = 1000/T
    j == 0 && return _WATER_R
    j == 1 && return -_WATER_R*1.544912
    return -_WATER_R*(tau-2.5)^(j-2)*(1.544912*(tau-2.5)+(j-1)*tau*(tau-1.544912))
end
@inline function _water_I(rho,j)
    factor = rho-(j == 0 ? 634.0 : 1000.0)
    total = zero(rho)
    @inbounds for i in 7:-1:1
        total = (total+_WATER_A[i+1,j+1])*factor
    end
    @inbounds total += _WATER_A[1,j+1]+exp(-.0048*rho)*(_WATER_A[9,j+1]+_WATER_A[10,j+1]*rho)
    return rho*total
end
@inline function _water_H(rho,j)
    factor = rho-(j == 0 ? 634.0 : 1000.0)
    total = zero(rho)
    @inbounds for i in 6:-1:1
        total = (total+_WATER_A[i+1,j+1]+rho*(i+1)*_WATER_A[i+2,j+1])*factor
    end
    @inbounds total += _WATER_A[1,j+1]+rho*_WATER_A[2,j+1]
    @inbounds total += exp(-.0048*rho)*((1-.0048*rho)*_WATER_A[9,j+1]+rho*(2-.0048*rho)*_WATER_A[10,j+1])
    @inbounds total += _WATER_A[8,j+1]*factor^7
    return rho*rho*total
end
function _water_pressure(T,rho)
    pressure = rho*_WATER_R*T
    for j in 0:6
        pressure += _water_C(T,j)*_water_H(rho,j)
    end
    return pressure
end
function _water_raw(T,rho)
    u = zero(T+rho)
    s = zero(T+rho)
    for j in 0:6
        integral = _water_I(rho,j)
        derivative = _water_Cprime(T,j)
        u += (_water_C(T,j)-T*derivative)*integral
        s -= derivative*integral
    end
    for i in 1:5
        u += _WATER_G[i+1]*(T^i-WATER_TMIN^i)/i
    end
    u += _WATER_G[1]*log(T/WATER_TMIN)+2375470.875
    for i in 2:5
        s += _WATER_G[i+1]*(T^(i-1)-WATER_TMIN^(i-1))/(i-1)
    end
    s += _WATER_G[2]*log(T/WATER_TMIN)-_WATER_G[1]*(1/T-1/WATER_TMIN)+6697.356635-_WATER_R*log(rho)
    return (;u,s)
end

function _water_psat_estimate(T)
    total = sum(_WATER_F[i]*(.01*(T-338.15))^(i-1) for i in 1:8)
    return WATER_PC*exp((WATER_TC/T-1)*total)
end
function _water_liquid_estimate(T)
    return WATER_RHOC*(1+sum(_WATER_D[i]*(1-T/WATER_TC)^(i/3) for i in 1:8))
end

# Stable-branch pressure inversion adapted from TPX set_TPp. The volume step
# limits and explicit brackets keep Newton off the unstable van der Waals loop.
function _water_tp_density(T,P,rho)
    v = 1/rho
    vmin,vmax,pmin,pmax = 0.0,1e30,1e30,0.0
    bracket = false
    damping = 1.0
    previous_v,previous_p,previous_residual = v,0.0,0.0
    have_previous = false
    for iteration in 1:220
        current = _water_pressure(T,1/v)
        residual = current-P
        iteration > 1 && abs(residual) < 1e-7*P && return 1/v
        if have_previous && !bracket && residual*previous_residual < 0
            if v < previous_v
                vmin,pmin,vmax,pmax = v,current,previous_v,previous_p
            else
                vmin,pmin,vmax,pmax = previous_v,previous_p,v,current
            end
            bracket = true
        end
        previous_v,previous_p,previous_residual = v,current,residual
        have_previous = true
        dv = .001v*(v <= 1/WATER_RHOC ? -1 : 1)
        slope = (_water_pressure(T,1/(v+dv))-current)/dv
        if current < 0 || slope > 0
            if !bracket
                dv = v < 1/WATER_RHOC ? -.05v : .2v
                vmin > 0 && (dv = .2v)
                vmax < 1e30 && (dv = -.05v)
            else
                slope = (pmax-pmin)/(vmax-vmin)
                v,current = vmax,pmax
                dv = damping*(P-current)/slope
                damping /= 2
            end
        else
            current > P && v > vmin && (vmin = v)
            current < P && v < vmax && (vmax = v)
            v == vmin && (pmin = current)
            v == vmax && (pmax = current)
            vmin < vmax || error("water density bracket collapsed")
            bracket = vmin > 0 && vmax < 1e30
            damping = 1.0
            if slope == 0
                dv = bracket ? .5*(P-current)*(vmax-vmin)/(pmax-pmin) : -.05v
            else
                dv = (P-current)/slope
            end
        end
        limit = .2v
        v < 2/WATER_RHOC && (limit *= .5)
        v < .7/WATER_RHOC && (limit *= .5)
        if bracket && !(vmin <= v+dv <= vmax)
            dv = vmin+(P-pmin)*(vmax-vmin)/(pmax-pmin)-v
        end
        dv = clamp(dv,-limit,limit)
        v += dv
        abs(dv) <= 4eps(Float64)*abs(v) && return 1/v
    end
    error("water density did not converge at $T K and $P Pa")
end

function _water_saturation!(water::PureWater,T)
    isfinite(T) && WATER_TMIN <= T <= WATER_TC ||
        throw(ArgumentError("water saturation temperature must lie in $(WATER_TMIN)–$(WATER_TC) K"))
    T == water.sat_T && return (P=water.sat_P,rhof=water.liquid_density,rhog=water.vapor_density)
    pressure = _water_psat_estimate(T)
    logpressure = log(pressure)
    rhof,rhog = _water_liquid_estimate(T),pressure*WATER_MW/(R*T)
    for iteration in 1:30
        rhof = _water_tp_density(T,pressure,rhof)
        rhog = _water_tp_density(T,pressure,rhog)
        liquid,vapor = _water_raw(T,rhof),_water_raw(T,rhog)
        gf = liquid.u+_water_pressure(T,rhof)/rhof-T*liquid.s
        gg = vapor.u+_water_pressure(T,rhog)/rhog-T*vapor.s
        dg = gg-gf
        if rhog > rhof
            rhof,rhog = rhog,rhof
            dg = -dg
        end
        if abs(dg) < .001
            water.sat_T,water.sat_P,water.liquid_density,water.vapor_density = T,pressure,rhof,rhog
            return (P=pressure,rhof=rhof,rhog=rhog)
        end
        dp = dg/(1/rhog-1/rhof)
        old = pressure
        if abs(dp) > pressure
            logpressure -= dp/pressure
            pressure = exp(logpressure)
        else
            pressure -= dp
            logpressure = log(pressure)
        end
        if pressure > WATER_PC
            pressure = old+.5*(WATER_PC-old)
            logpressure = log(pressure)
        elseif pressure <= 0
            pressure = old/2
            logpressure = log(pressure)
        end
    end
    error("water saturation did not converge at $T K")
end

function PureWater()
    water = PureWater(0.,0.,NaN,NaN,NaN,NaN)
    T = 298.15
    pressure = 1e-5*_water_saturation!(water,T).P
    rho = _water_tp_density(T,pressure,pressure/(_WATER_R*T))
    raw = _water_raw(T,rho)
    a = (4.19864056,-2.0364341e-3,6.52040211e-6,-5.48797062e-9,1.77197817e-12,-3.02937267e4,-.849032208)
    hRT = a[1]+a[2]*T/2+a[3]*T^2/3+a[4]*T^3/4+a[5]*T^4/5+a[6]/T
    sR = a[1]*log(T)+a[2]*T+a[3]*T^2/2+a[4]*T^3/3+a[5]*T^4/4+a[7]-log(pressure/one_atm)
    water.energy_offset = hRT*R*T/WATER_MW-(raw.u+_water_pressure(T,rho)/rho)
    water.entropy_offset = sR*R/WATER_MW-raw.s
    return water
end

"Return saturation pressure [Pa] and saturated liquid/vapor densities [kg/m³]."
water_saturation(water::PureWater,T) = _water_saturation!(water,Float64(T))

function water_saturation_temperature(water::PureWater,P)
    minimum_pressure = _water_saturation!(water,WATER_TMIN).P
    isfinite(P) && minimum_pressure <= P <= WATER_PC ||
        throw(ArgumentError("water saturation pressure must lie between the triple and critical pressures"))
    lower,upper = WATER_TMIN,WATER_TC
    T = lower+(upper-lower)*log(P/minimum_pressure)/log(WATER_PC/minimum_pressure)
    for iteration in 1:70
        pressure = _water_saturation!(water,T).P
        abs(pressure-P) <= 2e-9*P && return T
        pressure > P ? (upper = T) : (lower = T)
        delta = min(1e-3,(WATER_TC-WATER_TMIN)/1e5)
        Tprobe = T+delta <= WATER_TC ? T+delta : T-delta
        slope = (_water_saturation!(water,Tprobe).P-pressure)/(Tprobe-T)
        Tnew = T+(P-pressure)/slope
        T = lower < Tnew < upper ? Tnew : (lower+upper)/2
    end
    error("water saturation temperature did not converge")
end

function _water_state_density(water,T,rho)
    raw = _water_raw(T,rho)
    pressure = _water_pressure(T,rho)
    quality = rho > WATER_RHOC ? 0.0 : 1.0
    u,s = raw.u,raw.s
    h = u+pressure/rho
    if T < WATER_TC
        sat = _water_saturation!(water,T)
        if sat.rhog <= rho <= sat.rhof
            pressure = sat.P
            quality = clamp((1/rho-1/sat.rhof)/(1/sat.rhog-1/sat.rhof),0.,1.)
            liquid,vapor = _water_raw(T,sat.rhof),_water_raw(T,sat.rhog)
            u = (1-quality)*liquid.u+quality*vapor.u
            s = (1-quality)*liquid.s+quality*vapor.s
            h = (1-quality)*(liquid.u+_water_pressure(T,sat.rhof)/sat.rhof)+quality*(vapor.u+_water_pressure(T,sat.rhog)/sat.rhog)
        else
            quality = rho >= sat.rhof ? 0.0 : 1.0
        end
    end
    return WaterState(water,T,pressure,rho,quality,h+water.energy_offset,u+water.energy_offset,s+water.entropy_offset)
end

function _water_state_tq(water,T,Q)
    isfinite(Q) && 0 <= Q <= 1 || throw(ArgumentError("water vapor quality must lie in [0,1]"))
    sat = _water_saturation!(water,T)
    rho = Q == 0 ? sat.rhof : Q == 1 ? sat.rhog : 1/((1-Q)/sat.rhof+Q/sat.rhog)
    return _water_state_density(water,T,rho)
end
function _water_state_tp(water,T,P)
    isfinite(P) && P > 0 || throw(ArgumentError("water pressure must be finite and positive"))
    if T < WATER_TC
        sat = _water_saturation!(water,T)
        abs(P-sat.P) < 1e-8*P && throw(ArgumentError("T and P lie on saturation; specify vapor quality Q"))
        seed = P < sat.P ? sat.rhog : sat.rhof
    else
        seed = max(P/(_WATER_R*T),WATER_RHOC*1.1)
    end
    rho = _water_tp_density(T,P,seed)
    return _water_state_density(water,T,rho)
end

function _water_state_property_pressure(water,target,P,property)
    isfinite(target) || throw(ArgumentError("water target property must be finite"))
    isfinite(P) && P > 0 || throw(ArgumentError("water pressure must be finite and positive"))
    lower,upper = WATER_TMIN,WATER_TMAX
    satT = NaN
    satliquid,satvapor = nothing,nothing
    if _water_saturation!(water,WATER_TMIN).P <= P <= WATER_PC
        satT = water_saturation_temperature(water,P)
        satliquid,satvapor = _water_state_tq(water,satT,0.),_water_state_tq(water,satT,1.)
        f,g = getproperty(satliquid,property),getproperty(satvapor,property)
        if f <= target <= g
            return _water_state_tq(water,satT,(target-f)/(g-f))
        elseif target < f
            upper = satT
        else
            lower = satT
        end
    end
    function evaluate(T)
        T == satT && return upper == satT ? satliquid : satvapor
        return _water_state_tp(water,T,P)
    end
    first,last = evaluate(lower),evaluate(upper)
    fl,fu = getproperty(first,property)-target,getproperty(last,property)-target
    fl <= 0 <= fu || throw(ArgumentError("water target $property is outside the temperature domain at this pressure"))
    scale = property == :s ? max(abs(target),1e3) : max(abs(target),1e6)
    for iteration in 1:90
        width = upper-lower
        T = clamp((lower*fu-upper*fl)/(fu-fl),lower+.02width,upper-.02width)
        state = evaluate(T)
        residual = getproperty(state,property)-target
        (abs(residual) <= 2e-11*scale || width < 1e-8) && return state
        residual > 0 ? ((upper,fu) = (T,residual)) : ((lower,fl) = (T,residual))
    end
    error("water property/pressure inversion did not converge")
end

"""
    water_state(model::PureWater; T, P)
    water_state(model::PureWater; T, Q)
    water_state(model::PureWater; P, Q)
    water_state(model::PureWater; h, P)
    water_state(model::PureWater; s, P)

Construct an equilibrium liquid/vapor water state from exactly two independent
SI properties. `Q` is vapor mass fraction, `h` specific enthalpy and `s` specific
entropy. `(T,v)` and `(u,P)` are also supported. The EOS excludes ice and uses
the temperature limits of the Reynolds model; TP on saturation is ambiguous.
State fields include `T,P,Q,rho,v,h,u,s,cp,cv,g,a`. Two-phase cp is infinite and
cv is NaN, matching the convention used by Cantera's Reynolds water backend.
"""
function water_state(water::PureWater;T=nothing,P=nothing,Q=nothing,h=nothing,s=nothing,u=nothing,v=nothing)
    count(!isnothing,(T,P,Q,h,s,u,v)) == 2 || throw(ArgumentError("specify exactly two independent water properties"))
    if !isnothing(T)
        T = Float64(T)
        isfinite(T) && WATER_TMIN <= T <= WATER_TMAX || throw(ArgumentError("water temperature must lie in $(WATER_TMIN)–$(WATER_TMAX) K"))
        !isnothing(Q) && return _water_state_tq(water,T,Float64(Q))
        !isnothing(P) && return _water_state_tp(water,T,Float64(P))
        if !isnothing(v)
            isfinite(v) && v > 0 || throw(ArgumentError("specific volume must be finite and positive"))
            return _water_state_density(water,T,1/Float64(v))
        end
    elseif !isnothing(P)
        P = Float64(P)
        !isnothing(Q) && return _water_state_tq(water,water_saturation_temperature(water,P),Float64(Q))
        !isnothing(h) && return _water_state_property_pressure(water,Float64(h),P,:h)
        !isnothing(s) && return _water_state_property_pressure(water,Float64(s),P,:s)
        !isnothing(u) && return _water_state_property_pressure(water,Float64(u),P,:u)
    end
    throw(ArgumentError("unsupported water property pair"))
end

function water_cp(state::WaterState)
    water,T,P = state.model,state.T,state.P
    0 < state.Q < 1 && T < WATER_TC && return Inf
    delta = 1e-4*T
    lower,upper = max(WATER_TMIN,T-delta),min(WATER_TMAX,T+delta)
    Tsat = try
        water_saturation_temperature(water,P)
    catch exception
        exception isa ArgumentError || rethrow()
        NaN
    end
    if state.Q == 0 && T < WATER_TC
        upper = isnan(Tsat) ? upper : min(Tsat,upper)
    else
        lower = isnan(Tsat) ? lower : max(Tsat,lower)
    end
    first = lower == Tsat ? _water_state_tq(water,lower,1.) : _water_state_tp(water,lower,P)
    last = upper == Tsat ? _water_state_tq(water,upper,0.) : _water_state_tp(water,upper,P)
    return T*(last.s-first.s)/(upper-lower)
end
function water_cv(state::WaterState)
    water,T = state.model,state.T
    0 < state.Q < 1 && T < WATER_TC && return NaN
    delta = 1e-4*T
    lower,upper = max(WATER_TMIN,T-delta),min(WATER_TMAX,T+delta)
    first,last = _water_state_density(water,lower,state.rho),_water_state_density(water,upper,state.rho)
    if first.Q != state.Q
        lower,first = T,state
    end
    if last.Q != state.Q
        upper,last = T,state
    end
    return T*(last.s-first.s)/(upper-lower)
end

"""
    water_rankine(model; inlet_temperature=300, boiler_pressure=8e5,
                   pump_efficiency=0.6, turbine_efficiency=0.8)

Compute a Rankine cycle starting with saturated liquid, with an adiabatic pump,
constant-pressure heating to saturated vapor, and an adiabatic turbine. Returns
the four states, the two isentropic comparison states, work and heat in J/kg,
and thermal efficiency. Inputs and outputs use SI units.
"""
function water_rankine(water::PureWater;inlet_temperature=300.,boiler_pressure=8e5,
                       pump_efficiency=.6,turbine_efficiency=.8)
    0 < pump_efficiency <= 1 && 0 < turbine_efficiency <= 1 ||
        throw(ArgumentError("pump and turbine efficiencies must lie in (0,1]"))
    first = water_state(water;T=inlet_temperature,Q=0.)
    first.P < boiler_pressure < WATER_PC || throw(ArgumentError("boiler pressure must exceed the inlet pressure and remain subcritical"))
    pumpideal = water_state(water;s=first.s,P=boiler_pressure)
    pump_work = (pumpideal.h-first.h)/pump_efficiency
    second = water_state(water;h=first.h+pump_work,P=boiler_pressure)
    third = water_state(water;P=boiler_pressure,Q=1.)
    expansionideal = water_state(water;s=third.s,P=first.P)
    turbine_work = (third.h-expansionideal.h)*turbine_efficiency
    fourth = water_state(water;h=third.h-turbine_work,P=first.P)
    heat_added = third.h-second.h
    heat_added > 0 || throw(ArgumentError("cycle requires positive boiler heat input"))
    return (states=(first,second,third,fourth),ideal_states=(pumpideal,expansionideal),
            pump_work,turbine_work,heat_added,efficiency=(turbine_work-pump_work)/heat_added)
end

export PureWater,WaterState,water_state,water_saturation,water_saturation_temperature,
       water_cp,water_cv,water_rankine,WATER_TMIN,WATER_TMAX,WATER_TC,WATER_PC,WATER_RHOC,WATER_MW
