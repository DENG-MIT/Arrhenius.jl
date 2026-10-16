"""
    frozen_sound_speed(gas; T, P=one_atm, X)

Ideal-gas sound speed in m/s with composition held fixed, `sqrt((cp/cv)*R*T/MW)`.
"""
function frozen_sound_speed(gas::Solution; T, P=one_atm, X)
    gas.thermo isa IdealGasThermo || throw(ArgumentError("ideal-gas thermo required"))
    isfinite(T) && T > 0 && isfinite(P) && P > 0 || throw(ArgumentError("positive finite T and P required"))
    x = mole_fractions(gas,X)
    cp = cal_cp_mean(gas,T,P,x)
    cp > R || throw(ArgumentError("positive ideal-gas cv required"))
    return sqrt(cp/(cp-R)*R*T/dot(gas.MW,x))
end

"""
    isentropic_state(gas; T, P=one_atm, X, pressure,
                     temperature_bounds=(1,6000))

Find the state at `pressure` (Pa) with the entropy and composition of `(T,P,X)`.
Returns `(T,P,X,Y)` without changing the input. Species heat capacities may vary
with temperature. The explicit temperature bounds constrain the search; values
outside species fit ranges extrapolate their thermodynamic polynomials.
"""
function isentropic_state(gas::Solution; T, P=one_atm, X, pressure,
                          temperature_bounds=(1.,6000.))
    gas.thermo isa IdealGasThermo || throw(ArgumentError("ideal-gas thermo required"))
    all(v -> isfinite(v) && v > 0,(T,P,pressure)) || throw(ArgumentError("positive finite temperature and pressures required"))
    length(temperature_bounds) == 2 || throw(ArgumentError("provide lower and upper temperature bounds"))
    lo,hi = Float64.(temperature_bounds)
    isfinite(lo) && isfinite(hi) && 0 < lo < hi || throw(ArgumentError("temperature bounds must be finite, positive and ordered"))
    lo <= T <= hi || throw(ArgumentError("initial temperature lies outside temperature bounds"))
    x = mole_fractions(gas,X)
    mw = dot(gas.MW,x)
    target = cal_smass_mean(gas,T,P,x)
    residual(t) = cal_smass_mean(gas,t,pressure,x)-target
    lower = upper = Float64(T)
    fl = fu = residual(T)
    while fl > 0 && lower > lo
        lower = max(lo,lower/1.5)
        fl = residual(lower)
    end
    while fu < 0 && upper < hi
        upper = min(hi,upper*1.5)
        fu = residual(upper)
    end
    fl <= 0 <= fu || throw(ArgumentError("isentropic state is not bracketed in $(lo)–$(hi) K"))
    temperature = clamp(Float64(T),lower,upper)
    tolerance = 1e-13*max(abs(target),R/mw,1.)
    for _ in 1:60
        f = residual(temperature)
        if abs(f) <= tolerance
            return (T=temperature,P=Float64(pressure),X=x,Y=x.*gas.MW./mw)
        end
        if f > 0
            upper = temperature
        else
            lower = temperature
        end
        cp = cal_cpmass_mean(gas,temperature,pressure,x)
        cp > 0 || throw(ArgumentError("positive heat capacity required along the isentrope"))
        next = temperature - f*temperature/cp
        temperature = lower < next < upper ? next : (lower+upper)/2
    end
    error("isentropic temperature search did not converge")
end

"""
    equilibrium_sound_speeds(gas; T, P=one_atm, X, pressure_step=1e-4,
                            temperature_bounds=(200,6000))

Equilibrate at `(T,P)` and estimate equilibrium and frozen sound speeds from a
relative pressure perturbation at constant entropy. Returns speeds in m/s as
`equilibrium`, `frozen`, and `frozen_at_equilibrium`, plus the unperturbed `state`
and `perturbed_state`. The last speed uses the analytic frozen formula at the
perturbed equilibrium state. The pressure derivative has first-order truncation
error; `pressure_step` must also be large enough to resolve equilibrium errors.
"""
function equilibrium_sound_speeds(gas::Solution; T, P=one_atm, X, pressure_step=1e-4,
                                  temperature_bounds=(200.,6000.))
    isfinite(pressure_step) && pressure_step > 0 || throw(ArgumentError("positive finite pressure step required"))
    p1 = P*(1+pressure_step)
    p1 > P && isfinite(p1) || throw(ArgumentError("pressure step must produce a distinct finite pressure"))
    gas.thermo isa IdealGasThermo || throw(ArgumentError("ideal-gas thermo required"))
    isfinite(T) && T > 0 && isfinite(P) && P > 0 || throw(ArgumentError("positive finite T and P required"))
    length(temperature_bounds) == 2 || throw(ArgumentError("provide lower and upper temperature bounds"))
    lo,hi = Float64.(temperature_bounds)
    isfinite(lo) && isfinite(hi) && 0 < lo < hi || throw(ArgumentError("temperature bounds must be finite, positive and ordered"))
    system = _equilibrium_system(gas,mole_fractions(gas,X))
    state = _equilibrium_at!(system,Float64(T),Float64(P))
    rho0 = P*dot(gas.MW,state.X)/(R*state.T)
    frozen = isentropic_state(gas;T=state.T,P,X=state.X,pressure=p1,temperature_bounds)
    rho_frozen = p1*dot(gas.MW,frozen.X)/(R*frozen.T)
    # A density increment of order 1e-4 amplifies temperature-search error.
    # Rebase the conserved amounts to the actual unperturbed equilibrium state.
    # Keeping the inlet amounts instead would amplify its tiny balance residual
    # in the pressure derivative. Element potentials remain a valid warm start.
    system.X .= frozen.X
    mul!(system.b,system.A,view(system.X,system.species))
    @. system.logb = log(system.b)
    system.state[end] = 0.
    perturbed_state = _equilibrate(system,frozen.X,frozen.T,p1,:SP,lo,hi,1e-13)
    rho_eq = p1*dot(gas.MW,perturbed_state.X)/(R*perturbed_state.T)
    rho_eq > rho0 && rho_frozen > rho0 || error("pressure perturbation did not resolve positive isentropic compressibility")
    return (equilibrium=sqrt((p1-P)/(rho_eq-rho0)),frozen=sqrt((p1-P)/(rho_frozen-rho0)),
            frozen_at_equilibrium=frozen_sound_speed(gas;T=perturbed_state.T,P=p1,X=perturbed_state.X),
            state,perturbed_state)
end

export frozen_sound_speed, isentropic_state, equilibrium_sound_speeds
