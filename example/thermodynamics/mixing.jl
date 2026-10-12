using Arrhenius, LinearAlgebra

"Molar and extensive report properties of one ideal-gas mixture state."
function mixture_state_report(gas, T, P, X; moles)
    x = mole_fractions(gas, X)
    mw = dot(gas.MW, x)
    h = cal_h_mean(gas, T, P, x)
    return (T=Float64(T), P=Float64(P), density=P * mw / (R * T), meanMW=mw,
            h=h, u=cal_u_mean(gas, T, P, x), s=cal_s_mean(gas, T, P, x),
            g=cal_g_mean(gas, T, P, x), cp=cal_cp_mean(gas, T, P, x),
            cv=cal_cv_mean(gas, T, P, x), X=x, Y=x .* gas.MW ./ mw,
            mu_RT=cal_g(gas, T, P, x) ./ (R * T),
            moles=Float64(moles), mass=moles * mw, enthalpy=moles * h)
end

"""
    mixing_calculation(gas)

Native port of the Cantera `mixing.py` example. Mix 1 kmol of air
(O2:0.21, N2:0.78, AR:0.01) with the stoichiometric methane amount for
CH4 + 2 O2 -> CO2 + 2 H2O, both streams at 300 K and 1 atm, holding
enthalpy and pressure constant; then equilibrate the frozen mixture at
TP = 300 K, 1 atm. Returns report data `before` and `after` equilibration.
"""
function mixing_calculation(gas)
    air = mole_fractions(gas, "O2:0.21, N2:0.78, AR:0.01")
    nO2 = air[findfirst(==("O2"), gas.species_names)]
    streams = ((T=300.0, P=one_atm, X=air, moles=1.0),
               (T=300.0, P=one_atm, X="CH4:1", moles=0.5 * nO2))
    mixed = mix_constant_pressure(gas, streams)
    before = mixture_state_report(gas, mixed.T, mixed.P, mixed.X; moles=mixed.moles)
    equilibrated = equilibrate(gas; T=mixed.T, P=mixed.P, X=mixed.X, mode=:TP)
    after = mixture_state_report(gas, equilibrated.T, equilibrated.P, equilibrated.X;
                                 moles=mixed.mass / dot(gas.MW, equilibrated.X))
    return (; before, after)
end

function print_mixture_report(gas, label, state)
    println(label, ":")
    println("  T = ", state.T, " K, P = ", state.P, " Pa, density = ",
            state.density, " kg/m^3, mean MW = ", state.meanMW, " kg/kmol")
    println("  molar h = ", state.h, " J/kmol, u = ", state.u,
            " J/kmol, s = ", state.s, " J/(kmol K), g = ", state.g, " J/kmol")
    println("  molar cp = ", state.cp, " J/(kmol K), cv = ", state.cv, " J/(kmol K)")
    println("  moles = ", state.moles, " kmol, mass = ", state.mass,
            " kg, enthalpy = ", state.enthalpy, " J")
    println("  species mole/mass fractions and chemical potentials / RT:")
    for (k, name) in enumerate(gas.species_names)
        state.X[k] < 1e-14 && continue
        println("    ", rpad(name, 6), " X = ", state.X[k], ", Y = ", state.Y[k], ", mu/RT = ", state.mu_RT[k])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    mechanism = isempty(ARGS) ? joinpath(@__DIR__, "..", "..", "mechanism", "gri30.yaml") : ARGS[1]
    gas = CreateSolution(mechanism)
    result = mixing_calculation(gas)
    print_mixture_report(gas, "Mixed state (frozen composition)", result.before)
    print_mixture_report(gas, "Equilibrated at TP = 300 K, 1 atm", result.after)
end
