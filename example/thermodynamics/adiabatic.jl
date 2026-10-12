# Adiabatic flame temperature and equilibrium composition of a CH4/air
# mixture including graphite formation: the native counterpart of Cantera's
# samples/python/thermo/adiabatic.py example. Fifty equivalence ratios from
# 0.3 to 3.5, 1 kmol of gas with an initially absent graphite phase, at
# 300 K and 101325 Pa, equilibrated at constant enthalpy and pressure.
#
# Usage: julia example/thermodynamics/adiabatic.jl [gas.yaml] [phase.json]
# The phase JSON is written by mechanism/export_condensed_phase.py.

using Arrhenius

"""
    adiabatic_calculation(gas, phase; T=300.0, P=101325.0, points=50,
                          phi=range(0.3, 3.5; length=points),
                          fuel="CH4", oxidizer="O2:1,N2:3.76")

Adiabatic (HP) equilibrium temperatures and combined species amounts for a
fuel/air equivalence-ratio sweep with an initially absent condensed phase,
starting from 1 kmol of gas and zero kmol of condensed phase at `T`, `P`.
Returns `(; phi, T, species_moles)` where `species_moles[k, j]` is the amount
(kmol) of gas species `k` — followed by the condensed species — at `phi[j]`.
"""
function adiabatic_calculation(gas, phase; T=300.0, P=101325.0, points=50,
                               phi=range(0.3, 3.5; length=points),
                               fuel="CH4", oxidizer="O2:1,N2:3.76")
    ratios = collect(phi)
    temperatures = zeros(length(ratios))
    moles = zeros(gas.n_species + 1, length(ratios))
    for (j, equivalence) in enumerate(ratios)
        X = set_equivalence_ratio(gas, equivalence; fuel, oxidizer)
        result = equilibrate(gas, phase; T, P, X, mode=:HP)
        temperatures[j] = result.T
        moles[:, j] = result.species_moles
    end
    return (; phi=ratios, T=temperatures, species_moles=moles)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    gas_path = length(ARGS) >= 1 ? ARGS[1] :
        joinpath(@__DIR__, "..", "..", "mechanism", "gri30.yaml")
    phase_path = length(ARGS) >= 2 ? ARGS[2] :
        joinpath(@__DIR__, "..", "..", "mechanism", "graphite.condensed.json")
    gas = CreateSolution(gas_path)
    phase = StoichiometricCondensedPhase(phase_path)
    result = adiabatic_calculation(gas, phase)
    for (j, equivalence) in enumerate(result.phi)
        println("At phi = ", equivalence, ", Tad = ", result.T[j],
                " K, ", phase.species, " = ", result.species_moles[end, j], " kmol")
    end
end
