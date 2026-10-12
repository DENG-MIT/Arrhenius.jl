# Requires Arrhenius, SciMLBase, and OrdinaryDiffEqBDF in the active environment.
using Arrhenius
include("ode_solver.jl")

gas = CreateSolution(joinpath(@__DIR__, "..", "..", "mechanism", "h2o2.yaml"))
for constraint in (:constant_pressure, :constant_volume)
    reactor = IdealGasReactor(gas; temperature=1400.0, pressure=one_atm,
        mole_fractions=Dict("H2" => 2, "O2" => 1, "AR" => 4),
        constraint, energy=:isothermal)
    solution = solve_reactor(reactor, (0.0, 0.002); integrator=native_bdf,
        reltol=1e-8, abstol=1e-14, saveat=1e-5)
    final = reactor_properties(reactor, solution.u[end])
    println((; constraint, temperature_K=final.temperature, pressure_Pa=final.pressure,
             water_mass_fraction=final.mass_fractions[species_index(gas, "H2O")]))
end
