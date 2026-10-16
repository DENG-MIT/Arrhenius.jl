# Requires Arrhenius, SciMLBase, and OrdinaryDiffEqBDF in the active environment.
using Arrhenius
include("ode_solver.jl")

gas = CreateSolution(joinpath(@__DIR__, "..", "..", "mechanism", "h2o2.yaml"))
reactor = IdealGasReactor(gas; temperature=1001.0, pressure=one_atm,
    mole_fractions=Dict("H2" => 2, "O2" => 1, "AR" => 4),
    constraint=:constant_pressure, energy=:adiabatic)
solution = solve_reactor(reactor, (0.0, 0.001); integrator=native_bdf,
    reltol=1e-8, abstol=1e-14, saveat=1e-5)

initial = reactor_properties(reactor)
final = reactor_properties(reactor, solution.u[end])
println((temperature_K=final.temperature, pressure_Pa=final.pressure,
         enthalpy_drift_J_kg=final.enthalpy - initial.enthalpy))
