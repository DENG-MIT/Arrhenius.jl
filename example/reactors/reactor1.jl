# Native Julia calculation for Cantera's reactor1.py initial conditions.
# Usage: julia --project=SOLVER_ENV reactor1.jl PATH_TO_PREPARED_CANTERA_H2O2_YAML
# The mechanism and its .npz sidecar must include H2, O2, and N2. The repository's
# older nine-species h2o2.yaml instead contains Ar; see constant_pressure_ignition.jl.
# Dependencies: Arrhenius, SciMLBase, OrdinaryDiffEqBDF.
using Arrhenius
using Printf
include("ode_solver.jl")

length(ARGS) == 1 || error("usage: reactor1.jl PATH_TO_PREPARED_CANTERA_H2O2_YAML")
gas = CreateSolution(ARGS[1])
reactor = IdealGasReactor(gas; temperature=1001.0, pressure=one_atm,
    mole_fractions=Dict("H2" => 2, "O2" => 1, "N2" => 4),
    constraint=:constant_pressure, energy=:adiabatic)
times = collect(range(0.0, 0.001; length=101))
solution = solve_reactor(reactor, (first(times), last(times)); integrator=native_bdf,
    reltol=1e-8, abstol=1e-14, saveat=times, tstops=times)

println("       t [s]        T [K]        P [Pa]       u [J/kg]")
for (time, state) in zip(solution.t, solution.u)
    properties = reactor_properties(reactor, state)
    @printf("%12.5e %12.4f %13.5f %14.6f\n", time, properties.temperature,
            properties.pressure, properties.internal_energy)
end
