# Zero-dimensional reactors

IdealGasReactor defines a homogeneous closed reactor at constant pressure or volume, with adiabatic or isothermal energy treatment. It provides the state equations and Jacobian independently of any ODE package. Pass the problem to an integrator through solve_reactor.

## Constant-pressure ignition

From the repository root, activate example/reactors and run this example. The optional environment includes SciMLBase and OrdinaryDiffEqBDF; the helper selects QNDF and connects the Arrhenius Jacobian to the integrator.

~~~julia
using Arrhenius
include("example/reactors/ode_solver.jl")

gas = CreateSolution("mechanism/h2o2.yaml")
reactor = IdealGasReactor(gas; temperature=1001.0, pressure=one_atm,
    mole_fractions=Dict("H2" => 2, "O2" => 1, "AR" => 4))
solution = solve_reactor(reactor, (0.0, 0.001); integrator=native_bdf,
    reltol=1e-8, abstol=1e-14)
final = reactor_properties(reactor, solution.u[end])
println((temperature_K=final.temperature, pressure_Pa=final.pressure))
~~~

reactor_properties reports thermodynamic state and conserved quantities so callers can inspect integration drift. Arrhenius.jl does not select a mandatory integrator. See the [reactor examples](examples.md#reactors).
