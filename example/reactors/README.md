# Reactor examples

From the repository root, prepare the optional Julia solver environment:

```sh
julia --project=example/reactors -e 'using Pkg; Pkg.instantiate()'
julia --project=example/reactors example/reactors/constant_pressure_ignition.jl
```

This environment uses the local Arrhenius checkout. Ideal-gas mechanisms require a YAML
file and its prepared `.npz` sidecar, as described in the main README.

`ode_solver.jl` supplies the `native_bdf` integrator for `solve_reactor`.
For supported constant-pressure, adiabatic ideal-gas reactors, it uses a native
analytic Jacobian and a sparse bordered linear solve with KLU. Supported
mechanisms use Float64 NASA7 data, integer stoichiometry and ordinary reaction
orders; Arrhenius, three-body, falloff/Troe and fixed-pressure PLOG rates are
supported. Other mechanisms and constant-volume or isothermal reactors use
the dense Jacobian supplied by `reactor_problem`.

Pass `linear_solver=:dense` to `native_bdf` to select the guarded analytic
Jacobian and signed reactor RHS with QNDF's ordinary dense linear solver for a
supported reactor. `linear_solver=:auto` is the default and retains automatic
structured selection. Unsupported reactors continue to use the existing dense
`reactor_problem` fallback for either selector.

Each integration creates fresh solver workspaces. An already-loaded mechanism
can be reused for subsequent reactor calculations. Use separate workspaces for
concurrent integrations, and do not mutate a mechanism during a solve.

The structured solver uses the dependency versions in this environment.
Other environments retain dense QNDF integration. The core Arrhenius package
continues to accept caller-supplied integrators without requiring these solver
dependencies.

Run the optional solver regressions with:

```sh
julia --project=example/reactors test/structured_reactors.jl
```

Compare structured and dense integration of the detailed n-heptane/air reactor
using a prepared `n-heptane-NUIG-2016.yaml` mechanism and sidecar:

```sh
julia --project=example/reactors example/reactors/preconditioned_integration.jl /path/to/n-heptane-NUIG-2016.yaml
```

The example reuses the loaded mechanism, starts a fresh reactor for each mode,
and returns time, temperature, CO2 and fuel mass-fraction histories.

`plasma.jl` solves an isothermal oxygen glow discharge with separate gas and
electron temperatures. Its `PlasmaMechanism` reads the original plasma YAML
and imported species YAML directly, without a sidecar:

```sh
julia --project=example/reactors example/reactors/plasma.jl /path/to/oxygen-plasma-itikawa.yaml /path/to/species-data
```

The second path locates `nasa_gas.yaml`; it may be omitted when that file is
next to the mechanism. The callable `plasma(mechanism)` reuses the loaded
mechanism and creates fresh state and solver storage on every call. It returns
time and electron mole fraction together with the reactor and solution.
`native_bdf` uses a dense automatic-differentiation Jacobian for this
species-only plasma reactor. Gas and electron temperatures remain fixed.

`nanosecond_pulse_discharge.jl` integrates the methane/air pulse with gas-energy
coupling and a Boltzmann electron-energy distribution:

```sh
julia --project=example/reactors example/reactors/nanosecond_pulse_discharge.jl /path/to/methane-plasma-pavan-2023.yaml /path/to/species-data
```

The species-data directory supplies `gri30.yaml` and `nasa_gas.yaml` imported by
the mechanism. The function `nanosecond_pulse_discharge(mechanism)` reuses a loaded
mechanism and starts fresh physical state and solver storage. It returns time,
temperature, species histories, field values, snapshots and solver statistics;
check `completed` and `error` before using results. The original delayed field
updates and all 901 output times are retained. See the
[full pulse validation](../../validation/results/cantera4_wsl_nanosecond_pulse.json)
for accuracy and timing scope.
