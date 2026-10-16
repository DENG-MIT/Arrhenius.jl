# One-dimensional flames

Arrhenius.jl provides planar premixed flames, axisymmetric counterflow flames, and catalytic impinging flames. Flame objects are built from a Solution and solved with solve!; the API reference contains all constructors and diagnostics.

## Free and burner flames

FreeFlame solves a freely propagating planar premixed flame. Its inlet mass flux, and therefore flame speed, is an eigenvalue. BurnerFlame prescribes inlet mass flux. Both use adaptive spatial grids and support mixture-averaged or multicomponent diffusion, with optional Soret diffusion.

BurnerFlame also supports measured temperature profiles through set_temperature_profile!. The default full-knot policy inserts every supplied profile coordinate. On the conservative, non-Soret path, adaptive refinement adds cells where the profile differs from the cellwise linear chord by more than 1% of the largest supplied temperature. Adaptive profiles require conservative discretization without Soret; enabling Soret restores full knots.

## Counterflow and impinging flames

- CounterflowDiffusionFlame models an axisymmetric fuel/oxidizer diffusion flame. Optional radiation models optically thin CO₂/H₂O loss; two-point control supports continuation on a selected branch.
- CounterflowPremixedFlame models opposed premixed streams.
- CounterflowTwinPremixedFlame models a symmetric twin flame from two premixed streams.
- ImpingingJet models an axisymmetric premixed jet directed at an inert wall.

See the [examples index](examples.md#flames) for runnable entry points and required input data.

## Ionized flames

For mechanisms configured with ionized-gas transport, construct an IonizedFlame through FreeFlame or BurnerFlame. Solve first with charged diffusion frozen, then enable the electric field with set_electric_field!(flame, true) and solve again. The model includes charged species and Gauss's law. Soret diffusion, prescribed temperature profiles, and flame restart snapshots are not supported for this type.

## Transport

Multicomponent diffusion and Soret diffusion require collision-integral data exported for the same mechanism. Mixture-averaged diffusion supports mole- or mass-fraction gradients. See [solver options](solver_options.md#transport) for choices and constraints.

## Minimal flame calculation

Run from the repository root with the package project activated:

~~~julia
using Arrhenius

gas = CreateSolution("mechanism/h2o2.yaml")
flame = FreeFlame(gas; T=300.0, P=one_atm,
    X="H2:1.1,O2:1,AR:5", width=0.03)
solve!(flame)
println("Flame speed: ", flame_speed(flame), " m/s")
~~~
