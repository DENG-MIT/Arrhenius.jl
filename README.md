# Arrhenius

We are in an early-development. Expect some adventures and rough edges.

## Installation

> pkg> add https://github.com/Subsapce-Lab/Arrhenius.jl

## Mechanism preprocessing

`CreateSolution` reads a Cantera YAML file and a same-name `.yaml.npz`
sidecar. Generate a sidecar using Cantera 3.2 or 4.0 dev and NumPy:

```bash
python mechanism/export_sidecar.py path/to/mechanism.yaml
```

The exporter supports elementary and three-body Arrhenius reactions,
Lindemann and Troe falloff reactions, pressure-dependent Arrhenius (PLOG)
reactions, Blowers–Masel rates and explicit reaction orders. Unsupported Cantera rate models are
rejected during preprocessing. Native ideal-gas thermochemistry supports
NASA7, multi-region NASA9, Shomate and constant-cp species models.

The first gas phase may select local or imported species and reaction sections.
Regenerate the sidecar after changing any imported YAML file. Julia resolves
imports relative to the mechanism, then in optional search directories:

```julia
gas = CreateSolution("path/to/mechanism.yaml"; data_paths=["path/to/shared-data"])
```

## Native calculations

The Julia solvers provide ideal-gas equilibrium at TP, TV, HP, UV, SP and SV,
including charged species with conserved net charge;
isentropic states, frozen or equilibrium sound speeds, and constant-pressure stream mixing;
closed constant-pressure and constant-volume reactors; connected stirred
reactors with flow devices and heat-transfer walls; planar premixed free flames
and burner-stabilized flames; counterflow diffusion, opposed premixed and twin
premixed flames; and premixed impinging jets with an inert wall.
Gas/solid TP and HP equilibrium supports one initially absent, fixed-stoichiometry
condensed phase with NASA7 or constant-cp thermodynamics, including graphite formation.
Flame calculations support adaptive
grids, mixture-averaged and multicomponent diffusion, Soret diffusion, and
prescribed burner temperature profiles. Counterflow diffusion flames also
support optically thin CO₂/H₂O radiation. A native pure-water model provides
liquid/vapor states, saturation properties and Rankine-cycle calculations.
`critical_properties` returns critical temperature, pressure, density and
compressibility for eight TPX pure-fluid models.
Redlich–Kwong mixtures support gas and liquid cubic roots, caloric properties,
fugacity coefficients and partial molar properties.
Native Redlich–Kwong reaction rates and constant-volume adiabatic reactors
support nonideal shock-tube ignition calculations.
Reactor networks support prescribed, pressure-driven and inertial pistons,
including changing volumes, pressure work and wall heat transfer.
Ideal-surface chemistry supports elementary and sticking reactions, coverage
dependencies, fixed-stoichiometry solids and isothermal catalytic reactors.
`CoverageThermoModel` evaluates coverage-dependent standard enthalpy, entropy,
heat capacity and Gibbs energy with linear, polynomial, piecewise-linear or
interpolated self and cross interactions.
Isothermal catalytic plug flow supports a direct spatial DAE and a chain of
stirred reactors with species and elemental-flux diagnostics.
Catalytic impinging flames couple gas transport and chemistry to steady surface
coverages at a prescribed wall temperature.
Porous-media transport includes molecular and Knudsen diffusion and Darcy flow
through the native `DustyGasTransport` model.
Chemistry, thermodynamics, transport
evaluation and equation solves run in Julia. Cantera is used to preprocess
mechanisms and generate independent validation data.

```julia
using Arrhenius

gas = CreateSolution("mechanism/h2o2.yaml")
flame = FreeFlame(gas; T=300., P=one_atm, X="H2:1.1,O2:1,AR:5", width=.03)
solve!(flame)
println(flame_speed(flame))  # m/s
save_flame("flame.csv", flame; basis=:mole)
save_flame("flame.npz", flame)
```

For multicomponent or Soret diffusion, also export collision-integral data:

```bash
python mechanism/export_multicomponent.py mechanism/h2o2.yaml mechanism/h2o2.yaml.multicomponent.npz
```

```julia
data = MultiTransportData("mechanism/h2o2.yaml.multicomponent.npz", gas)
set_transport!(flame, :multicomponent; data, soret=true)
solve!(flame; slope=.02, curve=.04)
```

For a prescribed burner temperature, `set_temperature_profile!(burner, z, T;
relative=false, grid_policy=:adaptive)` stores every supplied pair and lets a
conservative, no-Soret flame refine its solution grid. Refinement also bounds the
piecewise-linear temperature's deviation from each cell chord to 1% of the largest
supplied temperature. The default `:full_knots` inserts all profile knots for
conservative flames. Enabling Soret restores those knots and interpolates the
existing state. Restart snapshots retain this policy. The fixed-temperature
example first solves a coarse mixture grid, then refines it using the final
multicomponent criteria before saving either transport stage.
See [complete premixed-flame validation](validation/premixed_flames.md) for
accuracy checks and measured runtimes.

See [premixed flames](example/flames/adiabatic_flame.jl),
[burner flames](example/flames/burner_flame.jl),
[counterflow diffusion flames](example/flames/counterflow_diffusion.jl),
[opposed premixed flames](example/flames/counterflow_premixed.jl),
[twin flames](example/flames/counterflow_twin.jl),
[inert-wall flames](example/flames/counterflow_stagnation.jl),
[closed reactors](example/reactors), and
[thermodynamics](example/thermodynamics) for runnable calculations.
The [nozzle example](example/thermodynamics/isentropic.jl) computes adiabatic
area–Mach curves, and the [sound-speed example](example/thermodynamics/sound_speed.jl)
compares frozen and equilibrium acoustic responses.
The [mixing example](example/thermodynamics/mixing.jl) conserves species and
enthalpy while combining streams, then evaluates the mixture at chemical equilibrium.
Transient reactor examples use a caller-supplied Julia ODE integrator.
The [optional reactor environment](example/reactors/README.md) provides QNDF
with an analytic Jacobian and sparse KLU solve for supported adiabatic,
constant-pressure reactors.

The [reactor mixer example](example/reactors/mix1.jl) solves the stationary
species and energy equations with separate air and fuel mechanisms. Prepare
both sidecars and run:

```bash
python mechanism/export_sidecar.py path/to/gri30.yaml
python mechanism/export_sidecar.py path/to/air.yaml
julia --project=. example/reactors/mix1.jl path/to/gri30.yaml path/to/air.yaml
```

The callable `solve_mixing_network(gas, air)` is provided by
[mixing_solver.jl](example/reactors/mixing_solver.jl).

The complete stationary calculation is checked against Cantera 4.0 on
[WSL](validation/results/cantera4_wsl_mix1.json), Windows, and
[Apple M4](validation/results/cantera4_m4_mix1.json). To reproduce the comparison
and [warm-call timings](validation/results/cantera4_m4_reactor_mixing_timing.json),
set `CANTERA` to the absolute path of Cantera commit
`726522be4e2a13454d8415b7ef799d621f665cf3` and use its Python runtime:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python validation/mixing_reference.py --cantera-source "$CANTERA" --output "$PWD/mixer-reference"
julia --threads=1 --project=. validation/mixing_case.jl "$PWD/mixer-reference/inputs" "$PWD/mixer-native"
python validation/mixing_compare.py --native "$PWD/mixer-native" --reference "$PWD/mixer-reference" --output "$PWD/mixer-comparison.json"
julia --threads=1 --project=. validation/reactor_mixing_timing.jl "$PWD" "$PWD/mixer-reference/inputs" "$PWD/mixer-native" "$PWD/mixer-native-timing"
python validation/reactor_mixing_timing.py --source-root "$PWD" --cantera-source "$CANTERA" --reference "$PWD/mixer-reference" --output "$PWD/mixer-reference-timing"
```

The timed calls include fresh network construction and the full steady solve;
mechanism loading, startup, output validation and diagram/report generation
are excluded. Use new output directories for each run.

The [inertial piston example](example/reactors/custom2.jl) computes ignition
coupled to a wall accelerated by the pressure difference. Its complete 0.5-second
calculation is checked against Cantera 4.0 on
[WSL](validation/results/cantera4_wsl_custom2.json) and
[Apple M4](validation/results/cantera4_m4_custom2.json).
Use a Julia environment containing Arrhenius, SciMLBase and OrdinaryDiffEqBDF.
With `CANTERA` pointing to Cantera source revision
`726522be4e2a13454d8415b7ef799d621f665cf3`, prepare the reference and run:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python validation/custom2_reference.py --cantera-source "$CANTERA" --output custom2-reference
julia --threads=1 --project=REACTOR_ENV validation/custom2_case.jl custom2-reference/inputs custom2-native
julia --threads=1 --project=REACTOR_ENV example/reactors/custom2.jl custom2-reference/inputs custom2-cli.csv
python validation/custom2_compare.py --native custom2-native --reference custom2-reference --cli custom2-cli.csv --output custom2-comparison.json
```

Use new output directories for each run. For the recorded WSL MKL configuration,
add `--lazy-library-receipt validation/results/cantera4_wsl_ic_engine.json` to
the comparison command. The receipt permits only its two hash-pinned numerical
libraries to load during the calculation.

The [engine example](example/reactors/ic_engine.jl) computes eight revolutions
of an n-dodecane engine with prescribed injection, valves and piston motion.
It uses QNDF from OrdinaryDiffEqBDF, with SciMLBase and ForwardDiff in the
caller's Julia environment. Prepare the mechanism and run:

```bash
python validation/ic_engine_case.py engine-input
julia --project=ENGINE_ENV example/reactors/ic_engine.jl engine-input/dodecane_IG.yaml engine.csv
```

The CSV contains crank-angle profiles; the companion TOML contains heat,
pressure work, efficiency and CO estimates integrated over accepted states.
The callable entry `solve_ic_engine_qndf` is provided by
[NativeEngineQNDF](example/reactors/ic_engine_qndf_solver.jl).
The complete calculation has been checked against Cantera 4.0 on
[WSL](validation/results/cantera4_wsl_ic_engine.json) and
[Apple M4](validation/results/cantera4_m4_ic_engine.json), including species and
thermal histories, conservation, and integrated outputs.

The [parallel transport example](example/transport/multiprocessing_viscosity.jl)
computes multicomponent thermal conductivity and viscosity over 5,000 temperatures
for a methane/oxygen/nitrogen mixture. It runs both serial and parallel sweeps,
using independent phase and transport storage for each Julia task:

```bash
python mechanism/export_sidecar.py path/to/gri30.yaml
python mechanism/export_multicomponent.py path/to/gri30.yaml path/to/gri30.yaml.multicomponent.npz
julia --threads=4 --project=. example/transport/multiprocessing_viscosity.jl path/to/gri30.yaml path/to/gri30.yaml.multicomponent.npz
```

The [gas/graphite example](example/thermodynamics/adiabatic.jl) computes adiabatic
equilibrium temperature and all gas/solid species amounts across 50 fuel/air mixtures.
Prepare the gas sidecar and condensed-phase parameters, then run:

```bash
python mechanism/export_condensed_phase.py graphite.yaml graphite --output graphite.condensed.json
julia --project=. example/thermodynamics/adiabatic.jl path/to/gri30.yaml graphite.condensed.json
```

The [CO2 equation-of-state example](example/thermodynamics/equations_of_state.jl)
computes ideal-gas, Redlich–Kwong and Span–Wagner properties over 1–100 bar at
300 K, including the stable vapor/liquid transition. It uses the native Julia
package Clapeyron 0.6.28 in a caller-supplied environment. Prepare Cantera's
`example_data/co2-thermo.yaml` and its sidecar, then export the full Helmholtz
parameters with CoolProp and run the Julia calculation:

```bash
julia --project=EOS_ENV -e 'using Pkg; Pkg.develop(path="."); Pkg.add(PackageSpec(name="Clapeyron", version="0.6.28"))'
python mechanism/export_sidecar.py path/to/co2-thermo.yaml
python mechanism/export_helmholtz.py CO2 carbon-dioxide.json
julia --project=EOS_ENV example/thermodynamics/equations_of_state.jl path/to/co2-thermo.yaml carbon-dioxide.json
```

For catalytic calculations, prepare an ideal-surface parameter archive:

```bash
python mechanism/export_surface.py diamond.yaml diamond_100 --output diamond.surface.npz
```

The [diamond-growth example](example/reactors/diamond_cvd.jl) uses native
surface rates and coverage integration. The [Blowers–Masel example](example/kinetics/blowers_masel.jl)
evaluates reaction rates and activation energies as temperature and enthalpy change.


The [electron-energy distribution example](example/thermodynamics/plasma_eedf.jl)
solves the temporal two-term Boltzmann equation on a prescribed energy grid.
`read_eedf_model` reads elastic/effective, excitation, ionization and attachment
cross sections directly from a Cantera-style plasma YAML file. `EEDFState`
supplies gas conditions and target molecular weights; `solve_eedf` returns
center and edge distributions, electron mobility and convergence diagnostics.
With the Phelps air dataset from Cantera's example data, run:

```bash
julia --project=. example/thermodynamics/plasma_eedf.jl path/to/air-plasma-Phelps.yaml eedf.csv
```

The [oxygen glow-discharge example](example/reactors/plasma.jl) integrates
charged-species kinetics with fixed gas and electron temperatures.
`PlasmaMechanism` reads isotropic plasma YAML and imported species compositions
directly, including two-temperature rates and electron-collision cross sections.
`PlasmaState` provides density-preserving electron-energy changes, and
`PlasmaReactor` supplies a constant-pressure, species-only ODE and Jacobian.
The optional reactor environment supplies its QNDF integrator. No mechanism
sidecar is required for this plasma path:

```bash
julia --project=example/reactors example/reactors/plasma.jl path/to/oxygen-plasma-itikawa.yaml path/to/species-data
```

The species-data directory must contain `nasa_gas.yaml` referenced by the
oxygen mechanism. The YAML files are read by Julia; Cantera is not required
at runtime.

For Boltzmann phases, `PlasmaEnergyReactor` supplies a constant-pressure
mass, enthalpy, and species RHS with cached EEDF and field updates. See
[plasma thermochemistry and energy equations](validation/plasma_thermochemistry.md)
for the API and current validation scope.
The [methane nanosecond-pulse example](example/reactors/nanosecond_pulse_discharge.jl)
integrates the coupled gas energy and charged-species equations through the full
pulse, retaining its prescribed field-update schedule.

For mechanisms with `transport: ionized-gas`, preprocess with
`mechanism/export_sidecar.py`, then evaluate native ionized transport:

```julia
gas = CreateSolution("gri30_ion.yaml")
data = IonTransportData(gas)
workspace = IonTransportWorkspace(data)
X = mole_fractions(gas, Dict("CH4"=>1.0, "O2"=>2.0, "N2"=>7.52))
viscosity, conductivity, electrical = ionized_transport!(workspace, data, one_atm, 1200.0, X)
# workspace.diffusion: m²/s; workspace.mobility: m²/(V s)
```

`ionized_flux!` evaluates diffusion and electric-field drift from midpoint
transport and endpoint mass fractions. This model uses gas-temperature
transport and a fixed electron mobility of 0.4 m²/(V s).

`FreeFlame` and `BurnerFlame` construct an `IonizedFlame` for these mechanisms.
Solve first with charged diffusion frozen, then enable the electric field:

```julia
flame = FreeFlame(gas; T=300.0, X="CH4:1,O2:2,N2:7.52", width=0.05,
                  discretization=:conservative)
solve!(flame; ratio=3.0, slope=0.05, curve=0.1)
set_electric_field!(flame, true)
solve!(flame; ratio=3.0, slope=0.05, curve=0.1)
E = electric_field(flame)  # V/m at each grid point
```

The planar model solves species, energy, continuity, and Gauss's law with
signed charged-species states and ionized-gas transport. Select
`discretization=:conservative` for conservative species and total-enthalpy
balances, or use the default `:finite_difference` formulation. Conservative
calculations initialize through a native frozen finite-difference solve.
Soret diffusion, prescribed temperature profiles, and restart
snapshots are not supported for this flame type.

## Publication

+ [Arrhenius.jl: A Differentiable Combustion Simulation Package](https://arxiv.org/pdf/2107.06172.pdf): overview of Arrhenius.jl and applications in deep mechanism reduction, uncertainty quantification, mechanism tuning and model discovery. [Slides in NCM21](https://www.slideshare.net/WeiqiJi/arrheniusjl-a-differentiable-combustion-simulation-package-248457895), [Vedio for NCM21](https://www.youtube.com/watch?v=X1mwpW78NvA).
+ [Machine Learning Approaches to Learn HyChem Models](https://www.researchgate.net/publication/350890609_Machine_Learning_Approaches_to_Learn_HyChem_Models): demonstrate 1000 times faster than genetic algorithms using commercial software for optimizing complex kinetic models.
+ [Neural Differential Equations for Inverse Modeling in Model Combustors](https://www.researchgate.net/publication/351223124_Neural_Differential_Equations_for_Inverse_Modeling_in_Model_Combustors)
+ [SGD-based Optimization in Modeling Combustion Kinetics: Case Studies in Tuning Mechanistic and Hybrid Kinetic Models](https://doi.org/10.1016/j.fuel.2022.124560)




## Applications

+ **Sensitivity analysis for auto-ignition** | [repo](https://github.com/DENG-MIT/ArrheniusActiveSubspace) | Features: auto-differentiation, multi-threading, sensitivity to all of three Arrhenius params A, b and Ea, active subspace based uncertainty quantification
+ **Sensitivity analysis for one-dimensional flames** | [repo](https://github.com/DENG-MIT/Arrhenius_Flame_1D) | Features: auto-differentiation, multi-threading, sensitivity to all of three Arrhenius params A, b and Ea.
+ **Automonous learn kinetic mechanism using neural network** | [repo](https://github.com/DENG-MIT/CRNN_HyChem) | Features: Chemical Reaction Neural Network (CRNN), Neural Ordinary Differential Equations.
+ **Deep Reduction** | [repo](https://github.com/DENG-MIT/DeepReduction) | Features: Two-stages mechanism reduction with deep learning.

**Examples**

> Note that some of the examples are in development and you can have early access by contacting [Weiqi Ji](mailto:weiqiji@mit.edu)
  + [Pyrolysis of JP10](./example/pyrolysis/pyrolysis.ipynb)
  + [Perfect Stirred Reactor](./example/perfect_stirred_reactor)
  + [Auto-ignition](https://github.com/DENG-MIT/NN-Ignition)
  + [Compute Jacobian using AD](https://gist.github.com/jiweiqi/21b8d149bd95b97d9ae948ab92e446df)

## Relevent packages
+ [ReactionMechanismSimulator.jl](https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl) The amazing Reaction Mechanism Simulator for simulating large chemical kinetic mechanisms
+ [Cantera](https://cantera.org/) A comprehensive C++ based combustion simulation package and with great python interface. Arrhenius relies on Cantera when it is applicable.
