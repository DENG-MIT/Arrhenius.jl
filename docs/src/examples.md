# Examples

The repository examples are executable entry points. Some require prepared mechanism sidecars or an optional Julia environment; each script describes its inputs near the top. Run commands from the repository root unless noted.

## Reactors

- [Constant-pressure ignition](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/reactors/constant_pressure_ignition.jl): adiabatic H₂/O₂/Ar ignition with QNDF.
- [Constant-volume ignition](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/reactors/constant_volume_ignition.jl): closed constant-volume ignition.
- [Isothermal reactor](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/reactors/isothermal.jl): compare constant-pressure and constant-volume reactors.
- [Reactor networks](https://github.com/DENG-MIT/Arrhenius.jl/tree/main/example/reactors): flow devices, stream mixing, heat-transfer walls, and pistons.
- [Engine](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/reactors/ic_engine.jl): a moving-piston model with prescribed injection and valves.
- [Plasma reactors](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/reactors/plasma.jl) and [nanosecond-pulse discharge](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/reactors/nanosecond_pulse_discharge.jl).

## Flames

- [Free premixed flame](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/adiabatic_flame.jl): transport-stage sequence for H₂/O₂/Ar.
- [Burner-stabilized flame](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/burner_flame.jl): prescribed inlet mass flux.
- [Fixed-temperature burner flame](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/flame_fixed_T.jl): measured temperature profile.
- [Counterflow diffusion](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/counterflow_diffusion.jl), [opposed premixed](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/counterflow_premixed.jl), and [twin premixed](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/counterflow_twin.jl).
- [Inert-wall impinging flame](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/counterflow_stagnation.jl) and [catalytic impinging flame](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/catalytic_combustion.jl).
- [Initial guess and restart workflow](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/flame_initial_guess.jl) and [ionized flames](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/flames/source_ion_flame_sequence.jl).

## Thermodynamics, transport, and chemistry

- [Thermodynamic calculations](https://github.com/DENG-MIT/Arrhenius.jl/tree/main/example/thermodynamics): equilibrium, sound speeds, stream mixing, and nozzles.
- [Gas and condensed-phase equilibrium](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/thermodynamics/adiabatic.jl): gas/graphite equilibrium.
- [Real-gas properties](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/thermodynamics/equations_of_state.jl): Redlich–Kwong and Helmholtz-based calculations.
- [Parallel transport sweep](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/transport/multiprocessing_viscosity.jl): viscosity and thermal conductivity over a temperature range.
- [Surface chemistry and catalytic flow](https://github.com/DENG-MIT/Arrhenius.jl/tree/main/example/reactors): surface kinetics, catalytic reactors, and diamond growth.
- [Electron-energy distribution](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/thermodynamics/plasma_eedf.jl) and [Blowers–Masel rates](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/example/kinetics/blowers_masel.jl).
- [Mechanism preprocessing](https://github.com/DENG-MIT/Arrhenius.jl/tree/main/mechanism): exporters for kinetics, transport, surface, condensed-phase, and real-gas data.
