# Arrhenius.jl

Arrhenius.jl provides Julia-native tools for chemical kinetics and reacting-flow calculations. It reads Cantera YAML mechanisms and uses preprocessed sidecar files for kinetic and transport data. Cantera is not required at calculation time.

## Start here

- [Installation](install.md) explains how to add the package and prepare a mechanism.
- [Mechanisms](mechanisms.md) describes supported mechanism data and sidecars.
- [Zero-dimensional reactors](reactors.md) introduces ignition calculations and ODE solver integration.
- [One-dimensional flames](flames.md) covers flame types, transport, and discretization.
- [Solver options](solver_options.md) documents refinement, time stepping, and transport controls.
- [Examples](examples.md) indexes runnable examples by application.
- [Validation methods](validation.md) explains code-parity and solution-accuracy comparisons.
- [API reference](api.md) lists exported types and functions.

## Capability overview

The package includes ideal-gas thermochemistry, kinetics, equilibrium and sound-speed calculations; constant-pressure and constant-volume reactors; reactor networks; planar and opposed-flow flames; catalytic and ionized models; transport, plasma, porous-media, pure-water, and real-gas tools. See the [examples index](examples.md) for entry points.

## Citation

Ji, W. et al., [“Arrhenius.jl: A Differentiable Combustion Simulation Package”](https://arxiv.org/abs/2107.06172).
