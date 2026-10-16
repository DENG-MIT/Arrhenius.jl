# Arrhenius.jl

Arrhenius.jl is a Julia package for reacting-flow and combustion calculations. It provides native chemistry, thermodynamics, transport, reactors, and one-dimensional flame solvers; Cantera YAML mechanisms are preprocessed into sidecar data files.

## Key capabilities

- Ideal-gas thermochemistry, kinetics, equilibrium, and sound speeds
- Zero-dimensional reactors and connected reactor networks
- Free, burner-stabilized, counterflow, impinging-jet, catalytic, and ionized flames
- Mixture-averaged and multicomponent transport, Soret diffusion, and radiation options
- Surface chemistry, plasma models, porous transport, pure-water properties, and real-gas models
- Julia-native calculations with differentiable components and caller-selected ODE integrators

## Install

In Julia's package manager, enter <code>] add https://github.com/DENG-MIT/Arrhenius.jl</code>. For a local checkout, activate the repository project and run <code>Pkg.instantiate()</code> before using its examples. See the [installation guide](https://deng-mit.github.io/Arrhenius.jl/dev/install/).

## Quick start: constant-pressure ignition

The optional reactor example environment provides SciML's QNDF integrator. From the repository root, run this with <code>julia --project=example/reactors</code>:

~~~julia
using Arrhenius
include("example/reactors/ode_solver.jl")

gas = CreateSolution("mechanism/h2o2.yaml")
reactor = IdealGasReactor(gas; temperature=1001.0, pressure=one_atm,
    mole_fractions=Dict("H2" => 2, "O2" => 1, "AR" => 4))
solution = solve_reactor(reactor, (0.0, 0.001); integrator=native_bdf,
    reltol=1e-8, abstol=1e-14)
@show reactor_properties(reactor, solution.u[end]).temperature
~~~

## Quick start: freely propagating flame

The repository includes the H₂/O₂ mechanism and its preprocessed sidecar. Run from the repository root with <code>julia --project=.</code>:

~~~julia
using Arrhenius

gas = CreateSolution("mechanism/h2o2.yaml")
flame = FreeFlame(gas; T=300.0, P=one_atm,
    X="H2:1.1,O2:1,AR:5", width=0.03)
solve!(flame)
@show flame_speed(flame)
~~~

## Documentation

- [Install and prepare mechanisms](https://deng-mit.github.io/Arrhenius.jl/dev/install/)
- [Flame models and discretization](https://deng-mit.github.io/Arrhenius.jl/dev/flames/)
- [Solver options](https://deng-mit.github.io/Arrhenius.jl/dev/solver_options/)
- [Reactor models](https://deng-mit.github.io/Arrhenius.jl/dev/reactors/)
- [Runnable examples](https://deng-mit.github.io/Arrhenius.jl/dev/examples/)
- [Validation methods](https://deng-mit.github.io/Arrhenius.jl/dev/validation/)
- [API reference](https://deng-mit.github.io/Arrhenius.jl/dev/api/)

## Citation and license

If you use Arrhenius.jl in published work, please cite:

W. Ji, X. Su, B. Pang, S. J. Cassady, A. M. Ferris, Y. Li, Z. Ren, R. Hanson, S. Deng, “Arrhenius.jl: A Differentiable Combustion Simulation Package,” arXiv:2107.06172 [physics.chem-ph], 2021. https://arxiv.org/abs/2107.06172

```bibtex
@misc{ji2021arrhenius,
  title         = {Arrhenius.jl: A Differentiable Combustion Simulation Package},
  author        = {Ji, Weiqi and Su, Xingyu and Pang, Bin and Cassady, Sean Joseph and Ferris, Alison M. and Li, Yujuan and Ren, Zhuyin and Hanson, Ronald and Deng, Sili},
  year          = {2021},
  eprint        = {2107.06172},
  archivePrefix = {arXiv},
  primaryClass  = {physics.chem-ph},
  url           = {https://arxiv.org/abs/2107.06172}
}
```

Related publications are listed in the [documentation](https://deng-mit.github.io/Arrhenius.jl/dev/concepts/#Publications). The package is released under the [MIT License](LICENSE).
