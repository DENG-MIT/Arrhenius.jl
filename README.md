# Arrhenius

Inspired by [ReactionMechanismSimulator.jl](https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl), this project aims at developing a mini package for interpreting combustion chemical kinetic models.

The package name of `Arrhenius.jl` is reflecting the fact that the distinction between combustion and other chemical reacting flow is temperature dependent kinetics and large activation energy.

## Installation

> pkg> add https://github.com/DENG-MIT/Arrhenius.jl

## Usage

Currently, the package rely on `Cantera` and `ReacTorch` for interpreting reaction mechanism. Therefore, you shall install Cantera and ReacTorch first, and run the python script `interpreter.py` under the folder of `python`.

Sample code in the test file
```Julia
using Arrhenius
using Arrhenius:one_atm,CreateSolution,CreateMSolution,set_states
using Test

gas = CreateSolution("../python/JP10skeletal.yaml")
mgas = CreateMSolution(gas)

set_states(gas, mgas, 1200., one_atm, ones(gas.n_species) ./ gas.n_species)

@show mgas.ρ_mass
@show mgas.wdot

```