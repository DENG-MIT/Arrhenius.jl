# Arrhenius

Inspired by [`ReactionMechanismSimulator.jl`](https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl), this project aims at developing a mini package for interpreting combustion chemical kinetic models and compute reaction source term.

The name of `Arrhenius.jl` is reflecting the fact that the distinction between combustion and other chemical reacting flow are temperature-dependent kinetics and large activation energy.

## Installation

> pkg> add https://github.com/DENG-MIT/Arrhenius.jl

## Usage

You can start from the example of pyrolysis of JP10 (an aviation fuel power the flight) under the folder of `example`. It will guide you on how to implement the governing equations with a couple of lines of code.

> Currently, the package relies on [`Cantera`](https://github.com/Cantera/cantera) and [`ReacTorch`](https://github.com/DENG-MIT/reactorch) for interpreting the reaction mechanism. If you want to have a try, you don't need to install Cantera and ReacTorch, since there are already some pre-compiled reaction mechanisms under the folder of `mechanism`. Otherwise, you can install [`Cantera`](https://github.com/Cantera/cantera) and [`ReacTorch`](https://github.com/DENG-MIT/reactorch) to compile it using the python script `interpreter.py` under the folder of `mechanism`. You can also ask for help in the discussion forum and our developers can compile the model for you.

## Validation with Cantera

In the example of pyrolysis.jl, we compare the results with Cantera. The example involves solving the equations of mass fractions and temperature under constant pressure.
![val](./example/JP10_pyrolysis.png)
