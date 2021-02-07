# Arrhenius

Inspired by [ReactionMechanismSimulator.jl](https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl), this project aims at developing a mini package for interpreting combustion chemical kinetic models.

The package name of `Arrhenius.jl` is reflecting the fact that the distinction between combustion and other chemical reacting flow is temperature dependent kinetics and large activation energy.

## Usage

Currently, the package rely on Cantera and reacTorch for interpreting reaction mechanism. Therefore, you shall install cantera and reactorch first, and run the python script `interpreter.py` under the folder of python.
