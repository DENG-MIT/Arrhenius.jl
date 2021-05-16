# Get Started

## Basics

We first have to make Arrhenius.jl available in our code (along with other packages you might use):
```julia
using Arrhenius
```
The next step is to create the solution structure similar to Cantera by providing a YAML file as the input. The complete definition of the solution structure can be found [here](https://github.com/DENG-MIT/Arrhenius.jl/blob/a6e77fa501f8e1bfb0b4bd244b507f0bc10f1f8c/src/DataStructure.jl#L29). 
```julia
gas = CreateSolution("../../mechanism/JP10skeletal.yaml")

```
Note that you have to provide the appropriate location of the YAML file. Here we have in the ```../mechanism``` folder. You can find the number of species and reactions in the mechanism you provided using: 
```julia

ns = gas.n_species
nr = gas.n_reactions

```
To view the species that your mechanism uses:
```julia
species_arr = gas.species_names
```
You might also want to access a particular species' data. For instance, N2's index can be accessed using: 
```julia
index_N2 = species_index(gas,"N2")
```
To get the molecular weights of each gas: 
```julia
mol_wt_arr = gas.MW
```
Hence the mean molecular weight can be obtained as: 
```julia
mean_MW = 1. / dot(Y, 1 ./ gas.MW)
#get the density using the ideal gas law
ρ_mass = P / R / T * mean_MW
```
Note that one has to include ```using LinearAlgebra``` in the code to be able to use ```dot()```. One can easily convert between the mass and molar fractions:
```julia
X = Y2X(gas, Y, mean_MW)
C = Y2C(gas, Y, ρ_mass)
```
Some other commonly used thermodynamic functions are given below: 
```julia
#molar specific heats
cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
#molar enthalpies
h_mole = get_H(gas, T, Y, X)
#entropy
S0 = get_S(gas, T, P, X)
```
One of the core functionalities of the Arrhenius.jl is its ability to compute the source term that frequently appears in the governing equations of chemical systems. One can simply compute this using : 
```julia
w_dot = wdot_func(gas.reaction, T, C, S0, h_mole)
```

## Input Files

Similar to Chemkin and Canetra, all calculations in Arrhenius.jl require an input file to describe the properties of the relevant phase(s) of mixture. We adopt the `YAML` format maintained in the Cantera community.

> Currently, the package relies on [`Cantera`](https://github.com/Cantera/cantera) and [`ReacTorch`](https://github.com/DENG-MIT/reactorch) for interpreting the reaction mechanism. The kinetic info in the mechanism files are inpreprested and saved into a `.npz` file with the same name as the `.yaml` file. If you want to have a try, you don't need to install Cantera and ReacTorch, since there are already some pre-compiled reaction mechanisms under the folder of `mechanism`. 

> Otherwise, you can install [`Cantera`](https://github.com/Cantera/cantera) and [`ReacTorch`](https://github.com/DENG-MIT/reactorch) to compile it using the python script `interpreter.py` under the folder of `mechanism`.

> You can also ask for help in the discussion forum and our developers can compile the model for you.

With Cantera and ReacTorch installed, you can use the following python command in your terminal to generate the `.npz` file.

`python interpreter.py -i gri30.yaml`