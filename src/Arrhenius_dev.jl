if pwd()[end - 2:end] != "src"
    cd("src")
end

module Arrhenius
    using YAML
    using NPZ
    using LinearAlgebra
    using SparseArrays

    include("Constants.jl")
    include("DataStructure.jl")
    include("Solution.jl")
    include("Magic.jl")
    include("Thermo.jl")
    include("Kinetics.jl")
end


```
Following codes are used during development phase only.
```
using LinearAlgebra
using DifferentialEquations
using Sundials
using Profile
using ForwardDiff
using SparseArrays

R = Arrhenius.R
one_atm = Arrhenius.one_atm

gas = Arrhenius.CreateSolution("../mechanism/gri30.yaml")
const ns = gas.n_species
@show gas.thermo.nasa_high[Arrhenius.species_index(gas, "O2"), :]
Y0 = zeros(ns);
Y0[Arrhenius.species_index(gas, "CH4")] = 0.1;
Y0[Arrhenius.species_index(gas, "O2")] = 0.2;
Y0[Arrhenius.species_index(gas, "N2")] = 0.7;

T0 = 1200.0
P = Arrhenius.one_atm
wdot = Arrhenius.set_states(gas, T0, P, Y0)
@show wdot

# T = T0
# Y = Y0
# mean_MW = 1. / dot(Y, 1 ./ gas.MW)
# ρ_mass = P / R / T * mean_MW
# X =  Arrhenius.Y2X(gas, Y, mean_MW)
# C =  Arrhenius.Y2C(gas, Y, ρ_mass)
# cp_mole, cp_mass =  Arrhenius.get_cp(gas, T, X, mean_MW)
# cv_mole, cv_mass = Arrhenius.get_cv(cp_mole, cp_mass, mean_MW)
# 
# h_mole =  Arrhenius.get_H(gas, T, Y, X)
# u_mole = Arrhenius.get_U(h_mole, Y)
# 
# S0 =  Arrhenius.get_S(gas, T, P, X)
# wdot = Arrhenius.wdot_func(gas.reaction, T, C, S0, h_mole; get_qdot=false)


# u0 = vcat(Y0, T0);
# function f(u)
#     # Arrhenius.set_states(gas, u[end], P, u[1:ns])[10]
#     T = u[end]
#     Y = @view(u[1:ns])
#     mean_MW = 1. / dot(Y, 1 ./ gas.MW)
#     ρ_mass = P / R / T * mean_MW
#     X =  Arrhenius.Y2X(gas, Y, mean_MW)
#     C =  Arrhenius.Y2C(gas, Y, ρ_mass)
#     cp_mole, cp_mass =  Arrhenius.get_cp(gas, T, X, mean_MW)
#     h_mole =  Arrhenius.get_H(gas, T, Y, X)
#     S0 =  Arrhenius.get_S(gas, T, P, X)
#     wdot = Arrhenius.wdot_func(gas.reaction, T, C, S0, h_mole; get_qdot=false)
#     # wdot = gas.reaction.vk * qdot
#     return wdot[end]
# end
# 
# @time f(u0)
# @time grad = ForwardDiff.gradient(f, u0)