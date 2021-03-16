module Arrhenius
    using YAML
    using NPZ
    using LinearAlgebra

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
# using LinearAlgebra
# using DifferentialEquations
# using Sundials
# using Profile
# using ForwardDiff
# using Zygote

# gas = Arrhenius.CreateSolution("./mechanism/gri30.yaml")
# ns = gas.n_species
# @show gas.thermo.nasa_high[Arrhenius.species_index(gas, "O2"), :]
# Y0 = zeros(ns)
# Y0[1] = 0.9
# Y0[Arrhenius.species_index(gas, "N2")] = 0.1
# T0 = 1200.0
# P = Arrhenius.one_atm
# wdot = Arrhenius.set_states(gas, T0, P, Y0)
# @show wdot

# function profile_test(n)
#     for i = 1:n
#         Arrhenius.set_states(gas, T0, P, Y0)
#     end
# end
# @time profile_test(1000)
#
# # Profile.clear
# # @profile profile_test(1000)
# # Juno.profiler(; C = false)
#
# u0 = vcat(Y0, T0);
# R = Arrhenius.R
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
#     wdot = Arrhenius.wdot_func(gas.reaction, T, C, S0, h_mole)
#     return wdot[end]
# end
# u0[end] = 1200 + rand()
# @time f(u0)
# @time f(u0)
#
# @time grad = ForwardDiff.gradient(f, u0)
# @time grad = ForwardDiff.gradient(f, u0)

# Profile.clear
# @profile grad = ForwardDiff.gradient(x -> f(x), u0)
# Juno.profiler(; C = false)
