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
#
# using LinearAlgebra
# using DifferentialEquations
# using Sundials
# using Profile
# using ForwardDiff
# using Zygote
#
# gas = Arrhenius.CreateSolution("./mechanism/JP10skeletal.yaml")
#
# ns = gas.n_species
# Y0 = zeros(ns)
# Y0[1] = 0.9
# Y0[Arrhenius.species_index("N2")] = 0.1
# T0 = 1200.0
# P = Arrhenius.one_atm
# wdot, qdot, kf, kr = Arrhenius.set_states(T0, P, Y0)
# @show wdot
#
# function profile_test(n)
#     for i = 1:n
#         Arrhenius.set_states(T0, P, Y0)
#     end
# end
# @time profile_test(1000)
#
# u0 = vcat(Y0, T0)
#
# T = u0[end]
# Y = u0[1:ns]
# mean_MW = 1. / dot(Y, 1 ./ gas.MW)
# ρ_mass = P / Arrhenius.R / T * mean_MW
# X = Arrhenius.Y2X(Y, mean_MW)
# C = Arrhenius.Y2C(Y, ρ_mass)
# cp_mole = Arrhenius.cp_mole_func(T, X)
# cp_mass = Arrhenius.cp_mass_func(cp_mole, mean_MW)
# H_mole, h_mole = Arrhenius.H_mole_func(T, X)
# S_mole, s_mole, S0 = Arrhenius.S_mole_func(T, P, X)
# wdot, qdot, kf, kr = Arrhenius.wdot_func(T, C, S0, h_mole)
#
# function f(u)
#     T = u[end]
#     Y = u[1:ns]
#     mean_MW = 1. / dot(Y, 1 ./ gas.MW)
#     ρ_mass = P / Arrhenius.R / T * mean_MW
#     X = Arrhenius.Y2X(Y, mean_MW)
#     C = Arrhenius.Y2C(Y, ρ_mass)
#     wdot, qdot, kf, kr = Arrhenius.set_states(u[end], P, u[1:ns])
#     return qdot[10]
# end
# grad = ForwardDiff.gradient(x -> f(x), u0)
# @show grad

# @time grad = ForwardDiff.gradient(x -> f(x), u0)
#
# Profile.clear
# # @profile profile_test(1000)
# @profile grad = ForwardDiff.gradient(x -> f(x), u0)
# Juno.profiler(; C = false)
