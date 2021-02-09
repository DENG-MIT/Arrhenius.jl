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
# Y0 = ones(ns) / ns
# T0 = 1200.0
# P = Arrhenius.one_atm
# wdot = Arrhenius.set_states(T0, P, Y0)
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
# function f(u)
#     wdot = Arrhenius.set_states(u[end], P, u[1:ns])
#     return wdot[1]
# end
# @time grad = ForwardDiff.gradient(x -> f(x), u0)
#
# Profile.clear
# # @profile profile_test(1000)
# @profile grad = ForwardDiff.gradient(x -> f(x), u0)
# Juno.profiler(; C = false)
