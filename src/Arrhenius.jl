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
#
# gas = Arrhenius.CreateSolution("./mechanism/JP10skeletal.yaml")
# mgas = Arrhenius.CreateMSolution(gas)
#
# ns = gas.n_species
#
# Y0 = zeros(ns)
# Y0[1] = 1.0
# T0 = 1200.0
# P = Arrhenius.one_atm
# Arrhenius.set_states(gas, mgas, T0, P, Y0)
# @show mgas.ρ_mass mgas.wdot
#
# function profile_test(n)
#     for i = 1:n
#         Arrhenius.set_states(gas, mgas, T0, P, Y0)
#     end
# end
# @time profile_test(1000)
#
# Profile.clear
# @profile profile_test(1000)
# Juno.profiler(; C = false)
