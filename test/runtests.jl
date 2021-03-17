using Arrhenius
using ForwardDiff
using LinearAlgebra
using Test

@testset "jl" begin
    # Write your tests here.

    gas = CreateSolution("../mechanism/gri30.yaml")
    ns = gas.n_species

    Y0 = ones(ns) ./ ns
    T0 = 1200.0
    P = one_atm
    wdot = set_states(gas, T0, P, Y0)
    @show size(wdot)
    @show wdot[1:3]

    u0 = vcat(Y0, T0)
    function f(u)
        T = u[end]
        Y = @view(u[1:ns])
        mean_MW = 1.0 / dot(Y, 1 ./ gas.MW)
        ρ_mass = P / R / T * mean_MW
        X = Y2X(gas, Y, mean_MW)
        C = Y2C(gas, Y, ρ_mass)
        cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
        h_mole = get_H(gas, T, Y, X)
        S0 = get_S(gas, T, P, X)
        wdot = wdot_func(gas.reaction, T, C, S0, h_mole)
        Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass
        return Tdot
    end
    grad = ForwardDiff.gradient(f, u0)
    @show size(grad)
    @show grad[1:3]

end
