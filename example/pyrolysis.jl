using Arrhenius
using Arrhenius:one_atm
using Arrhenius:CreateSolution, species_index, R, one_atm
using LinearAlgebra
using DifferentialEquations
using Sundials
using Plots
using DelimitedFiles

cantera_data = readdlm("pyrolysis.dat")
ct_ts= cantera_data[:, 1]
ct_T = cantera_data[:, 2]
ct_Y = cantera_data[:, 3:end]

gas = CreateSolution("../mechanism/JP10skeletal.yaml")
ns = gas.n_species

Y0 = zeros(ns)
Y0[species_index("C10H16")] = 0.05
Y0[species_index("N2")] = 0.95
T0 = 1200.0
P = one_atm

u0 = vcat(Y0, T0)

using Arrhenius:Y2X, Y2C,
                cp_mole_func, cp_mass_func,
                H_mole_func, S_mole_func, wdot_func

function dudt!(du, u, p, t)
    T = u[end]
    Y = u[1:ns]
    mean_MW = 1. / dot(Y, 1 ./ gas.MW)
    ρ_mass = P / R / T * mean_MW
    X = Y2X(Y, mean_MW)
    C = Y2C(Y, ρ_mass)
    cp_mole = cp_mole_func(T, X)
    cp_mass = cp_mass_func(cp_mole, mean_MW)
    H_mole, h_mole = H_mole_func(T, X)
    S_mole, s_mole, S0 = S_mole_func(T, P, X)
    wdot = wdot_func(T, C, S0, h_mole)
    Ydot = wdot / ρ_mass .* gas.MW
    Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass
    du .= vcat(Ydot, Tdot)
end

tspan = [0.0, 0.07];
prob = ODEProblem(dudt!, u0, tspan);
sol = solve(prob, Rosenbrock23(autodiff=false), reltol=1e-6, abstol=1e-9)

plt = plot(sol.t, sol[species_index("C10H16"), :], lw=2, label="Arrhenius.jl");
plot!(plt, ct_ts, ct_Y[:, species_index("C10H16")], label="Cantera")
ylabel!(plt, "Mass Fraction of C10H16")
xlabel!(plt, "Time [s]")
pltT = plot(sol.t, sol[end, :], lw=2, label="Arrhenius.jl");
plot!(pltT, ct_ts, ct_T, label="Cantera")
ylabel!(pltT, "Temperature [K]")
xlabel!(pltT, "Time [s]")
title!(plt, "JP10 pyrolysis @1200K/1atm")
pltsum = plot(plt, pltT, legend=true, framestyle=:box)
png(pltsum, "JP10_pyrolysis.png")

using ForwardDiff
using ForwardDiff: GradientConfig, Chunk, gradient!
using DiffEqSensitivity

function fsol(u0)
    _prob = remake(prob, u0=u0, tspan = [0.0, 0.07])
    sol = solve(_prob, KenCarp5(autodiff=false),
                reltol=1e-3, abstol=1e-6,
                sensealg=ForwardDiffSensitivity(convert_tspan=false))
    return sol[end, end]
end

out = similar(u0)
cfg = GradientConfig(fsol, u0);
@time gradient!(out, fsol, u0)


# using Profile
# Profile.clear
# @profile ForwardDiff.gradient(x -> fsol(x), u0)
# Juno.profiler(; C = false)

# using Arrhenius:set_states
# function f(u)
#     T = u[end]
#     Y = u[1:ns]
#     mean_MW = 1. / dot(Y, 1 ./ gas.MW)
#     ρ_mass = P / R / T * mean_MW
#     X = Y2X(Y, mean_MW)
#     C = Y2C(Y, ρ_mass)
#     cp_mole = cp_mole_func(T, X)
#     cp_mass = cp_mass_func(cp_mole, mean_MW)
#     H_mole, h_mole = H_mole_func(T, X)
#     S_mole, s_mole, S0 = S_mole_func(T, P, X)
#     # wdot = wdot_func(T, C, S0, h_mole)
#     # Ydot = wdot / ρ_mass .* gas.MW
#     # Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass
#     return S_mole
# end

# # grad = ForwardDiff.gradient(x -> f(x), u0)

# using Profile
# Profile.clear
# Profile.init(n = 10^7, delay = 0.01)
# @profile ForwardDiff.gradient(x -> f(x), u0)
# Juno.profiler(; C = false)
