using Arrhenius
using Arrhenius:one_atm
using Arrhenius:CreateSolution,CreateMSolution,set_states, species_index
using LinearAlgebra
using DifferentialEquations
using Sundials
using Plots
using DelimitedFiles

cantera_data = readdlm("pyrolysis.dat")
ct_ts= cantera_data[:, 1]
ct_T = cantera_data[:, 2]
ct_Y = cantera_data[:, 3:end]

gas = CreateSolution("../python/JP10skeletal.yaml")
mgas = CreateMSolution(gas)

ns = gas.n_species

Y0 = zeros(ns)
Y0[species_index(gas, "C10H16")] = 0.05
Y0[species_index(gas, "N2")] = 0.95
T0 = 1200.0
P = one_atm

set_states(gas, mgas, T0, P, Y0)

@show mgas.ρ_mass
@show mgas.wdot

u0 = vcat(Y0, T0)

function dudt!(du, u, p, t)
    set_states(gas, mgas, u[end], P, @view(u[1:ns]))
    Ydot = mgas.wdot / mgas.ρ_mass .* gas.molecular_weights
    Tdot = -dot(mgas.h_mole, mgas.wdot) / mgas.ρ_mass / mgas.cp_mass
    du .= vcat(Ydot, Tdot)
end

tspan = [0.0, 0.07];
prob = ODEProblem(dudt!, u0, tspan);

sol = solve(prob, CVODE_BDF(), reltol=1e-6, abstol=1e-9)

plt = plot(sol.t, sol[species_index(gas, "C10H16"), :], lw=2, label="Arrhenius.jl");
plot!(plt, ct_ts, ct_Y[:, species_index(gas, "C10H16")], label="Cantera")
ylabel!(plt, "Mass Fraction of C10H16")
xlabel!(plt, "Time [s]")
pltT = plot(sol.t, sol[end, :], lw=2, label="Arrhenius.jl");
plot!(pltT, ct_ts, ct_T, label="Cantera")
ylabel!(pltT, "Temperature [K]")
xlabel!(pltT, "Time [s]")
title!(plt, "JP10 pyrolysis @1200K/1atm")
pltsum = plot(plt, pltT, legend=true, framestyle=:box)
png(pltsum, "JP10_pyrolysis.png")
