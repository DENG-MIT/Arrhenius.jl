using Arrhenius
using LinearAlgebra
using DifferentialEquations
using ForwardDiff
using DiffEqSensitivity
using Sundials
using Plots
using DelimitedFiles
using Profile

cd("example")
cantera_data = readdlm("pyrolysis.dat")
ct_ts= cantera_data[:, 1]
ct_T = cantera_data[:, 2]
ct_Y = cantera_data[:, 3:end]

gas = CreateSolution("../mechanism/JP10skeletal.yaml")
ns = gas.n_species

Y0 = zeros(ns)
Y0[species_index(gas, "C10H16")] = 0.05
Y0[species_index(gas, "N2")] = 0.95
T0 = 1200.0
P = one_atm

u0 = vcat(Y0, T0)
@inbounds function dudt!(du, u, p, t)
    T = u[end]
    Y = @view(u[1:ns])
    mean_MW = 1. / dot(Y, 1 ./ gas.MW)
    ρ_mass = P / R / T * mean_MW
    X = Y2X(gas, Y, mean_MW)
    C = Y2C(gas, Y, ρ_mass)
    cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
    h_mole = get_H(gas, T, Y, X)
    S0 = get_S(gas, T, P, X)
    wdot = wdot_func(gas.reaction, T, C, S0, h_mole)
    Ydot = wdot / ρ_mass .* gas.MW
    Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass
    du .= vcat(Ydot, Tdot)
end

@inbounds function dudt(u)
    T = u[end]
    Y = u[1:end-1]
    mean_MW = 1. / dot(Y, 1 ./ gas.MW)
    ρ_mass = P / R / T * mean_MW
    X = Y2X(gas, Y, mean_MW)
    C = Y2C(gas, Y, ρ_mass)
    cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
    h_mole = get_H(gas, T, Y, X)
    S0 = get_S(gas, T, P, X)
    wdot = wdot_func(gas.reaction, T, C, S0, h_mole)
    Ydot = wdot / ρ_mass .* gas.MW
    Tdot = -dot(h_mole, wdot) / ρ_mass / cp_mass
    return Tdot
end
println("timing dudt ...")
@time dudt(u0)
@time dudt(u0)
@time du0 = ForwardDiff.gradient(dudt, u0)
@time du0 = ForwardDiff.gradient(dudt, u0)
# timing dudt ...
# 1.949944 seconds (5.23 M allocations: 262.966 MiB, 3.79% gc time)
# 0.000183 seconds (48 allocations: 22.797 KiB)
# 4.014816 seconds (9.36 M allocations: 458.699 MiB, 2.71% gc time)
# 0.001779 seconds (407 allocations: 981.609 KiB)

tspan = [0.0, 0.07];
prob = ODEProblem(dudt!, u0, tspan);
sol = solve(prob, TRBDF2(), reltol=1e-6, abstol=1e-9)

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

# sensealg = InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false))
# sensealg = BacksolveAdjoint(autojacvec=ReverseDiffVJP(false))
# sensealg = ForwardDiffSensitivity()
sensealg = SensitivityADPassThrough()
# sensealg = ForwardSensitivity(autojacvec=true)
# alg = Rosenbrock23(autodiff=true)
alg = TRBDF2()
# alg = Tsit5()
function fsol(u0)
    sol = solve(prob, u0=u0, alg, tspan = (0.0, 7.e-2),
                reltol=1e-3, abstol=1e-6, sensealg=sensealg)
    return sol[end, end]
end
u0[end] = 1200.0 + rand()
println("timing ode solver ...")
@time fsol(u0)
@time fsol(u0)
@time ForwardDiff.gradient(fsol, u0)
@time ForwardDiff.gradient(fsol, u0)
# timing ode solver ...
# 0.405083 seconds (614.32 k allocations: 45.126 MiB)
# 0.036229 seconds (16.72 k allocations: 11.618 MiB)
# 41.850114 seconds (90.09 M allocations: 4.771 GiB, 3.85% gc time)
# 1.517267 seconds (183.25 k allocations: 864.085 MiB, 7.46% gc time)

# Profile.clear
# @profile ForwardDiff.gradient(x -> fsol(x), u0)
# Juno.profiler(; C = false)
