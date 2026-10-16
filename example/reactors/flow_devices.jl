# Two fixed-volume vessels with a time-varying mass-flow controller, valve,
# and a conducting wall with prescribed heat flux. All transfer is internal.
using Arrhenius
include("network_cases.jl")
include("network_solver.jl")

gas = CreateSolution(isempty(ARGS) ? joinpath(@__DIR__,"..","..","mechanism","gri30.yaml") : ARGS[1])
network = closed_pair_network(gas)
times = collect(range(0.0,2.0;length=101))
solution = solve_network(network,(0.0,2.0);integrator=native_network_bdf,
    reltol=1e-10,abstol=1e-17,saveat=times,tstops=times)
initial = network_diagnostics(network)
final = network_diagnostics(network,solution.u[end],solution.t[end])
initial_mass = sum(node.mass for node in initial.nodes)
initial_energy = sum(node.total_internal_energy for node in initial.nodes)
println((hot_temperature_K=final.nodes.hot.temperature,cold_temperature_K=final.nodes.cold.temperature,
    mass_drift_kg=sum(node.mass for node in final.nodes)-initial_mass,
    energy_drift_J=sum(node.total_internal_energy for node in final.nodes)-initial_energy,
    final_mass_flow_rates_kg_s=final.mass_flow_rates,final_wall_heat_rates_W=final.wall_heat_rates))
