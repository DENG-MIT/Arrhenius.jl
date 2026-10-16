# Stationary mix1 calculation via the native Newton helper in mixing_solver.jl.
# Usage: julia mix1.jl path/to/gri30.yaml path/to/air.yaml
# Both yaml files must be prepared (parsed .npz sidecars present). A separate
# air mechanism is required; the gas model is never substituted for air.
using Arrhenius
include(joinpath(@__DIR__, "mixing_solver.jl"))

length(ARGS) == 2 ||
    error("usage: julia mix1.jl <gri30.yaml> <air.yaml> (prepared yaml paths with .npz sidecars)")
gas = CreateSolution(ARGS[1])
air = CreateSolution(ARGS[2])
result = solve_mixing_network(gas, air)
result.converged ||
    error("mixing network did not converge: $(result.message) " *
          "(fixed-scale residual $(result.residual), physical residual $(result.physical_residual))")
mixed = network_diagnostics(result.network, result.state).nodes.mixer
println((temperature_K=mixed.temperature, pressure_Pa=mixed.pressure,
         mass_fractions=Dict(zip(gas.species_names, mixed.mass_fractions)),
         residual=result.residual, physical_residual=result.physical_residual,
         iterations=result.iterations))
