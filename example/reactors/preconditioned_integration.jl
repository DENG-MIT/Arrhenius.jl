# Compare native structured and dense reactor integrations.
# Usage: julia --project=example/reactors preconditioned_integration.jl PREPARED_N_HEPTANE_YAML
using Arrhenius
include("ode_solver.jl")

"""
    preconditioned_integration(gas; linear_solver=:auto, t_end=0.1, saveat=nothing)

Integrate an adiabatic, constant-pressure n-heptane/air reactor at 1000 K and
one atmosphere. `gas` is an already-loaded mechanism containing NC7H16, O2,
N2 and CO2. Each call prepares an independent reactor and solver workspace.

`:auto` uses the supported structured linear solve; `:dense` uses a dense
linear solve with the same native chemistry and analytic Jacobian. Returns
time, temperature and CO2/fuel mass-fraction histories, plus the reactor and
solver solution for further analysis. `saveat` optionally selects output times.
"""
function preconditioned_integration(gas; linear_solver=:auto, t_end=0.1,
                                    saveat=nothing, callback=nothing)
    isfinite(t_end) && t_end > 0 || throw(ArgumentError("t_end must be positive and finite"))
    co2 = findfirst(==("CO2"), gas.species_names)
    fuel = findfirst(==("NC7H16"), gas.species_names)
    isnothing(co2) && throw(ArgumentError("mechanism must contain CO2"))
    isnothing(fuel) && throw(ArgumentError("mechanism must contain NC7H16"))
    reactor = IdealGasReactor(gas; temperature=1000.0, pressure=one_atm,
        mole_fractions=Dict("NC7H16"=>1.0, "O2"=>11.0, "N2"=>41.36),
        constraint=:constant_pressure, energy=:adiabatic)
    initial_mass = 0.1 * reactor.density
    species_atol = 1e-17 .* gas.MW ./ initial_mass
    output_options = isnothing(saveat) ? (;) : (; saveat)
    solution = solve_reactor(reactor, (0.0, Float64(t_end)); integrator=native_bdf,
        linear_solver=linear_solver, reltol=1e-11,
        abstol=vcat(species_atol, 1e-17), maxiters=1_000_000,
        callback=callback, output_options...)
    time = copy(solution.t)
    temperature = [u[end] for u in solution.u]
    co2_y = [u[co2] for u in solution.u]
    fuel_y = [u[fuel] for u in solution.u]
    return (; time, T=temperature, CO2_Y=co2_y, NC7H16_Y=fuel_y,
            reactor, solution)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 1 || error("usage: preconditioned_integration.jl PREPARED_N_HEPTANE_YAML")
    gas = CreateSolution(ARGS[1])
    for mode in (:auto, :dense)
        history = preconditioned_integration(gas; linear_solver=mode)
        println((linear_solver=mode, time_s=history.time[end],
                 temperature_K=history.T[end], CO2_Y=history.CO2_Y[end],
                 NC7H16_Y=history.NC7H16_Y[end],
                 saved_points=length(history.time),
                 solver_stats=history.solution.destats))
    end
end
