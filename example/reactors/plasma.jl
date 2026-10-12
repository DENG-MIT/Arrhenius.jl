# Isothermal oxygen glow discharge with an isotropic electron distribution.
# Usage: julia --project=example/reactors example/reactors/plasma.jl OXYGEN_YAML [SPECIES_DATA_DIR]
using Arrhenius
include("ode_solver.jl")

function _plasma_source_times(t_end, output_step)
    times = [0.0]
    t = 0.0
    while t < t_end
        next = t + output_step
        next > t || throw(ArgumentError("output step is too small to advance time"))
        t = next
        push!(times, t)
    end
    return times
end

function _plasma_electron_history(solution, mechanism)
    time = copy(solution.t)
    electron = Vector{Float64}(undef, length(time))
    ei = mechanism.electron_index
    for n in eachindex(time)
        y = solution.u[n]
        total = 0.0
        for i in eachindex(y)
            total += y[i] / mechanism.MW[i]
        end
        electron[n] = (y[ei] / mechanism.MW[ei]) / total
    end
    return time, electron
end

"""
    plasma(mechanism; t_end=1e-6, output_step=1e-9, integrator=native_bdf)

Run the oxygen glow-discharge calculation using a loaded isotropic mechanism.
Each call creates fresh state, reactor and solver workspaces. Set the initial
300 K, 0.01 atm composition before raising mean electron energy to 10 eV at fixed
density; the reactor retains the resulting pressure and both temperatures.

Output requests use successive additions of `output_step` until the requested
end is reached or passed, as in the Cantera example. Returns time and electron
mole fraction, together with the reactor and solver solution. The species-only
QNDF solve uses relative/absolute tolerances 1e-9/1e-15.
"""
function plasma(mechanism::PlasmaMechanism; t_end=1e-6, output_step=1e-9,
                integrator=native_bdf)
    all(x -> isfinite(x) && x > 0, (t_end, output_step)) ||
        throw(ArgumentError("t_end and output_step must be positive and finite"))
    state = PlasmaState(mechanism)
    set_plasma_state!(state; temperature=300.0, pressure=0.01 * one_atm,
        mole_fractions=Dict("O2" => 1.0, "e" => 0.005, "O2+" => 0.005))
    set_mean_electron_energy!(state, 10.0)
    reactor = PlasmaReactor(state)
    requested = _plasma_source_times(Float64(t_end), Float64(output_step))
    solution = solve_reactor(reactor, (0.0, last(requested)); integrator,
        reltol=1e-9, abstol=1e-15, saveat=requested,
        save_everystep=false, dense=false, maxiters=100_000)
    time, electron_mole_fraction = _plasma_electron_history(solution, mechanism)
    return (; time, electron_mole_fraction, reactor, solution)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    1 <= length(ARGS) <= 2 || error("usage: plasma.jl OXYGEN_YAML [SPECIES_DATA_DIR]")
    data_paths = length(ARGS) == 2 ? [ARGS[2]] : String[]
    mechanism = PlasmaMechanism(ARGS[1]; data_paths)
    history = plasma(mechanism)
    println((time_s=history.time[end], electron_mole_fraction=history.electron_mole_fraction[end],
        pressure_Pa=history.reactor.pressure, saved_points=length(history.time)))
end
