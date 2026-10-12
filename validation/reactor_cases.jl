# Usage: julia --project=SOLVER_ENV validation/reactor_cases.jl REFERENCE_DIR [MECHANISM_DIR]
# Requires Arrhenius, NPZ, SciMLBase, and OrdinaryDiffEqBDF.
using Arrhenius
using LinearAlgebra
using NPZ
using Printf
using Test
using TOML
include(joinpath(@__DIR__, "..", "example", "reactors", "ode_solver.jl"))

function validate_reactor_cases(reference_directory, mechanism_directory)
    reports = Dict{String,Any}[]
    references = sort(filter(name -> startswith(name, "reference_") && endswith(name, ".npz"),
                             readdir(reference_directory)))
    isempty(references) && error("no reference_*.npz trajectories found")
    @testset "native Julia reactor trajectories" begin
        for reference in references
            data = npzread(joinpath(reference_directory, reference))
            mechanism = String(vec(data["mechanism_utf8"]))
            gas = CreateSolution(joinpath(mechanism_directory, mechanism))
            constraint = only(data["constraint"]) == 0 ? :constant_pressure : :constant_volume
            energy = only(data["energy"]) == 0 ? :adiabatic : :isothermal
            times, expected = vec(data["time"]), data["state"]
            reactor = IdealGasReactor(gas; temperature=expected[end, 1],
                pressure=first(data["pressure"]), mass_fractions=expected[1:end-1, 1],
                constraint, energy)
            solution = solve_reactor(reactor, (first(times), last(times)); integrator=native_bdf,
                reltol=1e-9, abstol=1e-15, saveat=times, tstops=times, maxiters=1_000_000)
            states = reduce(hcat, solution.u)
            @test size(states) == size(expected)
            properties = [reactor_properties(reactor, state) for state in solution.u]
            initial = reactor_properties(reactor)
            temperature_error = maximum(abs, states[end, :] - expected[end, :])
            composition_error = maximum(abs, states[1:end-1, :] - expected[1:end-1, :])
            pressure_error = maximum(abs.([p.pressure for p in properties] ./ data["pressure"] .- 1))
            density_error = maximum(abs.([p.density for p in properties] ./ data["density"] .- 1))
            mass_drift = maximum(abs(p.mass_fraction_sum - 1) for p in properties)
            element_drift = maximum(norm(p.elemental_inventory - initial.elemental_inventory, Inf)
                                    for p in properties)
            invariant = constraint === :constant_pressure ? :enthalpy : :internal_energy
            energy_scale = max(abs(getproperty(initial, invariant)), 1e6)
            energy_drift = maximum(abs(getproperty(p, invariant) - getproperty(initial, invariant))
                                   for p in properties) / energy_scale
            rhs = reactor_rhs(reactor)
            derivative = similar(reactor_state(reactor))
            species_rhs_error, temperature_rhs_error = 0.0, 0.0
            for index in axes(expected, 2)
                rhs(derivative, view(expected, :, index), nothing, times[index])
                target = view(data["rhs"], :, index)
                species_rhs_error = max(species_rhs_error,
                    maximum(abs, derivative[1:end-1] - target[1:end-1]))
                temperature_rhs_error = max(temperature_rhs_error, abs(derivative[end] - target[end]))
            end
            # Net source terms approach zero at equilibrium despite large gross
            # reaction rates. Scale each equation family by its trajectory peak
            # to avoid dividing cancellation roundoff by a vanishing net rate.
            rhs_error = max(species_rhs_error / max(maximum(abs, data["rhs"][1:end-1, :]), 1),
                            temperature_rhs_error / max(maximum(abs, data["rhs"][end, :]), 1))
            @testset "$reference" begin
                @test temperature_error < 0.3
                @test composition_error < 2e-5
                @test pressure_error < 1e-4
                @test density_error < 1e-4
                @test mass_drift < 5e-10
                @test element_drift < 5e-11
                @test minimum(states[1:end-1, :]) > -1e-12
                @test rhs_error < 1e-8
                if energy === :adiabatic
                    @test energy_drift < 1e-6
                else
                    @test maximum(abs, states[end, :] .- initial.temperature) < 1e-7
                end
            end
            name = replace(reference, "reference_" => "", ".npz" => "")
            report = Dict{String,Any}(
                "case" => name, "cantera_version" => String(vec(data["cantera_version_utf8"])),
                "solver" => "OrdinaryDiffEqBDF.QNDF", "temperature_error_K" => temperature_error,
                "mass_fraction_error" => composition_error, "pressure_relative_error" => pressure_error,
                "density_relative_error" => density_error, "mass_drift" => mass_drift,
                "element_drift_kmol_kg" => element_drift, "rhs_relative_error" => rhs_error,
                "species_rhs_absolute_error_per_s" => species_rhs_error,
                "temperature_rhs_absolute_error_K_s" => temperature_rhs_error,
                "final_temperature_K" => states[end, end], "saved_points" => length(times),
                "minimum_mass_fraction" => minimum(states[1:end-1, :]))
            if energy === :adiabatic
                report["energy_relative_drift"] = energy_drift
            end
            push!(reports, report)
            npzwrite(joinpath(reference_directory, "native_$(name).npz"),
                     Dict("time" => solution.t, "state" => states))
            @printf("%-24s ΔT=%9.3g K  ΔY=%9.3g  mass=%9.3g  elements=%9.3g\n",
                    name, temperature_error, composition_error, mass_drift, element_drift)
        end
    end
    open(joinpath(reference_directory, "native_report.toml"), "w") do io
        TOML.print(io, Dict("julia_version" => string(VERSION), "cases" => reports))
    end
    return reports
end

if abspath(PROGRAM_FILE) == @__FILE__
    1 <= length(ARGS) <= 2 || error("usage: reactor_cases.jl REFERENCE_DIR [MECHANISM_DIR]")
    mechanism_directory = length(ARGS) == 2 ? ARGS[2] : joinpath(@__DIR__, "..", "mechanism")
    validate_reactor_cases(ARGS[1], mechanism_directory)
end
