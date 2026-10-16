# Usage: julia --project=SOLVER_ENV validation/reactor_network_cases.jl REFERENCE_DIRECTORY
using Arrhenius
using LinearAlgebra
using NPZ
using Printf
using Test
using TOML
include(joinpath(@__DIR__,"..","example","reactors","network_solver.jl"))
include(joinpath(@__DIR__,"..","example","reactors","network_cases.jl"))

function validate_network_cases(directory)
    gas = CreateSolution(joinpath(directory,"mechanisms","gri30.yaml"))
    air = CreateSolution(joinpath(directory,"mechanisms","air.yaml"))
    reports = Dict{String,Any}[]
    @testset "native reactor networks versus Cantera 4" begin
        for name in ("mix1","combustor","closed_pair")
            reference = npzread(joinpath(directory,"reference_network_$name.npz"))
            expected, times = reference["state"], vec(reference["time"])
            network = if name == "mix1"
                mixing_network(gas,air)
            elseif name == "combustor"
                Y = expected[1:end-1,1] / sum(expected[1:end-1,1])
                combustor_network(gas; burned_state=(temperature=expected[end,1],mass_fractions=Y))
            else
                closed_pair_network(gas)
            end
            rhs = network_rhs(network)
            solution = solve_network(network,(first(times),last(times));integrator=native_network_bdf,
                initial_state=expected[:,1],reltol=1e-11,abstol=1e-18,
                saveat=times,tstops=times,maxiters=1_000_000)
            states = reduce(hcat,solution.u)
            diagnostics = [network_diagnostics(rhs,state,time) for (state,time) in zip(solution.u,solution.t)]
            temperature_error, composition_error, mass_error, pressure_error = 0.0,0.0,0.0,0.0
            minimum_fraction = 0.0
            total_mass, total_energy = zeros(length(times)),zeros(length(times))
            reactor_index = 0
            for (node_index,node) in enumerate(network.nodes)
                offset = network.offsets[node_index]
                offset == 0 && continue
                reactor_index += 1
                ns = node.initial.gas.n_species
                native_masses = states[offset:offset+ns-1,:]
                reference_masses = expected[offset:offset+ns-1,:]
                native_mass, reference_mass = sum(native_masses;dims=1),sum(reference_masses;dims=1)
                native_Y, reference_Y = native_masses ./ native_mass,reference_masses ./ reference_mass
                temperature_error = max(temperature_error,maximum(abs,states[offset+ns,:]-expected[offset+ns,:]))
                composition_error = max(composition_error,maximum(abs,native_Y-reference_Y))
                mass_error = max(mass_error,maximum(abs,native_mass ./ reference_mass .- 1))
                node_pressures = [d.nodes[node_index].pressure for d in diagnostics]
                pressure_error = max(pressure_error,maximum(abs,node_pressures ./ reference["pressure"][reactor_index,:] .- 1))
                minimum_fraction = min(minimum_fraction,minimum(native_Y))
                total_mass .+= vec(native_mass)
                total_energy .+= [d.nodes[node_index].total_internal_energy for d in diagnostics]
            end
            native_flows = reduce(hcat,[d.mass_flow_rates for d in diagnostics])
            flow_error = maximum(abs,native_flows-reference["mass_flow"]) / max(maximum(abs,reference["mass_flow"]),1e-3)
            mass_balance_error = maximum(abs(d.total_mass_rate-d.external_mass_rate) for d in diagnostics)
            energy_balance_error = maximum(abs(d.total_energy_rate-d.external_energy_rate) for d in diagnostics)
            report = Dict{String,Any}("case"=>name,"cantera_version"=>String(vec(reference["cantera_version_utf8"])),
                "temperature_error_K"=>temperature_error,"mass_fraction_error"=>composition_error,
                "mass_relative_error"=>mass_error,"pressure_relative_error"=>pressure_error,
                "flow_relative_error"=>flow_error,"mass_balance_error_kg_s"=>mass_balance_error,
                "energy_balance_error_W"=>energy_balance_error,"minimum_mass_fraction"=>minimum_fraction,
                "saved_points"=>length(times))
            @testset "$name" begin
                @test size(states) == size(expected)
                @test temperature_error < 0.05
                @test composition_error < 5e-7
                @test mass_error < 1e-6
                @test pressure_error < 1e-6
                @test flow_error < 1e-4
                @test mass_balance_error < 1e-10
                @test energy_balance_error < 1e-5
                @test minimum_fraction > -2e-12
                if name == "closed_pair"
                    mass_drift = maximum(abs,total_mass .- first(total_mass)) / first(total_mass)
                    energy_drift = maximum(abs,total_energy .- first(total_energy)) / abs(first(total_energy))
                    wall_error = maximum(abs,reduce(hcat,[d.wall_heat_rates for d in diagnostics])-reference["wall_heat"])
                    @test mass_drift < 1e-9
                    @test energy_drift < 1e-8
                    @test wall_error < 0.05
                    report["closed_mass_relative_drift"] = mass_drift
                    report["closed_energy_relative_drift"] = energy_drift
                    report["wall_heat_error_W"] = wall_error
                else
                    steady = solve_network_steady(network;integrator=native_network_bdf,
                        initial_state=solution.u[end],interval=0.2,max_time=3.0,steady_tolerance=1e-7,
                        reltol=1e-11,abstol=1e-18,save_everystep=false,save_start=false,maxiters=1_000_000)
                    reference_steady = vec(reference["steady_state"])
                    @test abs(steady.state[end]-reference_steady[end]) < 0.01
                    @test maximum(abs,steady.state[1:end-1]-reference_steady[1:end-1]) < 1e-7
                    report["steady_temperature_K"] = steady.state[end]
                    report["steady_residual_per_s"] = steady.residual
                end
            end
            push!(reports,report)
            npzwrite(joinpath(directory,"native_network_$name.npz"),Dict("time"=>solution.t,"state"=>states))
            @printf("%-12s ΔT=%9.3g K  ΔY=%9.3g  mass balance=%9.3g kg/s  energy balance=%9.3g W\n",
                name,temperature_error,composition_error,mass_balance_error,energy_balance_error)
        end
    end
    open(joinpath(directory,"native_network_report.toml"),"w") do io
        TOML.print(io,Dict("julia_version"=>string(VERSION),"solver"=>"OrdinaryDiffEqBDF.QNDF",
            "relative_tolerance"=>1e-11,"absolute_tolerance"=>1e-18,"cases"=>reports))
    end
    return reports
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: reactor_network_cases.jl REFERENCE_DIRECTORY")
    validate_network_cases(ARGS[1])
end
