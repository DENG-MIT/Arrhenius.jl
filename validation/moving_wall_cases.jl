using Arrhenius
using NPZ
using Test
using TOML
include(joinpath(@__DIR__,"..","example","reactors","network_solver.jl"))
include(joinpath(@__DIR__,"..","example","reactors","moving_wall_setup.jl"))

function validate_moving_wall_cases(directory;run=true)
    metrics = Dict{String,Any}()
    @testset "Cantera 4 moving-wall source trajectories" begin
        for name in ("reactor2","piston","custom2")
            case = moving_wall_example(name,joinpath(directory,"mechanisms"))
            rhs = moving_wall_rhs(case.model)
            reference = npzread(joinpath(directory,name*".reference.npz"))
            output = joinpath(directory,name*".native.npz")
            if run
                result = run_moving_wall_example(name,joinpath(directory,"mechanisms"),
                    joinpath(directory,name*".native.csv");integrator=native_network_bdf)
                npzwrite(output,Dict("time"=>result.solution.t,"states"=>reduce(hcat,result.solution.u)))
            end
            native = npzread(output)
            times,states = native["time"],native["states"]
            @test times ≈ reference["time"] rtol=0 atol=1e-14
            temperature_error,pressure_error,volume_error,mass_fraction_error = 0.0,0.0,0.0,0.0
            energy_drift,mass_drift,element_drift,volume_drift,mechanical_drift = 0.0,0.0,0.0,0.0,0.0
            velocity_error = 0.0
            initial = moving_wall_diagnostics(rhs,@view(states[:,1]),times[1])
            scale = max(abs(initial.total_internal_energy),1.0)
            initial_mechanical = initial.total_internal_energy+initial.wall_kinetic_energy
            if name == "custom2"
                initial_mechanical += one_atm*initial.nodes.reactor.volume
            end
            for j in eachindex(times)
                d = moving_wall_diagnostics(rhs,@view(states[:,j]),times[j])
                vessel_index = 0
                for (i,node) in enumerate(case.model.network.nodes)
                    node isa WellStirredReactor || continue
                    vessel_index += 1
                    s = d.nodes[i]
                    temperature_error = max(temperature_error,abs(s.temperature-reference["temperature"][j,vessel_index]))
                    pressure_error = max(pressure_error,abs(s.pressure/reference["pressure"][j,vessel_index]-1))
                    volume_error = max(volume_error,abs(s.volume-reference["volume"][j,vessel_index]))
                    mass_fraction_error = max(mass_fraction_error,
                        maximum(abs,s.mass_fractions-reference["Y$vessel_index"][j,:]))
                end
                # CT may report the left or right limit for a discontinuous callback at t=0.1.
                if name != "piston" || abs(times[j]-0.1)>1e-12
                    velocity_error = max(velocity_error,abs(d.wall_velocities[1]-reference["velocity"][j,1]))
                end
                energy_drift = max(energy_drift,abs(d.energy_balance-initial.energy_balance)/scale)
                mass_drift = max(mass_drift,abs(d.mass_balance/initial.mass_balance-1))
                element_drift = max(element_drift,maximum(abs(d.element_inventories[k]-v)
                                    for (k,v) in initial.element_inventories))
                if name != "custom2"
                    volume_drift = max(volume_drift,abs(d.total_volume-initial.total_volume))
                else
                    mechanical = d.total_internal_energy+d.wall_kinetic_energy+one_atm*d.nodes.reactor.volume
                    mechanical_drift = max(mechanical_drift,abs(mechanical-initial_mechanical)/scale)
                end
            end
            @test temperature_error < 0.05
            @test pressure_error < 2e-5
            @test volume_error < 2e-6
            @test mass_fraction_error < 2e-6
            @test velocity_error < 2e-3
            @test energy_drift < 2e-7
            @test mass_drift < 1e-10
            @test element_drift < 1e-11
            @test volume_drift < 1e-10
            @test mechanical_drift < 2e-7
            metrics[name] = Dict("max_temperature_error_K"=>temperature_error,"max_pressure_relative_error"=>pressure_error,
                "max_volume_error_m3"=>volume_error,"max_mass_fraction_error"=>mass_fraction_error,
                "max_velocity_error_m_per_s"=>velocity_error,"energy_balance_relative_drift"=>energy_drift,
                "mass_balance_relative_drift"=>mass_drift,"element_inventory_drift_kmol"=>element_drift,
                "volume_sum_drift_m3"=>volume_drift,"inertial_energy_relative_drift"=>mechanical_drift,
                "points"=>length(times))
            println(name,": ",metrics[name])
        end
    end
    @testset "analytic moving-wall integration identities" begin
        gas = CreateSolution(joinpath(directory,"mechanisms","air.yaml"))
        chamber = WellStirredReactor(gas;temperature=900.0,pressure=2one_atm,
            mole_fractions=Dict("AR"=>1.0),volume=0.3,chemistry=false)
        ambient = Reservoir(gas;temperature=300.0,pressure=one_atm,mole_fractions=Dict("AR"=>1.0))
        model = MovingWallNetwork(ReactorNetwork((chamber=chamber,ambient=ambient));
            walls=(MovingWall(:chamber,:ambient;area=0.2,velocity=t->0.1cos(2t)),))
        times = collect(range(0,1.0;length=51))
        sol = solve_moving_wall(model,(0.0,1.0);integrator=native_network_bdf,
            reltol=1e-11,abstol=1e-18,saveat=times)
        rhs = moving_wall_rhs(model)
        volumes,temperatures,pressures,balances = Float64[],Float64[],Float64[],Float64[]
        for (t,u) in zip(sol.t,sol.u)
            d = moving_wall_diagnostics(rhs,u,t)
            push!(volumes,d.nodes.chamber.volume)
            push!(temperatures,d.nodes.chamber.temperature)
            push!(pressures,d.nodes.chamber.pressure)
            push!(balances,d.energy_balance)
        end
        analytic_volume = @. 0.3+0.01sin(2times)
        @test maximum(abs,volumes-analytic_volume) < 1e-9
        @test maximum(abs,temperatures-900 .* (0.3 ./ analytic_volume).^(2/3)) < 1e-5
        @test maximum(abs,pressures ./ (2one_atm .* (0.3 ./ analytic_volume).^(5/3)).-1) < 1e-8
        @test maximum(abs,balances.-balances[1]) < 1e-4
    end
    open(joinpath(directory,"moving_wall_validation.toml"),"w") do io
        TOML.print(io,metrics;sorted=true)
    end
    return metrics
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("provide the Cantera 4 moving-wall reference directory")
    validate_moving_wall_cases(ARGS[1];run=!("--cached" in ARGS))
end
