# Usage: julia --project=SOLVER_ENV real_gas_trajectories.jl REFERENCE_DIR [FILENAME_REGEX] [--saved]
# Requires SciMLBase and OrdinaryDiffEqSDIRK. Prepare references with
# real_gas_reactor_cases.py REFERENCE_DIR --sweep --refine.
using Arrhenius, NPZ, LinearAlgebra, Test, Printf, TOML
BLAS.set_num_threads(1)
for (name,file) in ((:RedlichKwongThermo,"RealGasThermo.jl"),(:RedlichKwongReactor,"RealGasReactors.jl"))
    isdefined(Arrhenius,name) || Base.include(Arrhenius,joinpath(@__DIR__,"..","src",file))
end
include(joinpath(@__DIR__,"..","example","reactors","real_gas_ode_solver.jl"))
directory = ARGS[1]
gas = CreateSolution(joinpath(directory,"dodecane_IG.yaml"))
model = RedlichKwongThermo(joinpath(directory,"dodecane_RK.yaml"))
files = sort(filter(f -> startswith(f,"trajectory_") && endswith(f,".npz"),readdir(directory)))
patterns = filter(!=("--saved"),ARGS[2:end])
!isempty(patterns) && filter!(f -> occursin(Regex(only(patterns)),f),files)
isempty(files) && error("no matching shock-tube references found")
reports = Dict{String,Any}[]
@testset "Native nonideal shock-tube trajectories" begin
    for file in files
        data = npzread(joinpath(directory,file))
        haskey(data,"refined_state") || error("prepare converged references with real_gas_reactor_cases.py --refine")
        expected,times = data["refined_state"],data["time"]
        temperature,pressure,Y = expected[end,1],data["refined_P"][1],expected[1:end-1,1]
        if occursin("_RK_",file)
            reactor = RedlichKwongReactor(gas,model;temperature,pressure,mass_fractions=Y)
        else
            reactor = IdealGasReactor(gas;temperature,pressure,mass_fractions=Y,constraint=:constant_volume)
        end
        if "--saved" in ARGS
            saved = npzread(joinpath(directory,"native_"*file))
            saved["time"] == times || error("saved native sample times do not match")
            solution = (;u=collect(eachcol(saved["state"])),t=times)
            elapsed = NaN
        else
            println("Integrating ",file," at ",length(times)," source sample times"); flush(stdout)
            elapsed = @elapsed solution = solve_reactor(reactor,(0.,max(.005,last(times)));integrator=shocktube_sdirk,
                reltol=1e-13,abstol=1e-26,saveat=times,tstops=times,save_end=last(times)>=.005,maxiters=1_000_000)
        end
        states = reduce(hcat,solution.u)
        properties = [reactor_properties(reactor,u) for u in solution.u]
        initial = first(properties)
        dT = maximum(abs,states[end,:]-expected[end,:])
        dY = maximum(abs,states[1:end-1,:]-expected[1:end-1,:])
        dP = maximum(abs.([p.pressure for p in properties]./data["refined_P"].-1))
        dm = maximum(abs(p.mass_fraction_sum-1) for p in properties)
        de = maximum(norm(p.elemental_inventory-initial.elemental_inventory,Inf) for p in properties)
        denergy = maximum(abs(p.internal_energy-initial.internal_energy) for p in properties)/max(abs(initial.internal_energy),1e6)
        rhs = reactor_rhs(reactor)
        derivative = zero(first(solution.u))
        species_error,temp_error = 0.,0.
        for i in axes(expected,2)
            rhs(derivative,view(expected,:,i),nothing,times[i])
            species_error = max(species_error,maximum(abs,derivative[1:end-1]-data["refined_rhs"][1:end-1,i]))
            temp_error = max(temp_error,abs(derivative[end]-data["refined_rhs"][end,i]))
        end
        rhs_error = max(species_error/max(maximum(abs,data["refined_rhs"][1:end-1,:]),1),temp_error/max(maximum(abs,data["refined_rhs"][end,:]),1))
        oh = only(data["oh_index"])
        delay = times[argmax(states[oh,:])]
        reference_delay = times[argmax(expected[oh,:])]
        @testset "$file" begin
            @test dT < 0.3
            @test dY < 2e-5
            @test dP < 1e-4
            @test dm < 5e-10
            @test de < 5e-11
            @test denergy < 1e-6
            @test minimum(states[1:end-1,:]) > -1e-12
            @test rhs_error < 1e-8
            @test delay == reference_delay
        end
        @printf("%s: ΔT=%.4g K, ΔY=%.4g, energy drift=%.3g, RHS=%.3g, source-grid delay=%.9g s, elapsed incl JIT=%.3f s\n",file,dT,dY,denergy,rhs_error,delay,elapsed)
        flush(stdout)
        npzwrite(joinpath(directory,"native_"*file),Dict("state"=>states,"time"=>times))
        push!(reports,Dict("case"=>file,"temperature_error_K"=>dT,"mass_fraction_error"=>dY,
            "pressure_relative_error"=>dP,"mass_drift"=>dm,"element_drift"=>de,"energy_relative_drift"=>denergy,
            "rhs_relative_error"=>rhs_error,"source_sample_ignition_delay_s"=>delay,
            "published_ignition_delay_s"=>only(data["ignition_delay"]),
            "source_default_temperature_difference_K"=>maximum(abs,states[end,:]-data["state"][end,:]),
            "reference_rtol"=>only(data["reference_rtol"]),"reference_atol"=>only(data["reference_atol"]),
            "native_solver"=>"OrdinaryDiffEqSDIRK.KenCarp4","native_rtol"=>1e-13,"native_atol"=>1e-26,
            "native_end_time_s"=>max(.005,last(times)),
            "elapsed_including_JIT_s"=>elapsed,"qualification"=>"correctness only; not a controlled WSL benchmark"))
    end
end
open(joinpath(directory,"native_shocktube_report.toml"),"w") do io
    TOML.print(io,Dict("cases"=>reports))
end
