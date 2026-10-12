using Arrhenius
using NPZ
using Test
using TOML

include(joinpath(@__DIR__, "..", "example", "reactors", "surface_solver.jl"))

function validate_surfaces(directory)
    reports = Dict{String,Any}()
    for name in ("pt_h2", "diamond")
        m = SurfaceMechanism(joinpath(directory, "$name.surface.npz"))
        ref = npzread(joinpath(directory, "$name.reference.npz"))
        errors = Dict{String,Float64}()
        for (key, getter) in (("forward", w->w.forward_rates), ("reverse", w->w.reverse_rates),
                              ("forward_constants", w->w.forward_rate_constants),
                              ("reverse_constants", w->w.reverse_rate_constants),
                              ("production", w->w.production_rates))
            predicted = similar(ref[key])
            for j in eachindex(ref["temperatures"])
                s = IdealSurface(m; temperature=ref["temperatures"][j],
                    pressure=only(ref["pressure"]), gas_temperature=ref["gas_temperatures"][j],
                    mole_fractions=ref["mole_fractions"], coverages=ref["coverages"][:,j])
                predicted[:,j] = getter(surface_rates!(SurfaceWorkspace(m), s))
            end
            scale = max.(abs.(ref[key]), 1e-30)
            errors[key*"_relative_error"] = maximum(abs.(predicted-ref[key])./scale)
            @test isapprox(predicted, ref[key]; rtol=2e-11, atol=1e-25)
        end
        s = IdealSurface(m; temperature=only(ref["temperature"]), pressure=only(ref["pressure"]),
                         mole_fractions=ref["mole_fractions"], coverages=ref["initial_coverages"])
        times = vec(ref["times"])
        solution = solve_surface(s, (0.0,last(times)); integrator=native_surface_bdf,
            reltol=1e-10, abstol=1e-18, saveat=times, tstops=times[2:end], maxiters=1_000_000)
        trajectory = reduce(hcat, solution.u)
        errors["trajectory_absolute_error"] = maximum(abs, trajectory-ref["trajectory"])
        errors["site_balance_drift"] = maximum(abs, vec(sum(trajectory; dims=1)).-1)
        @test errors["trajectory_absolute_error"] < 2e-7
        @test errors["site_balance_drift"] < 1e-9
        @test minimum(trajectory) >= -1e-13
        stationary = steady_coverages(s; integrator=native_surface_bdf,
            reltol=1e-11, abstol=1e-19, maxiters=1_000_000, initial_coverages=solution.u[end])
        errors["steady_coverage_error"] = maximum(abs, stationary.coverages-ref["steady_coverages"])
        errors["steady_residual_per_second"] = stationary.residual
        @test errors["steady_coverage_error"] < 2e-7
        @test abs(sum(stationary.coverages)-1) < 1e-9
        rates = surface_rates(s, stationary.coverages)
        @test maximum(abs, rates.elemental_rates) < 1e-13
        @test isapprox(vcat(rates.surface,rates.gas,rates.bulk), ref["steady_production"]; rtol=1e-5, atol=1e-13)
        reports[name] = errors
        @show name errors
        npzwrite(joinpath(directory, "$name.native.npz"), Dict("trajectory"=>trajectory,
            "steady_coverages"=>stationary.coverages, "times"=>times))
    end
    open(joinpath(directory, "surface_validation.toml"), "w") do io
        TOML.print(io, reports; sorted=true)
    end
    return reports
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: julia surface_cases.jl REFERENCE_DIRECTORY")
    @testset "Cantera surface kinetics and coverage references" begin
        validate_surfaces(ARGS[1])
    end
end
