using Arrhenius
using NPZ
using Test
using TOML

include(joinpath(@__DIR__,"..","example","reactors","catalytic_reactor.jl"))
include(joinpath(@__DIR__,"..","example","reactors","diamond_cvd.jl"))

function validate_surface_reactors(directory)
    reports = Dict{String,Any}()
    for case in ("closed","flow")
        system = catalytic_reactor(joinpath(directory,"pt_h2.surface.npz");flow=case=="flow")
        ref = npzread(joinpath(directory,"surface_reactor_$case.reference.npz"))
        times = vec(ref["times"])
        solution = solve_catalytic(system,(0.0,last(times));integrator=native_surface_bdf,
            reltol=1e-10,abstol=1e-19,saveat=times,tstops=times[2:end],maxiters=1_000_000)
        u = reduce(hcat,solution.u)
        ns = system.network.nodes.reactor.initial.gas.n_species
        mass = vec(sum(u[1:ns,:];dims=1))
        Y = u[1:ns,:]./transpose(mass)
        theta = u[system.offsets[1]:end,:]
        reports[case] = Dict(
            "mass_relative_error"=>maximum(abs.(mass-ref["mass"])./ref["mass"]),
            "mass_fraction_absolute_error"=>maximum(abs,Y-ref["mass_fractions"]),
            "coverage_absolute_error"=>maximum(abs,theta-ref["coverages"]),
            "site_balance_drift"=>maximum(abs,vec(sum(theta;dims=1)).-1))
        @test reports[case]["mass_relative_error"] < 1e-7
        @test reports[case]["mass_fraction_absolute_error"] < 2e-7
        @test reports[case]["coverage_absolute_error"] < 2e-7
        @test reports[case]["site_balance_drift"] < 1e-8
        rhs = catalytic_rhs(system)
        initial = catalytic_diagnostics(rhs,solution.u[1])
        if case == "closed"
            diagnostics = [catalytic_diagnostics(rhs,v,t) for (v,t) in zip(solution.u,times)]
            mass_drift = maximum(abs(d.mass-initial.mass)/initial.mass for d in diagnostics)
            element_drift = maximum(maximum(abs,d.element_inventory-initial.element_inventory) for d in diagnostics)
            reports[case]["total_mass_drift"] = mass_drift
            reports[case]["element_inventory_drift"] = element_drift
            @test mass_drift < 1e-8
            @test element_drift < 1e-17
            @test maximum(abs(d.mass_rate) for d in diagnostics) < 1e-15
            @test maximum(maximum(abs,d.element_rates) for d in diagnostics) < 1e-15
        end
        @show case reports[case]
        npzwrite(joinpath(directory,"surface_reactor_$case.native.npz"),Dict("state"=>u,"times"=>times))
    end
    ref = npzread(joinpath(directory,"diamond_continuation.reference.npz"))
    result = diamond_cvd(joinpath(directory,"diamond.surface.npz"))
    growth = [r.growth for r in result]
    hydrogen = [r.hydrogen for r in result]
    coverage = reduce(hcat,[r.coverages for r in result])
    growth_error = maximum(abs.((growth-ref["growth"])./ref["growth"]))
    reports["diamond_continuation"] = Dict("growth_relative_error"=>growth_error,
        "coverage_absolute_error"=>maximum(abs,coverage-ref["coverages"]))
    @test hydrogen ≈ ref["hydrogen"] rtol=1e-13
    @test growth_error < 1e-5
    @test maximum(abs,coverage-ref["coverages"]) < 2e-7
    open(joinpath(directory,"surface_reactor_validation.toml"),"w") do io
        TOML.print(io,reports;sorted=true)
    end
    return reports
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("supply the surface reference directory")
    @testset "Cantera coupled surface reactors" begin
        validate_surface_reactors(ARGS[1])
    end
end
