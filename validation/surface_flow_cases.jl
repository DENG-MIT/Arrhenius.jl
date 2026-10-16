using Arrhenius
using NPZ
using Test
using TOML

include(joinpath(@__DIR__,"..","example","reactors","surf_pfr.jl"))
include(joinpath(@__DIR__,"..","example","reactors","surf_pfr_chain.jl"))

function validate_surface_flows(directory;run=true)
    path = joinpath(directory,"pox.surface.npz")
    flow = methane_surface_pfr(path)
    if run
        direct = surf_pfr(path;output=joinpath(directory,"flow.native.csv"))
        npzwrite(joinpath(directory,"flow.native.npz"),Dict("state"=>reduce(hcat,direct.solution.u),
            "mole_fractions"=>direct.mole_fractions,"pressure"=>direct.pressure,"distance"=>direct.solution.t))
        chain = surf_pfr_chain(path;output=joinpath(directory,"chain.native.csv"))
        npzwrite(joinpath(directory,"chain.native.npz"),Dict("mole_fractions"=>chain.mole_fractions,
            "mass_fractions"=>chain.mass_fractions,"coverages"=>chain.coverages,
            "pressure"=>chain.pressure,"distance"=>chain.distance,"residuals"=>chain.residuals))
    end
    reports = Dict{String,Any}()
    inlet_elements = flow.gas.ele_matrix*(flow.mass_fractions./flow.gas.MW)
    for case in ("flow","chain")
        ref = npzread(joinpath(directory,case=="flow" ? "flow.reference.npz" : "chain.refined.npz"))
        native = npzread(joinpath(directory,"$case.native.npz"))
        ng = flow.gas.n_species
        Y = case=="flow" ? native["state"][1:ng,:] : native["mass_fractions"]
        theta = case=="flow" ? native["state"][ng+1:end,:] : native["coverages"]
        elements = flow.gas.ele_matrix*(Y./flow.gas.MW)
        reference_elements = flow.gas.ele_matrix*(ref["mass_fractions"]./flow.gas.MW)
        reports[case] = Dict(
            "mole_fraction_absolute_error"=>maximum(abs,native["mole_fractions"]-ref["mole_fractions"]),
            "coverage_absolute_error"=>maximum(abs,theta-ref["coverages"]),
            "pressure_relative_error"=>maximum(abs.(native["pressure"]-ref["pressure"])./ref["pressure"]),
            "element_inventory_drift"=>maximum(abs,elements.-inlet_elements),
            "reference_element_inventory_drift"=>maximum(abs,reference_elements.-inlet_elements),
            "mass_flux_relative_drift"=>maximum(abs,vec(sum(Y;dims=1)).-1),
            "site_balance_drift"=>maximum(abs,vec(sum(theta;dims=1)).-1))
        @test native["distance"] ≈ ref["distance"] rtol=1e-14
        @test reports[case]["mole_fraction_absolute_error"] < 2e-7
        @test reports[case]["coverage_absolute_error"] < 2e-7
        @test reports[case]["pressure_relative_error"] < 1e-8
        @test reports[case]["element_inventory_drift"] < 1e-9
        @test reports[case]["mass_flux_relative_drift"] < 1e-9
        @test reports[case]["site_balance_drift"] < 1e-9
        @test minimum(Y) >= -1e-12
        @test minimum(theta) >= -1e-12
        if case == "flow"
            diagnostics = [surface_flow_diagnostics(flow,u) for u in eachcol(native["state"])]
            stationarity = maximum(d.stationary_relative_residual for d in diagnostics)
            momentum_error = maximum(abs(d.pressure+d.density*d.speed^2-flow.momentum_flux)/flow.momentum_flux for d in diagnostics)
            reports[case]["stationary_relative_residual"] = stationarity
            reports[case]["momentum_flux_relative_error"] = momentum_error
            @test stationarity < 1e-7
            @test momentum_error < 1e-13
            @test [d.density for d in diagnostics] ≈ ref["density"] rtol=1e-7
            @test [d.speed for d in diagnostics] ≈ ref["speed"] rtol=1e-7
            npzwrite(joinpath(directory,"flow.fluxes.npz"),Dict(
                "species_mass_flux"=>reduce(hcat,[d.species_mass_flux for d in diagnostics]),
                "elemental_flux"=>reduce(hcat,[d.elemental_flux for d in diagnostics])))
        else
            source = npzread(joinpath(directory,"chain.reference.npz"))
            reports[case]["source_default_mole_fraction_difference"] = maximum(abs,native["mole_fractions"]-source["mole_fractions"])
            reports[case]["source_default_element_inventory_drift"] = maximum(abs,
                flow.gas.ele_matrix*(source["mass_fractions"]./flow.gas.MW).-inlet_elements)
            @test maximum(native["residuals"]) < 1e-9
            npzwrite(joinpath(directory,"chain.fluxes.npz"),Dict(
                "species_mass_flux"=>flow.mass_flux.*Y,"elemental_flux"=>flow.mass_flux.*elements))
        end
        @show case reports[case]
    end
    open(joinpath(directory,"surface_flow_validation.toml"),"w") do io
        TOML.print(io,reports;sorted=true)
    end
    return reports
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("supply the prepared surface-flow reference directory")
    @testset "Cantera surface flow and 201-reactor chain" begin
        validate_surface_flows(ARGS[1])
    end
end
