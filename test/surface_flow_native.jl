using Arrhenius
using Test
using NPZ
using LinearAlgebra

@testset "Native stationary surface flow equations" begin
    mktempdir() do directory
        gas = CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
        ng = gas.n_species
        kH,kH2 = findfirst(==("H"),gas.species_names),findfirst(==("H2"),gas.species_names)
        names = vcat(["PT"],gas.species_names)
        elements = vcat(["Pt"],gas.elements)
        elemental = zeros(length(elements),ng+1)
        elemental[1,1] = 1
        elemental[2:end,2:end] = gas.ele_matrix
        reactants,products = zeros(ng+1,1),zeros(ng+1,1)
        reactants[1,1],reactants[kH2+1,1] = 1,1
        products[1,1],products[kH+1,1] = 1,2
        coefficients = zeros(ng+1,15); coefficients[:,1] .= 298.15
        b(s) = collect(codeunits(s))
        path = joinpath(directory,"surface.npz")
        npzwrite(path,Dict{String,Any}(
            "surface_format_utf8"=>b("arrhenius-surface-v1"),"phase_name_utf8"=>b("toy"),
            "species_names_utf8"=>b(join(names,'\n')),"element_names_utf8"=>b(join(elements,'\n')),
            "n_surface"=>[1],"n_gas"=>[ng],"site_density"=>[2e-5],"site_sizes"=>[1.0],
            "molecular_weights"=>vcat([195.084],gas.MW),"elemental_matrix"=>elemental,
            "reactants"=>reactants,"products"=>products,"orders"=>reactants,
            "arrhenius"=>reshape([1.0,0.0,0.0],1,3),"reversible"=>[false],
            "coverage_a"=>zeros(1,1),"coverage_m"=>zeros(1,1),"coverage_energy"=>zeros(1,1,4),
            "sticking_species"=>[0],"sticking_order"=>[0.0],"sticking_factor"=>[0.0],
            "motz_wise"=>[false],"thermo_type"=>fill(2,ng+1),"thermo_coefficients"=>coefficients,
            "reference_pressure"=>fill(one_atm,ng+1),"initial_coverages"=>[1.0],
            "initial_mole_fractions"=>ones(ng)./ng,"initial_temperature"=>[900.0],
            "initial_pressure"=>[one_atm]))
        m = SurfaceMechanism(path)
        flow = SurfaceFlowReactor(m;gas,temperature=900,mole_fractions=ones(ng),
            area=0.05,mass_flow_rate=0.01,surface_area_per_length=2.0)
        props = surface_flow_properties(flow,flow.mass_fractions)
        @test props.pressure ≈ one_atm rtol=1e-14
        @test props.density*props.speed ≈ flow.mass_flux
        @test props.pressure+props.density*props.speed^2 ≈ flow.momentum_flux
        p = surface_flow_problem(flow,0.1;initial_coverages=[1.0])
        @test diag(p.mass_matrix) == vcat(ones(ng),[0.0])
        f = similar(p.u0);p.f(f,p.u0,nothing,0.0)
        @test abs(sum(f[1:ng])) < 1e-10*maximum(abs,f)
        @test f[end] == 0
        @test maximum(abs,gas.ele_matrix*(f[1:ng]./gas.MW)) < 1e-10*maximum(abs,f)
        @test !p.isoutofdomain(p.u0,nothing,0)
        @test p.isoutofdomain(fill(-1.,ng+1),nothing,0)
        J = zeros(ng+1,ng+1)
        surface_flow_jacobian!(J,p.u0,p.f)
        direction = p.u0.*sin.(0.5.*eachindex(p.u0))
        fp,fm = similar(f),similar(f);h=1e-5
        p.f(fp,p.u0+h*direction,nothing,0)
        p.f(fm,p.u0-h*direction,nothing,0)
        @test J*direction ≈ (fp-fm)/(2h) rtol=2e-6
        diagnostics = surface_flow_diagnostics(flow,p.u0)
        @test sum(diagnostics.species_mass_flux) ≈ flow.mass_flux
        @test diagnostics.coverage_sum == 1
        @test diagnostics.stationary_relative_residual == 0
        @test_throws ArgumentError SurfaceFlowReactor(m;gas,temperature=900,mole_fractions=ones(ng),
            area=0.05,mass_flow_rate=-0.01,surface_area_per_length=2)
        @test_throws ArgumentError SurfaceFlowReactor(m;gas,temperature=900,mole_fractions=ones(ng),
            area=0.05,mass_flow_rate=1e5,surface_area_per_length=2)
        @test_throws DimensionMismatch surface_flow_properties(flow,[1.0])
        @test_throws ArgumentError surface_flow_problem(flow,-1;initial_coverages=[1.0])
        @test_throws ArgumentError surface_reactor_chain(flow,0.1;reactors=1,integrator=identity)
        push!(m.bulk_molar_volumes,1.0)
        @test_throws ArgumentError SurfaceFlowReactor(m;gas,temperature=900,mole_fractions=ones(ng),
            area=0.05,mass_flow_rate=0.01,surface_area_per_length=2)
    end
end
