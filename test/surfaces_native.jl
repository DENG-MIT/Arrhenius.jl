using Test
using NPZ
using Arrhenius

# H(g) + PT(s) <=> H(s), with analytic constant-cp reference thermodynamics.
function simple_surface_fixture(path; sticking=false, reversible=true, site_size=1.0,
                                coverage=false, motz_wise=false)
    textbytes(s) = collect(codeunits(s))
    cov_a, cov_m, cov_e = zeros(1,2), zeros(1,2), zeros(1,2,4)
    if coverage
        cov_a[1,2], cov_m[1,2] = 0.2, 0.3
        cov_e[1,2,:] = [2e6, -3e6, 4e6, -1e6]
    end
    coeffs = zeros(3,15)
    coeffs[:,1] .= 298.15
    payload = Dict{String,Any}(
        "surface_format_utf8"=>textbytes("arrhenius-surface-v1"),
        "phase_name_utf8"=>textbytes("simple"),
        "species_names_utf8"=>textbytes("PT\nHS\nH"),
        "element_names_utf8"=>textbytes("Pt\nH"), "gas_file_utf8"=>textbytes("unused.yaml"),
        "n_surface"=>[2],"n_gas"=>[1], "site_density"=>[2e-5],
        "site_sizes"=>[site_size,site_size], "molecular_weights"=>[195.084,196.092,1.008],
        "elemental_matrix"=>[1.0 1.0 0.0; 0.0 1.0 1.0],
        "reactants"=>reshape([1.0,0.0,1.0],3,1),
        "products"=>reshape([0.0,1.0,0.0],3,1),
        "orders"=>reshape([1.0,0.0,1.0],3,1),
        "arrhenius"=>reshape([sticking ? 0.2 : 2.0,0.0,0.0],1,3),
        "reversible"=>[reversible], "coverage_a"=>cov_a,"coverage_m"=>cov_m,
        "coverage_energy"=>cov_e,"sticking_species"=>[sticking ? 3 : 0],
        "sticking_order"=>[1.0],
        "sticking_factor"=>[sqrt(Arrhenius.R/(2pi*1.008))*site_size],
        "motz_wise"=>[motz_wise],"thermo_type"=>[2,2,2],
        "thermo_coefficients"=>coeffs,"reference_pressure"=>fill(101325.0,3),
        "initial_coverages"=>[0.7,0.3],"initial_mole_fractions"=>[1.0],
        "initial_temperature"=>[900.0],"initial_pressure"=>[101325.0])
    npzwrite(path,payload)
    return SurfaceMechanism(path)
end

@testset "Native gas and surface coupling" begin
    mktempdir() do directory
        path = joinpath(directory,"coupled.npz")
        simple_surface_fixture(path;reversible=false)
        data = npzread(path)
        gas = CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
        ng = gas.n_species
        kH = findfirst(==("H"),gas.species_names)
        data["n_gas"] = [ng]
        data["species_names_utf8"] = collect(codeunits(join(vcat(["PT","HS"],gas.species_names),'\n')))
        elements = vcat(["Pt"],gas.elements)
        data["element_names_utf8"] = collect(codeunits(join(elements,'\n')))
        elemental = zeros(length(elements),2+ng)
        elemental[1,1:2] .= 1
        elemental[findfirst(==("H"),elements),2] = 1
        elemental[2:end,3:end] = gas.ele_matrix
        data["elemental_matrix"] = elemental
        data["molecular_weights"] = vcat([195.084,195.084+gas.MW[kH]],gas.MW)
        reactants,products = zeros(2+ng,1),zeros(2+ng,1)
        reactants[1,1],reactants[kH+2,1],products[2,1] = 1,1,1
        data["reactants"],data["products"],data["orders"] = reactants,products,reactants
        data["thermo_type"] = fill(2,2+ng)
        coefficients = zeros(2+ng,15)
        coefficients[:,1] .= 298.15
        data["thermo_coefficients"] = coefficients
        data["reference_pressure"] = fill(one_atm,2+ng)
        data["initial_mole_fractions"] = ones(ng)./ng
        npzwrite(path,data)
        m = SurfaceMechanism(path)
        initial = IdealGasReactor(gas;temperature=900,mole_fractions=ones(ng),
            constraint=:constant_volume,energy=:isothermal)
        vessel = WellStirredReactor(initial;volume=0.1,chemistry=false)
        network = ReactorNetwork((reactor=vessel,))
        surface = ReactorSurface(:reactor,m;area=0.2)
        system = CatalyticNetwork(network;surfaces=(surface,))
        u = catalytic_state(system)
        rhs = catalytic_rhs(system)
        report = catalytic_diagnostics(rhs,u)
        @test abs(report.mass_rate) < 1e-18
        @test maximum(abs,report.element_rates) < 1e-18
        @test report.derivative[kH] < 0
        @test report.derivative[end] > 0
        @test report.derivative[ng+1] == 0
        @test report.boundary_power == 0
        @test report.thermostat_power == report.internal_energy_rate
        @test report.mass > initial.density*vessel.volume
        before = copy(u)
        J = zeros(length(u),length(u))
        catalytic_jacobian!(J,u,rhs)
        direction = u.*sin.(0.7.*eachindex(u))
        plus,minus = similar(u),similar(u)
        h = 1e-5
        rhs(plus,u+h*direction,nothing,0)
        rhs(minus,u-h*direction,nothing,0)
        @test J*direction ≈ (plus-minus)/(2h) rtol=2e-6
        @test u == before
        problem = catalytic_problem(system,(0,1))
        @test !problem.isoutofdomain(u,nothing,0)
        invalid = copy(u); invalid[end] = -1e-4
        @test problem.isoutofdomain(invalid,nothing,0)
        @test_throws ArgumentError ReactorSurface(:reactor,m;area=-1)
        @test_throws ArgumentError CatalyticNetwork(network;surfaces=())
        @test_throws ArgumentError CatalyticNetwork(network;surfaces=(ReactorSurface(:missing,m;area=1),))
        hot = ReactorNetwork((reactor=WellStirredReactor(gas;temperature=900,mole_fractions=ones(ng)),))
        @test_throws ArgumentError CatalyticNetwork(hot;surfaces=(surface,))
        @test_throws ArgumentError catalytic_problem(system,(1,0))
        @test_throws DimensionMismatch catalytic_problem(system,(0,1);initial_state=[1,2])
    end
end

@testset "Native ideal surfaces" begin
    mktempdir() do directory
        path = joinpath(directory,"surface.npz")
        m = simple_surface_fixture(path)
        s = IdealSurface(m)
        Cg = s.pressure/(Arrhenius.R*s.temperature)
        rates = surface_rates(s)
        @test rates.forward[1] ≈ 2*m.site_density*0.7*Cg
        @test rates.reverse[1] ≈ 2*m.site_density*0.3*Cg
        @test rates.gas[1] ≈ -only(rates.net)
        @test sum(rates.coverages) == 0
        @test maximum(abs,rates.elemental_rates) == 0
        @test rates.heat_release == 0
        @test IdealSurface(m;coverages=Dict("PT"=>7,"HS"=>3)).coverages ≈ s.coverages
        @test_throws ArgumentError IdealSurface(m;coverages=[-0.1,1.1])
        @test_throws ArgumentError IdealSurface(m;coverages=Dict("bad"=>1))
        @test_throws ArgumentError IdealSurface(m;temperature=0)
        @test_throws DimensionMismatch IdealSurface(m;mole_fractions=[0.2,0.8])
        rhs = surface_rhs(s)
        J = zeros(2,2)
        surface_jacobian!(J,s.coverages,rhs)
        @test J ≈ 2Cg*[-1 1;1 -1] rtol=1e-8
        p = surface_problem(s,(0,1))
        @test p.u0 == s.coverages
        @test p.isoutofdomain([-1e-4,1.0001],nothing,0)
        @test !p.isoutofdomain([0.0,1.0],nothing,0)
        @test_throws ArgumentError surface_problem(s,(1,0))
        @test_throws ArgumentError surface_problem(s,(0,1);initial_coverages=[-1e-3,1.001])
        @test surface_problem(s,(0,1);initial_coverages=[-1e-15,1.0]).u0 == [0,1]
        # Equal forward/reverse rates give the exact stationary coverage.
        eq = surface_rates(s,[0.5,0.5])
        @test maximum(abs,eq.coverages) < 1e-15
        fake_integrator(p;kwargs...) = (u=[p.u0,[0.5,0.5]],)
        steady = steady_coverages(s;integrator=fake_integrator)
        @test steady.coverages == [0.5,0.5]
        @test steady.residual < 1e-15
        for sites in (1.0,2.0), motz in (false,true)
            m = simple_surface_fixture(path;sticking=true,reversible=false,
                                       site_size=sites,coverage=true,motz_wise=motz)
            s = IdealSurface(m)
            q = 0.3
            probability = 0.2*10^(0.2q)*q^0.3 *
                exp(-(2e6q-3e6q^2+4e6q^3-1e6q^4)/(Arrhenius.R*900))
            motz && (probability /= 1-0.5probability)
            flux = probability*sqrt(Arrhenius.R*900/(2pi*1.008))*Cg*0.7
            r = surface_rates(s)
            @test r.forward[1] ≈ flux rtol=2e-14
            @test r.coverages[2] ≈ sites*flux/m.site_density
            @test r.reverse == [0.0]
            @test sum(r.coverages) == 0
            @test maximum(abs,r.elemental_rates) == 0
        end
        # Loading rejects mechanisms violating elemental or site conservation.
        data = npzread(path)
        data["products"][2,1] = 2
        npzwrite(path,data)
        @test_throws ArgumentError SurfaceMechanism(path)
    end
end
