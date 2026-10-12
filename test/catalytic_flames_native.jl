using Arrhenius, Test, NPZ, LinearAlgebra
if !isdefined(Arrhenius,:CatalyticImpingingJet)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","CatalyticFlames.jl"))
end

@testset "native catalytic impinging wall" begin
    mktempdir() do directory
        gas=CreateSolution(joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml"))
        ng=gas.n_species; kH=findfirst(==("H"),gas.species_names); kH2=findfirst(==("H2"),gas.species_names)
        # One catalytic site: H2 + PT -> 2 H + PT. The site coverage is exactly
        # one, allowing independent checks of wall flux signs and mass balance.
        elements=vcat(["Pt"],gas.elements)
        elemental=zeros(length(elements),ng+1); elemental[1,1]=1.
        elemental[2:end,2:end]=gas.ele_matrix
        reactants=zeros(ng+1,1); products=zeros(ng+1,1)
        reactants[1,1]=products[1,1]=1.
        reactants[kH2+1,1]=1.; products[kH+1,1]=2.
        coefficients=zeros(ng+1,15); coefficients[:,1].=298.15
        b(s)=collect(codeunits(s))
        path=joinpath(directory,"catalyst.npz")
        npzwrite(path,Dict{String,Any}(
            "surface_format_utf8"=>b("arrhenius-surface-v1"),"phase_name_utf8"=>b("site"),
            "species_names_utf8"=>b(join(vcat(["PT"],gas.species_names),'\n')),
            "element_names_utf8"=>b(join(elements,'\n')),"n_surface"=>[1],"n_gas"=>[ng],
            "site_density"=>[2e-5],"site_sizes"=>[1.],"molecular_weights"=>vcat([195.084],gas.MW),
            "elemental_matrix"=>elemental,"reactants"=>reactants,"products"=>products,"orders"=>reactants,
            "arrhenius"=>reshape([10.,0.,0.],1,3),"reversible"=>[false],"coverage_a"=>zeros(1,1),
            "coverage_m"=>zeros(1,1),"coverage_energy"=>zeros(1,1,4),"sticking_species"=>[0],
            "sticking_order"=>[0.],"sticking_factor"=>[0.],"motz_wise"=>[false],
            "thermo_type"=>fill(2,ng+1),"thermo_coefficients"=>coefficients,
            "reference_pressure"=>fill(one_atm,ng+1),"initial_coverages"=>[1.],
            "initial_mole_fractions"=>ones(ng)./ng,"initial_temperature"=>[900.],"initial_pressure"=>[one_atm]))
        m=SurfaceMechanism(path)
        f=Arrhenius.CatalyticImpingingJet(gas,m;reactants="H2:1,AR:4",mdot=.06,T_inlet=900.,T_surface=900.)
        @test_throws ArgumentError Arrhenius.set_catalytic_reactions!(f;gas_multiplier=-1.)
        @test_throws ArgumentError Arrhenius.set_catalytic_inlet!(f,"H2:1";mdot=0.)
        @test_throws ArgumentError Arrhenius.set_catalytic_inlet!(f,"H2:1";temperature=0.)
        @test_throws ArgumentError Arrhenius.solve!(f;ratio=1.)
        Arrhenius.set_catalytic_reactions!(f;gas_multiplier=0.,surface_multiplier=0.,coverage_enabled=false)
        Arrhenius.solve!(f)
        inert=similar(f.state); Arrhenius.counterflow_residual!(inert,f)
        @test f.converged
        @test maximum(abs.(temperature(f).-900.))<1e-6
        @test abs(velocity(f)[end])<1e-12
        @test f.coverages==[1.]
        Arrhenius.set_catalytic_reactions!(f;surface_multiplier=1.,coverage_enabled=true)
        reactive=similar(f.state); w=Arrhenius.CounterflowWorkspace(f)
        Arrhenius.counterflow_residual!(reactive,f,f.state,w)
        expected=m.stoichiometry[:,1].*w.surface.net_rates[1]
        @test reactive[1,end]==inert[1,end]
        @test reactive[kH+1,end]-inert[kH+1,end] ≈ expected[kH+1]*gas.MW[kH] rtol=1e-12
        @test reactive[kH2+1,end]-inert[kH2+1,end] ≈ expected[kH2+1]*gas.MW[kH2] rtol=1e-12
        Arrhenius.solve!(f)
        d=Arrhenius.catalytic_wall_diagnostics(f)
        @test f.converged
        @test d.diffusive_mass_flux[kH]<0
        @test d.diffusive_mass_flux[kH2]>0
        @test abs(d.diffusive_mass_flux[kH]+d.diffusive_mass_flux[kH2])<1e-10
        @test norm(d.wall_species_residual,Inf)<1e-10
        @test norm(d.elemental_production,Inf)<1e-15
        @test abs(d.total_mass_production)<1e-15
        @test d.coverage_sum==1.
        @test d.coverage_rates==[0.]
        @test maximum(abs.(temperature(f).-900.))<1e-6
        @test maximum(abs.(sum(mass_fractions(f);dims=1).-1))<1e-12
        Arrhenius.set_catalytic_inlet!(f,"H2:.8,AR:4")
        @test !isnothing(f.inlet_start) && !f.converged
    end
end
