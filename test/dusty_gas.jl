using Test, Arrhenius, LinearAlgebra
if !isdefined(Arrhenius,:DustyGasTransport)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","DustyGasTransport.jl"))
end
@testset "native dusty gas invariants" begin
    mechanism=get(ENV,"DUSTY_GAS_MECHANISM",joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml"))
    gas=CreateSolution(mechanism)
    args=(;porosity=.2,tortuosity=4.,mean_pore_radius=1.5e-7,mean_particle_diameter=1.5e-6)
    @test_throws ArgumentError DustyGasTransport(gas;args...,porosity=0.)
    @test_throws ArgumentError DustyGasTransport(gas;args...,tortuosity=-1.)
    @test_throws ArgumentError DustyGasTransport(gas;args...,permeability=-1.)
    w=DustyGasTransport(gas;args...)
    X=mole_fractions(gas,"H2:2,O2:1,AR:7")
    Y=X.*gas.MW./dot(X,gas.MW)
    rho=one_atm*dot(X,gas.MW)/(R*500)
    D=copy(dusty_gas_diffusion!(w,one_atm,500.,X))
    @test minimum(D)>=0
    @test w.resistance*D ≈ Matrix{Float64}(I,gas.n_species,gas.n_species) rtol=1e-12
    @test dusty_gas_diffusion!(w,one_atm,500.,X)==D
    @test w.knudsen.*sqrt.(gas.MW) ≈ fill(w.knudsen[1]*sqrt(gas.MW[1]),gas.n_species) rtol=1e-14
    set_porous_medium!(w;porosity=.4)
    @test dusty_gas_diffusion!(w,one_atm,500.,X) ≈ 2D rtol=1e-13
    set_porous_medium!(w;porosity=.2,tortuosity=8.)
    @test dusty_gas_diffusion!(w,one_atm,500.,X) ≈ .5D rtol=1e-13
    set_porous_medium!(w;tortuosity=4.)
    @test dusty_gas_permeability(w) ≈ .2^3*(1.5e-6)^2/(72*4*.8^2)
    @test all(iszero,dusty_gas_molar_fluxes(w,500.,500.,rho,rho,Y,Y,.001))
    flux=dusty_gas_molar_fluxes(w,500.,600.,rho,1.2rho,Y,Y,.001)
    reversed=dusty_gas_molar_fluxes(w,600.,500.,1.2rho,rho,Y,Y,.001)
    @test flux ≈ -reversed rtol=1e-13
    @test dusty_gas_molar_fluxes(w,500.,600.,rho,1.2rho,Y,Y,.002) ≈ .5flux rtol=1e-13
    @test sum(flux)<0 # net gas flow toward lower pressure
    @test_throws ArgumentError dusty_gas_diffusion!(w,one_atm,500.,2X)
    @test_throws ArgumentError dusty_gas_molar_fluxes(w,500.,500.,rho,rho,Y,Y,0.)
    @test_throws DimensionMismatch dusty_gas_molar_fluxes!(zeros(1),w,500.,500.,rho,rho,Y,Y,.001)
    @test_throws ArgumentError dusty_gas_thermal_conductivity(w,one_atm,500.,X)

    # A pure gas reduces to scalar Knudsen diffusion plus Darcy flow.
    pure=mole_fractions(gas,"AR:1"); index=findfirst(==("AR"),gas.species_names)
    density=one_atm*gas.MW[index]/(R*500)
    set_porous_medium!(w;permeability=1e-13)
    pureflux=dusty_gas_molar_fluxes(w,500.,500.,density,1.2density,pure,pure,.001)
    gradient=.2*one_atm/(R*500*.001)
    mean_concentration=1.1*one_atm/(R*500)
    exact=-w.knudsen[index]*gradient-mean_concentration*1e-13*.2*one_atm/(.001*w.gas_viscosity)
    @test pureflux[index] ≈ exact rtol=1e-13
    @test maximum(abs.(pureflux[setdiff(1:gas.n_species,[index])]))<1e-20
end
