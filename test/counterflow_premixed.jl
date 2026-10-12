using Test, Arrhenius, LinearAlgebra
if !isdefined(Arrhenius,:CounterflowPremixedFlame)
    Base.include(Arrhenius,joinpath(@__DIR__,"..","src","PremixedCounterflowFlames.jl"))
end

@testset "premixed axisymmetric boundaries" begin
    gas=CreateSolution(joinpath(dirname(pathof(Arrhenius)),"..","mechanism","h2o2.yaml"))
    composition="H2:1.6,O2:1,AR:7"
    kwargs=(;reactants=composition,mdot_reactants=.12,mdot_products=.06,
        T_reactants=373.,P=.05*one_atm,width=.2)
    @test_throws ArgumentError Arrhenius.CounterflowPremixedFlame(gas;kwargs...,products="H2O:1")
    @test_throws ArgumentError Arrhenius.CounterflowTwinPremixedFlame(gas;reactants=composition,mdot=-.1)
    @test_throws ArgumentError Arrhenius.ImpingingJet(gas;reactants=composition,mdot=.06,T_surface=0.)
    f=Arrhenius.CounterflowPremixedFlame(gas;kwargs...)
    solve!(f;slope=.1,curve=.2)
    @test f.converged
    @test 2000<maximum(temperature(f))<2100
    @test velocity(f)[1]>0 && velocity(f)[end]<0
    @test f.state[end,[1,end]] ≈ [.12,-.06] atol=1e-9
    @test maximum(abs.(sum(mass_fractions(f);dims=1).-1))<1e-10
    r=similar(f.state); Arrhenius.counterflow_residual!(r,f)
    @test norm(r,Inf)<1e-8

    # Inert wall and symmetry are checked directly with signed trial traces.
    twin=Arrhenius.CounterflowTwinPremixedFlame(gas;reactants=composition,mdot=.12,T=373.,P=.05*one_atm,width=.2)
    wall=Arrhenius.ImpingingJet(gas;reactants=composition,mdot=.06,T_inlet=373.,T_surface=500.,P=.05*one_atm,width=.2)
    for flow in (twin,wall)
        rr=similar(flow.state)
        Arrhenius.counterflow_residual!(rr,flow)
        @test rr[end,end]==0
        @test maximum(abs.(rr[2:gas.n_species+1,end]))<1e-12
        @test_throws ArgumentError Arrhenius.set_mass_flux!(flow;products=.01)
        Arrhenius.set_mass_flux!(flow;reactants=.07)
        @test flow.fuel_mass_flux==.07 && !flow.converged
        @test Arrhenius._premixed_reinitialize(flow,flow.grid).fuel_mass_flux==.07
    end
    rt=similar(twin.state); Arrhenius.counterflow_residual!(rt,twin)
    @test rt[1,end]==0
    @test rt[end-2,end] ≈ twin.state[end-2,end]-twin.state[end-2,end-1]
    rw=similar(wall.state); Arrhenius.counterflow_residual!(rw,wall)
    @test rw[1,end] ≈ 0 atol=1e-14
    @test rw[end-2,end]==0
    @test_throws ArgumentError solve!(wall;prune=.5)
    @test_throws ArgumentError solve!(wall;grid_min=0.)

    cold=Arrhenius.ImpingingJet(gas;reactants="AR:1",mdot=.06,T_inlet=373.,
        T_surface=500.,P=.05*one_atm,width=.2,initial_products=:inlet)
    solve!(cold)
    @test cold.converged
    @test Arrhenius.extinct(cold)
    @test temperature(cold)[[1,end]] ≈ [373.,500.] atol=1e-6
    @test minimum(temperature(cold)) >= 373. - 1e-6
    @test maximum(temperature(cold)) <= 500. + 1e-6
    @test abs(velocity(cold)[end])<1e-10
    @test maximum(abs.(mass_fractions(cold)[findfirst(==("AR"),gas.species_names),:].-1))<1e-10
end
