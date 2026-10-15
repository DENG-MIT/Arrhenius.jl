using Test, Arrhenius

struct _Fix39NonIdealThermo <: Arrhenius.Thermo end
struct _Fix39NonIdealFlow
    gas::Any
end
struct _Fix39NonIdealCounterflow <: Arrhenius.AbstractPremixedCounterflow
    flow::_Fix39NonIdealFlow
end

function _fix39_assert_refined_mass_flux(f,old_grid,old_state,old_velocity)
    B=size(old_state,1)
    for i in eachindex(f.grid)
        j=findfirst(==(f.grid[i]),old_grid)
        isnothing(j) || continue
        left=searchsortedlast(old_grid,f.grid[i])
        @test left>0 && left<length(old_grid)
        midpoint=(old_state[:,left]+old_state[:,left+1])/2
        rho=Arrhenius._flame_density(f,1000*midpoint[1],@view(midpoint[2:f.gas.n_species+1])).density
        expected=rho*(old_velocity[left]+old_velocity[left+1])/2
        @test f.state[B,i] ≈ expected rtol=2e-14 atol=1e-15
    end
    properties=Arrhenius.CounterflowWorkspace(f).properties
    Arrhenius._flame_properties!(properties,f,f.state)
    @test properties.rho ≈ density(f) rtol=2e-14 atol=1e-15
end

@testset "axisymmetric midpoint refinement matches Cantera" begin
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    ref=Arrhenius.npzread(joinpath(@__DIR__,"data","cantera_counterflow_refine.npz"))
    composition="H2:1.6,O2:1,AR:7"
    make_flows()=(
        "premixed"=>CounterflowPremixedFlame(gas;reactants=composition,mdot_reactants=.12,
            mdot_products=.06,T_reactants=373.,P=.05*one_atm,width=.2),
        "twin"=>CounterflowTwinPremixedFlame(gas;reactants=composition,mdot=.12,
            T=373.,P=.05*one_atm,width=.2),
        "stagnation"=>ImpingingJet(gas;reactants=composition,mdot=.06,
            T_inlet=373.,T_surface=500.,P=.05*one_atm,width=.2))
    flows=make_flows()
    nonideal=Arrhenius.Solution(gas.n_species,gas.n_reactions,gas.MW,gas.species_names,
        gas.elements,gas.ele_matrix,_Fix39NonIdealThermo(),gas.trans,gas.reaction)
    err=try
        solve!(_Fix39NonIdealCounterflow(_Fix39NonIdealFlow(nonideal)))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("ideal-gas thermodynamics",sprint(showerror,err))
    @test_throws ArgumentError solve!(flows[2].second;refinement_interpolation=:invalid)
    for (kind,f) in flows
        # Same starting state as the Cantera fixture (recorded before Sim1D.refine).
        @test f.grid ≈ ref["$kind/pre/grid"] atol=1e-15
        @test temperature(f) ≈ ref["$kind/pre/T"] rtol=1e-9
        @test mass_fractions(f) ≈ ref["$kind/pre/Y"] atol=1e-9
        @test velocity(f) ≈ ref["$kind/pre/velocity"] rtol=1e-9
        @test spread_rate(f) ≈ ref["$kind/pre/spreadRate"] rtol=1e-9 atol=1e-12
        @test pressure_curvature(f) ≈ ref["$kind/pre/Lambda"] rtol=1e-9 atol=1e-12
        z=copy(f.grid); state=copy(f.state); uold=velocity(f)
        @test Arrhenius._refine_axisymmetric!(f;ratio=3.,slope=.03,curve=.06,prune=.02,
            grid_min=1e-10,max_points=1200)
        _fix39_assert_refined_mass_flux(f,z,state,uold)
        @test f.grid ≈ ref["$kind/grid"] atol=1e-15
        @test temperature(f) ≈ ref["$kind/T"] rtol=1e-9
        @test mass_fractions(f) ≈ ref["$kind/Y"] atol=1e-9
        @test velocity(f) ≈ ref["$kind/velocity"] rtol=1e-9
        @test spread_rate(f) ≈ ref["$kind/spreadRate"] rtol=1e-9
        @test pressure_curvature(f) ≈ ref["$kind/Lambda"] rtol=1e-9
        # Retained nodes preserve all stored components exactly.
        for j in eachindex(z)
            i=findfirst(==(z[j]),f.grid)
            isnothing(i) || @test f.state[:,i]==state[:,j]
        end
    end
    # Legacy rho*u interpolation misses Cantera's refined velocity (measured 0.21-0.64 m/s).
    for (kind,f) in make_flows()
        @test Arrhenius._refine_axisymmetric!(f;ratio=3.,slope=.03,curve=.06,prune=.02,
            grid_min=1e-10,max_points=1200,interpolate_velocity=false)
        @test length(f.grid)==length(ref["$kind/grid"])
        @test maximum(abs.(velocity(f)-ref["$kind/velocity"]))>1e-2
    end
end
