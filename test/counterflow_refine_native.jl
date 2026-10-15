using Test, Arrhenius

@testset "axisymmetric midpoint refinement matches Cantera" begin
    gas=CreateSolution(joinpath(@__DIR__,"..","mechanism","h2o2.yaml"))
    ref=Arrhenius.npzread(joinpath(@__DIR__,"data","cantera_counterflow_refine.npz"))
    composition="H2:1.6,O2:1,AR:7"
    flows=(
        "premixed"=>CounterflowPremixedFlame(gas;reactants=composition,mdot_reactants=.12,
            mdot_products=.06,T_reactants=373.,P=.05*one_atm,width=.2),
        "twin"=>CounterflowTwinPremixedFlame(gas;reactants=composition,mdot=.12,
            T=373.,P=.05*one_atm,width=.2),
        "stagnation"=>ImpingingJet(gas;reactants=composition,mdot=.06,
            T_inlet=373.,T_surface=500.,P=.05*one_atm,width=.2))
    @test_throws ArgumentError solve!(flows[2].second;refinement_interpolation=:invalid)
    for (kind,f) in flows
        z=copy(f.grid); state=copy(f.state)
        @test Arrhenius._refine_axisymmetric!(f;ratio=3.,slope=.03,curve=.06,prune=.02,
            grid_min=1e-10,max_points=1200)
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
end
