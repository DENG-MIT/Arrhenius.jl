# Platinum catalytic combustion: the hydrogen-to-methane sequence from
# https://cantera.org/dev/examples/python/onedim/catalytic_combustion.html
# Prepare parameters once:
# python mechanism/export_surface.py ptcombust.yaml Pt_surf --output pt.surface.npz
# Run with a Julia environment containing Arrhenius, SciMLBase and OrdinaryDiffEqBDF:
# julia catalytic_combustion.jl pt.surface.npz
using Arrhenius
include(joinpath(@__DIR__,"..","reactors","surface_solver.jl"))

surface=SurfaceMechanism(ARGS[1])
gas=CreateSolution(surface.gas_file)
f=CatalyticImpingingJet(gas,surface;reactants="H2:.05,O2:.21,N2:.78,AR:.01",
    mdot=.06,T_inlet=300.,T_surface=900.,P=one_atm,width=.1,
    coverages=Dict("PT(S)"=>.5,"O(S)"=>.5))
initialize_catalytic_coverages!(f;integrator=native_surface_bdf)
set_catalytic_reactions!(f;gas_multiplier=0.,surface_multiplier=0.,coverage_enabled=false)
solve!(f)
for multiplier in 10.0.^(-5:0)
    set_catalytic_reactions!(f;gas_multiplier=multiplier,surface_multiplier=multiplier,coverage_enabled=true)
    solve!(f)
end
set_catalytic_inlet!(f,"CH4:.095,O2:.21,N2:.78,AR:.01")
solve!(f;ratio=100.,slope=.15,curve=.2,prune=0.)
println("Grid points: ",length(f.grid),"; wall temperature: ",temperature(f)[end]," K")
for (name,theta) in zip(surface.species_names[1:surface.n_surface],f.coverages)
    println(name,": ",theta)
end
