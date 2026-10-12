using Arrhenius

# Detached hydrogen flame against an inert wall, with mass-flow continuation.
# https://cantera.org/dev/examples/python/onedim/stagnation_flame.html
gas=CreateSolution(joinpath(@__DIR__,"..","..","mechanism","h2o2.yaml"))
f=ImpingingJet(gas;reactants="H2:1.8,O2:1,AR:7",mdot=.06,T_inlet=373.,
    T_surface=500.,P=.05*one_atm,width=.2)
for mdot in (.06,.07,.08,.09,.10,.11,.12)
    set_mass_flux!(f;reactants=mdot)
    solve!(f;ratio=3,slope=.1,curve=.2,prune=.06,grid_min=1e-4,loglevel=1)
    println("Mass flux: ",mdot," kg/(m² s); peak temperature: ",maximum(temperature(f))," K")
end
