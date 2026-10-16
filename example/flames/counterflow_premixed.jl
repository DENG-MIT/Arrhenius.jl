using Arrhenius

# Lean hydrogen/oxygen opposed by equilibrium products.
# https://cantera.org/dev/examples/python/onedim/premixed_counterflow_flame.html
gas=CreateSolution(joinpath(@__DIR__,"..","..","mechanism","h2o2.yaml"))
f=CounterflowPremixedFlame(gas;reactants="H2:1.6,O2:1,AR:7",T_reactants=373.,
    P=.05*one_atm,mdot_reactants=.12,mdot_products=.06,width=.2)
solve!(f;ratio=3,slope=.1,curve=.2,prune=.02,loglevel=1)
println("Peak temperature: ",maximum(temperature(f))," K")
println("Pressure curvature: ",pressure_curvature(f)[1]," Pa/m²")
