using Arrhenius, LinearAlgebra

# Methane/air, equivalence ratio .75, inlet velocity 2 m/s; width is half-domain.
# https://cantera.org/dev/examples/python/onedim/premixed_counterflow_twin_flame.html
gas=CreateSolution(joinpath(@__DIR__,"..","..","mechanism","gri30.yaml"))
reactants="CH4:.75,O2:2,N2:7.52"
X=mole_fractions(gas,reactants)
mdot=2*one_atm*dot(X,gas.MW)/(R*300.)
f=CounterflowTwinPremixedFlame(gas;reactants,mdot,width=.025)
solve!(f;ratio=2,slope=.3,curve=.3,prune=.05,loglevel=1)
println("Peak temperature: ",maximum(temperature(f))," K")
println("Symmetry-plane axial velocity: ",velocity(f)[end]," m/s")
diagnostics=twin_flame_diagnostics(f)
println("Consumption speed: ",100*diagnostics.consumption_speed," cm/s")
println("Characteristic strain rate: ",diagnostics.characteristic_strain_rate," 1/s")
