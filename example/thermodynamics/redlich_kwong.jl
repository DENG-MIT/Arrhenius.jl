using Arrhenius

# Initial thermodynamic state of Cantera's non_ideal_shock_tube.py example.
# This example evaluates the native EOS; it does not integrate ignition.
model = RedlichKwongThermo(joinpath(@__DIR__,"..","..","mechanism","nDodecane_Reitz.yaml"))
state = redlich_kwong_state(model;T=1000.,P=40one_atm,X="c12h26:1,o2:18.5,n2:69.56")
println("Redlich–Kwong density: ",state.rho," kg/m³; compressibility factor: ",state.Z)
println("Specific internal energy: ",state.u_mass," J/kg; cv: ",state.cv_mass," J/kg/K")
println("Ideal-gas density at the same T/P/X: ",state.P*state.MW/(R*state.T)," kg/m³")

# Reuse arrays when a reactor repeatedly evaluates T and composition at fixed density.
work = RedlichKwongWorkspace(model)
properties = redlich_kwong_properties!(work,model,state.T,state.rho,state.X)
println("Constant-volume species energy coefficient for fuel: ",properties.u_TV[1]," J/kmol")
