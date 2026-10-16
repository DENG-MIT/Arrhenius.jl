using Arrhenius, LinearAlgebra

# https://cantera.org/dev/examples/python/transport/dusty_gas.html
function dusty_gas_calculation(g,X)
    gas=g.gas
    T=500.; P=one_atm
    diffusion=copy(dusty_gas_diffusion!(g,P,T,X))
    conductivity=dusty_gas_thermal_conductivity(g,P,T,X)
    Y=X.*gas.MW./dot(X,gas.MW)
    rho1=P*dot(X,gas.MW)/(R*T)
    rho2=1.2rho1
    uniform=dusty_gas_molar_fluxes(g,T,T,rho1,rho1,Y,Y,.001)
    pressure_gradient=dusty_gas_molar_fluxes(g,T,T,rho1,rho2,Y,Y,.001)
    return (;diffusion,conductivity,uniform,pressure_gradient)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS)==2 || error("usage: julia --project example/transport/dusty_gas.jl mechanism.yaml multicomponent.npz")
    gas=CreateSolution(ARGS[1])
    data=MultiTransportData(ARGS[2],gas)
    g=DustyGasTransport(gas;porosity=.2,tortuosity=4.,mean_pore_radius=1.5e-7,
        mean_particle_diameter=1.5e-6,multicomponent_data=data)
    X=mole_fractions(gas,Dict("OH"=>1.,"H"=>2.,"O2"=>3.,"O"=>1e-8,"H2"=>1e-8,
        "H2O"=>1e-8,"H2O2"=>1e-8,"HO2"=>1e-8,"AR"=>1e-8))
    result=dusty_gas_calculation(g,X)
    display(result.diffusion)
    println(result.conductivity)
    println(result.uniform)
    println(result.pressure_gradient)
end
