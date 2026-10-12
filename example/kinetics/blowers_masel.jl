using Arrhenius

# The reaction equations and parameters in Cantera's blowers_masel.py example.
function reaction_enthalpies(gas,T,index)
    x = mole_fractions(gas,"H2:1")
    h = cal_h_RT(gas,T,one_atm,x).*(R*T)
    return (h[index["H"]]+h[index["OH"]]-h[index["O"]]-h[index["H2"]],
            h[index["CH3"]]+h[index["H2"]]-h[index["H"]]-h[index["CH4"]])
end
function blowers_masel_calculation(gas)
    rate = BlowersMaselRate(3.87e1,2.7,6260*1000*4.184,1e9)
    index = Dict(name=>species_index(gas,name) for name in ("H","O","H2","OH","CH4","CH3"))
    temperatures = collect(300.:100.:3400.)
    rates = zeros(3,length(temperatures))
    for (j,T) in enumerate(temperatures)
        h2,methane = reaction_enthalpies(gas,T,index)
        rates[:,j] = [rate.A*T^rate.b*exp(-rate.Ea0/(R*T)),
                      rate_constant(rate,T,h2),rate_constant(rate,T,methane)]
    end
    # Sweep ±5 times the intrinsic barrier, matching the example's rate object
    # before its reaction enthalpy is assigned explicitly.
    enthalpies = collect(range(-5rate.Ea0,5rate.Ea0;length=100))
    barriers = activation_energy.(Ref(rate),enthalpies)
    return (;temperatures,rates,enthalpies,barriers)
end
if abspath(PROGRAM_FILE) == @__FILE__
    gas = CreateSolution(joinpath(@__DIR__,"..","..","mechanism","gri30.yaml"))
    result = blowers_masel_calculation(gas)
    println("Forward rate constants at 300 K: ",result.rates[:,1])
    println("Calculated ",length(result.temperatures)," temperature states and ",length(result.enthalpies)," enthalpy shifts.")
end
