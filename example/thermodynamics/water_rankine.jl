using Arrhenius

water = PureWater()
cycle = water_rankine(water;inlet_temperature=300.,boiler_pressure=8e5,
                       pump_efficiency=.6,turbine_efficiency=.8)
for (i,state) in enumerate(cycle.states)
    println("State $i: T = ",state.T," K, P = ",state.P," Pa, Q = ",state.Q,
            ", h = ",state.h," J/kg, s = ",state.s," J/kg/K")
end
println("Pump work = ",cycle.pump_work," J/kg")
println("Turbine work = ",cycle.turbine_work," J/kg")
println("Heat added = ",cycle.heat_added," J/kg")
println("Efficiency = ",cycle.efficiency)
