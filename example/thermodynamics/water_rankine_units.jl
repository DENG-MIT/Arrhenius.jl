using Arrhenius

# SI-equivalent inputs from Cantera's rankine_units.py (80.33 °F, 116.03 psi).
temperature_K = (80.33-32)*5/9+273.15
pressure_Pa = 116.03*6894.757293168364
cycle = water_rankine(PureWater();inlet_temperature=temperature_K,
    boiler_pressure=pressure_Pa,pump_efficiency=.6,turbine_efficiency=.8)
println("Pump work = ",cycle.pump_work," J/kg")
println("Turbine work = ",cycle.turbine_work," J/kg")
println("Heat added = ",cycle.heat_added," J/kg")
println("Efficiency = ",cycle.efficiency)
