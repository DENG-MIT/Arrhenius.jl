import json
import sys

import cantera as ct


spec = json.load(open(sys.argv[1]))["spec"]
kind = sys.argv[2]
gas = ct.Solution(spec["mechanism"] + ".yaml")
pressure = spec["P_atm"] * ct.one_atm
gas.TPX = spec["T"], pressure, spec["X"]

if kind == "premixed":
    flame = ct.CounterflowPremixedFlame(gas, width=spec["width"])
    flame.reactants.mdot = spec["mdot"]
    gas.equilibrate("HP")
    flame.products.X = gas.X
    flame.products.T = gas.T
    flame.products.mdot = spec["mdot_products"]
elif kind == "twin":
    flame = ct.CounterflowTwinPremixedFlame(gas, width=spec["width"])
    flame.reactants.mdot = spec["mdot"]
else:
    raise ValueError(f"unsupported kind: {kind}")

gas.TPX = spec["T"], pressure, spec["X"]
flame.transport_model = spec["transport"]
flame.solve(loglevel=2, auto=False, refine_grid=False)
print(f"done points={len(flame.grid)} Tmax={max(flame.T):.8f}")
