"""Generate Cantera midpoint-refinement fixtures without solving a flame.

Run with Cantera 4.0 using the repository mechanism/h2o2.yaml. The initial
conditions match counterflow_initial_guess_native.jl; the output records
physical flow components before and after one Sim1D.refine call.
"""
from pathlib import Path
import sys
import cantera as ct
import numpy as np

root = Path(__file__).resolve().parents[2]
gas = ct.Solution(str(root / "mechanism/h2o2.yaml"))
records = {}
for kind in ("premixed", "twin", "stagnation"):
    gas.TPX = 373., .05 * ct.one_atm, "H2:1.6,O2:1,AR:7"
    if kind == "premixed":
        f = ct.CounterflowPremixedFlame(gas, width=.2)
        f.reactants.mdot = .12
        gas.equilibrate("HP")
        f.products.X, f.products.T, f.products.mdot = gas.X, gas.T, .06
        gas.TPX = 373., .05 * ct.one_atm, "H2:1.6,O2:1,AR:7"
    elif kind == "twin":
        f = ct.CounterflowTwinPremixedFlame(gas, width=.2)
        f.reactants.mdot = .12
    else:
        f = ct.ImpingingJet(gas, width=.2)
        f.inlet.mdot = .06
        f.surface.T = 500.
    f.set_initial_guess("equil") if kind == "stagnation" else f.set_initial_guess()
    f.energy_enabled = True
    # pre-refinement state: the native test asserts it starts from this same state
    records.update({f"{kind}/pre/{key}": value for key, value in {
        "grid": f.grid, "T": f.T, "Y": f.Y, "velocity": f.velocity,
        "spreadRate": f.spread_rate, "Lambda": f.L,
    }.items()})
    f.set_refine_criteria(ratio=3., slope=.03, curve=.06, prune=.02)
    f.refine(loglevel=0)
    records.update({f"{kind}/{key}": value for key, value in {
        "grid": f.grid, "T": f.T, "Y": f.Y, "velocity": f.velocity,
        "spreadRate": f.spread_rate, "Lambda": f.L,
    }.items()})
    print(kind, "points", len(f.grid), "velocity", f.velocity.tolist())
np.savez(sys.argv[1] if len(sys.argv)>1 else Path(__file__).with_name("cantera_counterflow_refine.npz"), **records)
