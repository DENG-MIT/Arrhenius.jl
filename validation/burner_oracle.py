"""Reference burner calculations using the published low-pressure H2 example."""
import argparse
import json
from pathlib import Path
import cantera as ct
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("mechanism")
parser.add_argument("output",type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
gas = ct.Solution(args.mechanism)
gas.TPX = 373.,.05*ct.one_atm,"H2:1.5,O2:1,AR:7"
f = ct.BurnerFlame(gas,width=.5)
f.burner.mdot = .06
f.set_refine_criteria(ratio=3,slope=.05,curve=.1)
f.solve(loglevel=0,auto=True)
np.savez(args.output/"burner.npz",grid=f.grid,T=f.T,Y=f.Y,velocity=f.velocity,
         heat_release_rate=f.heat_release_rate)
print(json.dumps({"version":ct.__version__,"Tmax":float(max(f.T)),"points":len(f.grid)}))
# Reuse a prescribed profile to test the species BVP separately from energy.
positions = np.array([0.,.005,.01,.02,.05,.1,1.])
temperatures = np.array([373.,650.,1000.,1350.,1650.,1750.,1750.])
f2 = ct.BurnerFlame(gas,width=.5)
f2.burner.T = 373.
f2.burner.X = "H2:1.5,O2:1,AR:7"
f2.burner.mdot = .06
f2.energy_enabled = False
f2.flame.set_fixed_temp_profile(positions,temperatures)
f2.set_refine_criteria(ratio=3,slope=.05,curve=.1)
f2.solve(loglevel=0,auto=False)
np.savez(args.output/"burner-fixed.npz",grid=f2.grid,T=f2.T,Y=f2.Y,velocity=f2.velocity,
         positions=positions,temperatures=temperatures)
print(json.dumps({"fixed_Tmax":float(max(f2.T)),"points":len(f2.grid)}))
