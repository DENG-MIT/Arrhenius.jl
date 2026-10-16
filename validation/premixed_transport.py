"""Cantera reference for premixed diffusion models, using the supplied mechanism."""
import argparse
import json
from pathlib import Path
import cantera as ct
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("mechanism")
p.add_argument("output", type=Path)
p.add_argument("--slope",type=float,default=.06)
p.add_argument("--curve",type=float,default=.12)
args = p.parse_args()
assert ct.__version__.startswith("4.0"), ct.__version__
args.output.mkdir(parents=True, exist_ok=True)
gas = ct.Solution(args.mechanism)
gas.TPX = 300., ct.one_atm, "H2:1.1,O2:1,AR:5"
f = ct.FreeFlame(gas, width=.03)
f.set_refine_criteria(ratio=3, slope=args.slope, curve=args.curve)
for model, basis, soret in (("mixture-averaged","mole",False),
                           ("mixture-averaged","mass",False),
                           ("mixture-averaged","mass",True),
                           ("multicomponent","mass",False),
                           ("multicomponent","mass",True)):
    f.transport_model = model
    f.flux_gradient_basis = "molar" if basis == "mole" else basis
    f.soret_enabled = soret
    f.solve(loglevel=0, auto=False)
    label = f"{model}-{basis}-{int(soret)}"
    arrays = dict(grid=f.grid,T=f.T,Y=f.Y,mdot=f.density*f.velocity,
                  fixed_temperature=np.array([f.fixed_temperature]),inlet_Y=f.inlet.Y,
                  speed=np.array([f.velocity[0]]),heat_release_rate=f.heat_release_rate)
    if soret:
        dt = []
        for j in range(len(f.grid)-1):
            gas.TPY = .5*(f.T[j]+f.T[j+1]),ct.one_atm,.5*(f.Y[:,j]+f.Y[:,j+1])
            dt.append(gas.thermal_diff_coeffs)
        arrays["thermal_diffusion"] = np.array(dt).T
    np.savez(args.output/f"{label}.npz", **arrays)
    print(json.dumps(dict(case=label,speed=f.velocity[0],Tmax=max(f.T),points=len(f.grid))))
