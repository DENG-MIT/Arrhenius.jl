"""Official inert-wall mass-flow continuation with independent CT initialization."""
import sys
from pathlib import Path
import cantera as ct
import numpy as np

root,output=map(Path,sys.argv[1:3])
gas=ct.Solution(str(root/"mechanism/h2o2.yaml"))
gas.TPX=373.,.05*ct.one_atm,"H2:1.8,O2:1,AR:7"
f=ct.ImpingingJet(gas,width=.2)
f.inlet.mdot=.06
f.surface.T=500.
f.set_grid_min(1e-4)
f.set_refine_criteria(ratio=3,slope=.1,curve=.2,prune=.06)
f.set_initial_guess(products="equil")
for mdot in (.06,.07,.08,.09,.10,.11,.12):
    f.inlet.mdot=mdot
    f.solve(loglevel=0,auto=mdot==.06)
    np.savez(output/f"counterflow-wall-{round(100*mdot)}-reference.npz",grid=f.grid,T=f.T,Y=f.Y,
        velocity=f.velocity,spread_rate=f.spread_rate,Lambda=f.L,mdot=mdot)
    print(dict(mdot=mdot,points=len(f.grid),Tmax=max(f.T),Lambda=f.L[0]),flush=True)

    # Same-grid source comparison initialized solely from the independent CT
    # profile. Native state variables never enter this calculation.
    native=np.load(output/f"counterflow-wall-{round(100*mdot)}-native.npz")
    gas.TPX=373.,.05*ct.one_atm,"H2:1.8,O2:1,AR:7"
    probe=ct.ImpingingJet(gas,grid=native["grid"])
    probe.inlet.mdot=mdot
    probe.surface.T=500.
    probe.set_initial_guess(products="equil")
    locations=f.grid/f.grid[-1]
    for name,values in [("T",f.T),("velocity",f.velocity),("spreadRate",f.spread_rate),("Lambda",f.L)]:
        probe.flame.set_profile(name,locations,values)
    for name,values in zip(gas.species_names,f.Y):
        probe.flame.set_profile(name,locations,values)
    probe.flame.set_steady_tolerances(default=(1e-9,1e-14))
    probe.solve(loglevel=0,refine_grid=False)
    np.savez(output/f"counterflow-wall-{round(100*mdot)}-reference-samegrid.npz",grid=probe.grid,T=probe.T,Y=probe.Y,
        velocity=probe.velocity,spread_rate=probe.spread_rate,Lambda=probe.L,mdot=mdot)
