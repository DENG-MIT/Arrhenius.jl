"""Independent CT4 example cases; only a comparison grid may be imported."""
import argparse
from pathlib import Path
import cantera as ct
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("project",type=Path)
p.add_argument("kind",choices=["premixed","twin","wall"])
p.add_argument("output",type=Path)
p.add_argument("--grid-from",type=Path)
p.add_argument("--tight",action="store_true")
a=p.parse_args()
gas=ct.Solution(str(a.project/"mechanism"/("gri30.yaml" if a.kind=="twin" else "h2o2.yaml")))
grid_kwargs={"grid":np.load(a.grid_from)["grid"]} if a.grid_from else {"width":.025 if a.kind=="twin" else .2}
if a.kind=="premixed":
    gas.TPX=373.,.05*ct.one_atm,"H2:1.6,O2:1,AR:7"
    f=ct.CounterflowPremixedFlame(gas,**grid_kwargs)
    f.reactants.mdot=.12
    f.products.mdot=.06
    f.set_initial_guess()
elif a.kind=="twin":
    gas.TPX=300.,ct.one_atm,"CH4:.75,O2:2,N2:7.52"
    mdot=2*gas.density
    f=ct.CounterflowTwinPremixedFlame(gas,**grid_kwargs)
    f.reactants.mdot=mdot
else:
    gas.TPX=373.,.05*ct.one_atm,"H2:1.8,O2:1,AR:7"
    f=ct.ImpingingJet(gas,**grid_kwargs)
    f.inlet.mdot=.06
    f.surface.T=500.
    f.set_initial_guess(products="equil")
f.set_refine_criteria(ratio=2 if a.kind=="twin" else 3,
    slope=.3 if a.kind=="twin" else .1,curve=.3 if a.kind=="twin" else .2,
    prune=.05 if a.kind=="twin" else .06 if a.kind=="wall" else .02)
if a.kind=="wall":
    f.set_grid_min(1e-4)
f.solve(loglevel=0,auto=True,refine_grid=a.grid_from is None)
if a.tight:
    f.flame.set_steady_tolerances(default=(1e-9,1e-14))
    f.solve(loglevel=0,refine_grid=False)
np.savez(a.output,grid=f.grid,T=f.T,Y=f.Y,velocity=f.velocity,
    spread_rate=f.spread_rate,Lambda=f.L)
print(dict(kind=a.kind,points=len(f.grid),Tmax=max(f.T),Lambda=f.L[0]))
if a.kind=="twin":
    derivative=np.r_[np.diff(f.velocity)/np.diff(f.grid),(f.velocity[-1]-f.velocity[-2])/(f.grid[-1]-f.grid[-2])]
    location=np.abs(derivative).argmax()
    minimum=f.velocity[:location].argmin()
    point=np.abs(derivative[:minimum]).argmax()
    speed=np.trapezoid(f.heat_release_rate/f.cp,f.grid)/(max(f.T)-min(f.T))/max(f.density)
    np.savez(str(a.output).replace(".npz","-diagnostics.npz"),consumption_speed=speed,
        characteristic_strain_rate=abs(derivative[point]),strain_rate_point=point+1,
        strain_rate_profile=derivative,density=f.density,heat_release_rate=f.heat_release_rate)
