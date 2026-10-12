"""Independent Cantera 4 counterflow validation; no data enters native solves."""
import argparse
from pathlib import Path
import cantera as ct
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("project",type=Path)
parser.add_argument("output",type=Path)
parser.add_argument("--kind",choices=["h2","methane","ethane"],default="h2")
parser.add_argument("--grid-from",type=Path,help="Use only the native grid; initialize all state fields independently")
parser.add_argument("--radiation",action="store_true",help="Repeat the solve with radiation on and grid refinement off")
parser.add_argument("--emissivities",default="0,0")
args=parser.parse_args()
mechanism=args.project/"mechanism"/("h2o2.yaml" if args.kind=="h2" else "gri30.yaml")
gas=ct.Solution(str(mechanism))
f=ct.CounterflowDiffusionFlame(gas,grid=np.load(args.grid_from)['grid']) if args.grid_from else ct.CounterflowDiffusionFlame(gas,width=.02)
f.fuel_inlet.X="H2:1,AR:1" if args.kind=="h2" else "CH4:1" if args.kind=="methane" else "C2H6:1"
f.oxidizer_inlet.X="O2:.2,AR:.8" if args.kind=="h2" else "O2:.21,N2:.78,AR:.01"
f.fuel_inlet.T=f.oxidizer_inlet.T=300
f.fuel_inlet.mdot=.24
f.oxidizer_inlet.mdot=.72
f.P=ct.one_atm
f.radiation_enabled=False
f.boundary_emissivities=tuple(map(float,args.emissivities.split(',')))
f.set_refine_criteria(ratio=4,slope=.2,curve=.3,prune=0)
f.solve(loglevel=0,auto=True,refine_grid=args.grid_from is None)
nonradiating_T=f.T.copy()
nonradiating_Lambda=f.L.copy()
if args.radiation:
    f.radiation_enabled=True
    f.solve(loglevel=0,refine_grid=False)
args.output.parent.mkdir(parents=True,exist_ok=True)
np.savez(args.output,grid=f.grid,T=f.T,Y=f.Y,velocity=f.velocity,spread_rate=f.spread_rate,
         Lambda=f.L,cantera_version=np.frombuffer(ct.__version__.encode(),dtype=np.uint8),
         radiative_heat_loss=f.flame.radiative_heat_loss,nonradiating_T=nonradiating_T,
         nonradiating_Lambda=nonradiating_Lambda,boundary_emissivities=f.boundary_emissivities)
print(dict(kind=args.kind,Tmax=max(f.T),grid_points=len(f.grid),Lambda=f.L[0]))
