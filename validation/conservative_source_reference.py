"""Execute official example inputs/transport sequence, then independently refine."""
import argparse,json
from pathlib import Path
import numpy as np
import cantera as ct
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("parameters",type=Path);p.add_argument("output",type=Path)
p.add_argument("case",choices=["free","burner","fixed"])
p.add_argument("--levels",type=int,default=4)
p.add_argument("--free-width",type=float,default=.03)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
free=a.case=="free";fixed=a.case=="fixed";g=ct.Solution(str(a.parameters/("gri30.yaml" if fixed else "h2o2.yaml")))
g.TPX=(300.,ct.one_atm,"H2:1.1,O2:1,AR:5") if free else (373.7,ct.one_atm,"CH4:.65,O2:1,N2:3.76") if fixed else (373.,.05*ct.one_atm,"H2:1.5,O2:1,AR:7")
inlet=g.Y.copy();E=np.array([[g.n_atoms(k,e)*g.atomic_weights[e]/g.molecular_weights[k] for k in range(g.n_species)] for e in range(g.n_elements)])
f=ct.FreeFlame(g,width=a.free_width) if free else ct.BurnerFlame(g,width=.01 if fixed else .5)
if not free:f.burner.mdot=.04 if fixed else .06
if fixed:
    profile=np.load(a.parameters/"fixed-profile.npz")
    f.flame.set_fixed_temp_profile(profile["positions"]/.01,profile["temperatures"])
    f.energy_enabled=False
f.max_grid_points=20000
# An explicitly recorded Jacobian choice avoids the previously observed CT4
# analytic-Jacobian stall on fine burner grids; physical residuals are unchanged.
f.flame.jacobian_mode="finite-difference"
f.flame.set_steady_tolerances(default=(1e-9,1e-14))
results=[]
for stage,mode in enumerate(["mass","mass-soret","multi","multi-soret"] if free else ["mole","multi"]):
    multi=mode.startswith("multi")
    f.transport_model="multicomponent" if multi else "mixture-averaged"
    f.flux_gradient_basis="mass" if free else "molar"
    f.soret_enabled=mode.endswith("soret")
    slope=.06 if free else (.1 if multi else .3) if fixed else .05
    curve=.12 if free else (.2 if multi else 1.) if fixed else .1
    f.set_refine_criteria(ratio=3.,slope=slope,curve=curve)
    print("Official sequence",a.case,mode,"species",g.n_species,flush=True)
    f.solve(loglevel=0,auto=stage==0 and not fixed)
    source_state=f.to_array();source_anchor=f.fixed_temperature if free else None
    for level in range(a.levels+1):
        if level:
            old=f.grid.copy();new=np.sort(np.r_[old,.5*(old[:-1]+old[1:])])
            anchor=f.fixed_temperature if free else None
            restart=ct.SolutionArray(g,len(new),extra={"grid":new,"velocity":np.interp(new,old,f.velocity)})
            restart.TPY=np.interp(new,old,f.T),f.P,np.stack([np.interp(new,old,y) for y in f.Y],axis=1)
            if free:f.set_initial_guess(data=restart);f.fixed_temperature=anchor
            else:f.from_array(restart)
            f.solve(loglevel=0,refine_grid=False,auto=False)
        d=dict(case=a.case,mode=mode,level=level,points=len(f.grid),domain=[float(f.grid[0]),float(f.grid[-1])],
            speed=float(f.velocity[0]),Tmax=float(max(f.T)),element_mass_drift=float(np.max(np.abs(E@(f.Y[:,-1]-inlet)))))
        print(json.dumps(d),flush=True);results.append(d)
        np.savez(a.output/f"{a.case}-{mode}-{level}.npz",grid=f.grid,T=f.T,Y=f.Y,inlet_Y=inlet,
            velocity=f.velocity,density=f.density,P=[f.P],element_matrix=E)
        (a.output/(a.case+"-results.json")).write_text(json.dumps(dict(cantera_version=ct.__version__,jacobian_mode="finite-difference",
            source_commit="726522be4e2a13454d8415b7ef799d621f665cf3",results=results),indent=2))
    f.from_array(source_state)
    if free:f.fixed_temperature=source_anchor
