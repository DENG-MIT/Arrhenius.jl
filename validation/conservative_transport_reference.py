"""Independent CT4 direct initialization and whole-grid refinement; no timing claims."""
import argparse,json
from pathlib import Path
import cantera as ct
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("mechanism",type=Path);p.add_argument("output",type=Path)
p.add_argument("--cases",nargs="+",default=["free"])
p.add_argument("--modes",nargs="+",default=["mole","mass","mole-soret","mass-soret","multi","multi-soret"])
p.add_argument("--levels",type=int,default=4)
p.add_argument("--burner-mdot",type=float,default=.06)
a=p.parse_args();a.output.mkdir(exist_ok=True,parents=True)
results=[]
for case in a.cases:
    for mode in a.modes:
        free=case=="free";g=ct.Solution(str(a.mechanism))
        g.TPX=(300.,ct.one_atm,"H2:1.1,O2:1,AR:5") if free else (373.,.05*ct.one_atm,"H2:1.5,O2:1,AR:7")
        inlet=g.Y.copy()
        E=np.array([[g.n_atoms(k,e)*g.atomic_weights[e]/g.molecular_weights[k] for k in range(g.n_species)] for e in range(g.n_elements)])
        f=ct.FreeFlame(g,width=.06) if free else ct.BurnerFlame(g,width=.5)
        if not free:f.burner.mdot=a.burner_mdot
        f.transport_model="multicomponent" if mode.startswith("multi") else "mixture-averaged"
        f.flux_gradient_basis="mass" if mode.startswith("mass") else "molar"
        f.soret_enabled=mode.endswith("soret")
        f.flame.jacobian_mode="finite-difference"
        f.flame.set_steady_tolerances(default=(1e-9,1e-14))
        f.set_refine_criteria(ratio=3.,slope=.06 if free else .05,curve=.12 if free else .1)
        f.max_grid_points=20000
        print("Direct",case,mode,flush=True)
        # auto=False keeps the physical domain fixed. Cantera's own native
        # initial profile and pseudo-transient Newton solve initialize the case.
        f.solve(loglevel=0,auto=False)
        if case=="fixed":
            f.flame.set_fixed_temp_profile([0.,.005,.01,.02,.05,.1,1.],[373.,650.,1000.,1350.,1650.,1750.,1750.])
            f.energy_enabled=False
            f.solve(loglevel=0,auto=False)
        for level in range(a.levels+1):
            if level:
                old=f.grid.copy();new=np.sort(np.r_[old,.5*(old[:-1]+old[1:])])
                fixed=f.fixed_temperature if free else None
                restart=ct.SolutionArray(g,len(new),extra={"grid":new,"velocity":np.interp(new,old,f.velocity)})
                restart.TPY=np.interp(new,old,f.T),f.P,np.stack([np.interp(new,old,y) for y in f.Y],axis=1)
                if free:
                    f.set_initial_guess(data=restart);f.fixed_temperature=fixed
                else:f.from_array(restart)
                f.solve(loglevel=0,refine_grid=False,auto=False)
            assert np.isclose(f.grid[-1]-f.grid[0],.06 if free else .5)
            d=dict(case=case,mode=mode,level=level,points=len(f.grid),speed=float(f.velocity[0]),Tmax=float(max(f.T)),
                element_mass_drift=float(np.max(np.abs(E@(f.Y[:,-1]-inlet)))))
            print(json.dumps(d),flush=True);results.append(d)
            np.savez(a.output/f"{case}-{mode}-{level}.npz",grid=f.grid,T=f.T,Y=f.Y,inlet_Y=inlet,
                velocity=f.velocity,density=f.density,element_matrix=E,P=[f.P])
            (a.output/"results.json").write_text(json.dumps(dict(version=ct.__version__,source_commit="726522be4e2a13454d8415b7ef799d621f665cf3",
                jacobian_mode="finite-difference",steady_tolerances=[1e-9,1e-14],results=results),indent=2))
