"""Independent exact catalytic_combustion.py sequence; no native state seeding."""
import argparse
from pathlib import Path
import cantera as ct
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("output",type=Path)
p.add_argument("--tight",action="store_true")
a=p.parse_args()
a.output.mkdir(exist_ok=True,parents=True)
surf=ct.Interface("ptcombust.yaml","Pt_surf")
surf.TP=900.,ct.one_atm
gas=surf.adjacent["gas"]
gas.TPX=300.,ct.one_atm,"H2:.05,O2:.21,N2:.78,AR:.01"
surf.coverages="PT(S):.5,O(S):.5"
surf.advance_coverages_to_steady_state()
np.savez(a.output/"initial-coverages.npz",coverages=surf.coverages)
f=ct.ImpingingJet(gas,width=.1,surface=surf)
f.inlet.mdot=.06
f.inlet.T=300.
f.surface.T=900.
f.set_initial_guess(products="inlet")
if a.tight:
    f.flame.set_steady_tolerances(default=(1e-9,1e-14))
    f.surface.set_steady_tolerances(default=(1e-9,1e-14))

def save(name):
    # set the phase state explicitly before retrieving rates
    gas.TPY=f.T[-1],f.P,f.Y[:,-1]
    theta=np.array([f.surface.value(name) for name in surf.species_names])
    surf.coverages=theta
    rates=surf.net_production_rates
    gas.TPY=f.T[-2],f.P,f.Y[:,-2]
    xleft=gas.X
    gas.TPY=f.T[-1],f.P,f.Y[:,-1]
    xright=gas.X
    gas.TPY=.5*(f.T[-2]+f.T[-1]),f.P,.5*(f.Y[:,-2]+f.Y[:,-1])
    flux=-gas.density*gas.molecular_weights/gas.mean_molecular_weight*gas.mix_diff_coeffs*(xright-xleft)/(f.grid[-1]-f.grid[-2])
    flux-=f.Y[:,-2]*sum(flux)
    gas.TPY=f.T[-1],f.P,f.Y[:,-1]
    np.savez(a.output/(name+".npz"),grid=f.grid,T=f.T,Y=f.Y,velocity=f.velocity,
             spread_rate=f.spread_rate,Lambda=f.L,coverages=theta,
             surface_production_rates=rates,
             diffusive_mass_flux=flux,
             coverage_rates=rates[:surf.n_species]*np.array([s.size for s in surf.species()])/surf.site_density)
    print(name,dict(points=len(f.grid),Tmax=max(f.T),coverages=theta.tolist()),flush=True)

f.surface.coverage_enabled=False
surf.set_multiplier(0.)
gas.set_multiplier(0.)
f.solve(loglevel=0,auto=True)
save("inert")
f.surface.coverage_enabled=True
for exponent in range(-5,1):
    surf.set_multiplier(10.**exponent)
    gas.set_multiplier(10.**exponent)
    f.solve(loglevel=0)
    save("h2-"+str(exponent))
gas.TPX=300.,ct.one_atm,"CH4:.095,O2:.21,N2:.78,AR:.01"
f.inlet.X=gas.X
f.set_refine_criteria(ratio=100.,slope=.15,curve=.2,prune=0.)
f.solve(loglevel=0)
if a.tight:
    f.flame.set_steady_tolerances(default=(1e-9,1e-14))
    f.surface.set_steady_tolerances(default=(1e-9,1e-14))
    f.solve(loglevel=0,refine_grid=False)
save("methane")
