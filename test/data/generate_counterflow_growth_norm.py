"""Generate Cantera 4.0 counterflow residual/growth regression values.

Reuses the lockstep ct_dump builder/writer. Run with a Cantera-enabled Python:
    python test/data/generate_counterflow_growth_norm.py
"""
import json
from pathlib import Path
import sys
import tempfile
import numpy as np
import cantera as ct

root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'test'/'cantera_parity'))
from ct_dump import build, write_state
spec=dict(kind='diffusion',mechanism=str(root/'mechanism'/'h2o2'),
          P_atm=15e6/ct.one_atm,T=800.,width=.03,fuel='H2:1',oxidizer='O2:1',
          mdot_fuel=.3,mdot_oxidizer=14.284,transport='mixture-averaged',
          ratio=2.,slope=.06,curve=.08)
f=build(spec); f.set_initial_guess(); f.energy_enabled=True
lam=np.array(f.flame.values('Lambda'))
with tempfile.TemporaryDirectory() as directory:
    # Lambda is algebraic. All transient terms vanish between these states.
    # Its gradient competes with the energy residual in physical units, while
    # the native scaled Lambda row dominates the energy row incorrectly.
    f.flame.set_profile('Lambda',np.linspace(0,1,len(f.grid)),lam+1e13*np.linspace(0,1,len(f.grid)))
    _,before,_=write_state(f,directory,'before',jacobian=False)
    f.flame.set_profile('Lambda',np.linspace(0,1,len(f.grid)),lam)
    _,after,_=write_state(f,directory,'after',jacobian=False)
names=['T']+f.gas.species_names+['spreadRate','Lambda','velocity']
out=dict(grid=f.grid,pressure=np.array([15e6]),mdot_fuel=np.array([.3]),
         mdot_oxidizer=np.array([14.284]),
         species_names_utf8=np.frombuffer('\n'.join(f.gas.species_names).encode(),dtype=np.uint8),
         cantera_version_utf8=np.frombuffer(ct.__version__.encode(),dtype=np.uint8))
for tag,data in [('before',before),('after',after)]:
    out[tag+'_values']=np.array([data['values/'+name] for name in names])
    out[tag+'_residual']=np.array([data['residual/'+name] for name in names])
np.savez(root/'test'/'data'/'cantera_counterflow_growth_norm.npz',**out)
print('reference growth:',(abs(out['before_residual']).max()/abs(out['after_residual']).max())**.2)
