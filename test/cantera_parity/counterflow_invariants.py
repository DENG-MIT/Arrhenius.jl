"""Audit native counterflow continuity, closure and elemental rate balances.

Usage: counterflow_invariants.py MECHANISM_YAML FINAL_PROPERTIES_NPZ
The NPZ contains the native state, residual, transport rates and fluxes. Cantera
supplies only element/species stoichiometry. Element transport residuals are in
kg m^-3 s^-1 at interior points.
"""
import json,sys
import numpy as np
import cantera as ct

g=ct.Solution(sys.argv[1]);d=np.load(sys.argv[2]);x=d['state'];z=d['grid'];mw=d['MW'];n=len(mw)
a=np.array([[g.n_atoms(k,e)*g.atomic_weight(e)/mw[k] for k in range(n)] for e in range(g.n_elements)])
src=mw[:,None]*d['source'];element_src=a@src
scale=np.maximum(np.abs(a)@np.abs(src),1.)
physical=[]
for j in range(1,len(z)-1):
    cell=.5*(z[j+1]-z[j-1]); jm,jp=(j-1,j) if x[-1,j]>0 else (j,j+1)
    physical.append(src[:,j]-(d['flux'][:,j]-d['flux'][:,j-1])/cell-x[-1,j]*(x[1:n+1,jp]-x[1:n+1,jm])/(z[jp]-z[jm]))
element_residual=a@np.array(physical).T
result={'sumY_max':float(np.max(np.abs(np.sum(x[1:n+1],axis=0)-1))),
    'mass_flux_continuity_max':float(np.max(np.abs(d['residual'][-1]))),
    'native_residual_max':float(np.max(np.abs(d['residual']))),
    'elements_max':float(np.max(np.abs(element_residual))),
    'element_source_relative_max':float(np.max(np.abs(element_src)/scale)),
    'element_transport_units':'kg m^-3 s^-1',
    'by_element':{name:float(np.max(np.abs(element_residual[i]))) for i,name in enumerate(g.element_names)}}
print(json.dumps(result,indent=2))
