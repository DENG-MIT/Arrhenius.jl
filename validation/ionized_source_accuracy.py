"""Complete ionized profile/physical judge; timings are a separate gate."""
from pathlib import Path
import json,sys,numpy as np,cantera as ct
from flame_benchmarks import element_conservation
FIELDS=('T','Y','X','E','velocity','rho','qdot')
LIMITS=dict(T=.01,Y=.05,X=.05,E=.05,velocity=.01,rho=.01,qdot=.05)
def load(path):
 with np.load(path) as a:d={k:a[k] for k in a.files}
 if 'grid' not in d:d['grid']=d.pop('z')
 return d
def coords(d,case):
 z=d['grid'].copy()
 if case=='free':z-=np.interp(700.,d['T'],z)
 return z
def interp(values,z,q):
 a=np.asarray(values)
 return np.array([np.interp(q,z,v) for v in a]) if a.ndim==2 else np.interp(q,z,a)
def extrapolate(a,b,case):
 za,zb=coords(a,case),coords(b,case);q=np.union1d(za,zb)
 q=q[(q>=max(za[0],zb[0]))&(q<=min(za[-1],zb[-1]))]
 return dict(grid=q,**{k:2*interp(b[k],zb,q)-interp(a[k],za,q) for k in FIELDS})
def discrepancies(a,b,case,gas):
 za,zb=coords(a,case),coords(b,case);q=np.union1d(za,zb)
 q=q[(q>=max(za[0],zb[0]))&(q<=min(za[-1],zb[-1]))]
 errors={}
 for key in FIELDS:
  x,y=interp(a[key],za,q),interp(b[key],zb,q)
  scale=np.max(np.abs(x),axis=-1,keepdims=True)
  if key in ('Y','X'):
   floors=np.full(gas.n_species,1e-7)
   for i,(name,charge) in enumerate(zip(gas.species_names,gas.charges)):
    if charge:floors[i]=20*(1e-20 if name=='E' else 1e-16)
   if key=='X':floors=floors*float(np.max(1/np.sum(a['Y']/gas.molecular_weights[:,None],axis=0)))/gas.molecular_weights
   scale=np.maximum(scale,floors[:,None])
  elif np.max(scale)==0:
   errors[key]=0. if np.max(np.abs(y))<=1e-9 else float('inf');continue
  errors[key]=float(np.max(np.abs(x-y)/scale))
 return errors

def evaluate(mechanism,native,reference,case,out):
 if case not in ('ion-free','ion-burner'):
  raise ValueError('case must be ion-free or ion-burner')
 source_case=case.removeprefix('ion-')
 spec={'mechanism':str(mechanism),'cases':{source_case:{stage:str(Path(native)/f'{case}-{stage}-0.npz') for stage in ('frozen','field')}}}
 g=ct.Solution(spec['mechanism']);results={}
 for case,paths in spec['cases'].items():
  g.TPX=(300. if case=='free' else 600.),ct.one_atm,'CH4:1,O2:2,N2:7.52';inlet=g.Y.copy()
  for stage,path in paths.items():
   d=load(path);d['inlet_Y']=inlet;d['species_names']=g.species_names
   if 'X' not in d:d['X']=d['Y']/g.molecular_weights[:,None];d['X']/=d['X'].sum(axis=0)
   levels=[load(Path(reference)/f'{case}-{stage}-{i}.npz') for i in range(3)]
   assert all(len(b['grid'])==2*len(a['grid'])-1 and np.allclose(b['grid'][::2],a['grid'],rtol=0,atol=1e-14) for a,b in zip(levels,levels[1:]))
   fine=extrapolate(*levels[-2:],case);prior=extrapolate(*levels[:2],case)
   errors=discrepancies(fine,d,case,g);uncertainty=discrepancies(fine,prior,case,g)
   elements=element_conservation(g,d)
   speed=abs(d['velocity'][0]-fine['velocity'][0])/abs(fine['velocity'][0])
   speed_unc=abs(prior['velocity'][0]-fine['velocity'][0])/abs(fine['velocity'][0])
   rho=d['P'][0]/(ct.gas_constant*d['T']*np.sum(d['Y']/g.molecular_weights[:,None],axis=0)) if 'P' in d else ct.one_atm/(ct.gas_constant*d['T']*np.sum(d['Y']/g.molecular_weights[:,None],axis=0))
   mdot=rho*d['velocity'];continuity=float(np.max(np.abs(mdot-mdot[0]))/abs(mdot[0]))
   rhs=ct.faraday*rho*np.sum(g.charges[:,None]*d['Y']/g.molecular_weights[:,None],axis=0)/ct.epsilon_0
   gauss=float(np.max(np.abs(np.diff(d['E'])-np.diff(d['grid'])*rhs[1:]))/max(1.,float(np.max(np.abs(d['E']))))) if stage=='field' else float(np.max(np.abs(d['E'])))
   gates={k:errors[k]<=v for k,v in LIMITS.items()}
   gates.update(elements=all(e['pass'] for e in elements.values()),normalization=float(np.max(np.abs(d['Y'].sum(axis=0)-1)))<=1e-7,continuity=continuity<=1e-7,gauss=gauss<=(1e-6 if stage=='field' else 1e-9),finite=all(np.isfinite(d[k]).all() for k in FIELDS),positive_density=bool(np.min(rho)>0),rho_eos=bool(np.allclose(rho,d['rho'],rtol=1e-12,atol=0)),reference_uncertainty=all(uncertainty[k]<=.2*v for k,v in LIMITS.items()))
   if case=='free':gates.update(speed=speed<=.01,reference_speed_uncertainty=speed_unc<=.002)
   key=case+'-'+stage;results[key]=dict(passed=all(gates.values()),gates=gates,profile_errors=errors,reference_extrapolation_change=uncertainty,elements=elements,speed_error=float(speed),speed_reference_uncertainty=float(speed_unc),gauss_scaled_defect=gauss,continuity_relative_error=continuity,native_points=len(d['grid']),reference_points=[len(x['grid']) for x in levels])
   print(key,'PASS' if results[key]['passed'] else 'FAIL',[k for k,v in gates.items() if not v],flush=True)
 result=dict(passed=all(v['passed'] for v in results.values()),cases=results,limits=LIMITS,element_absolute_limit=1e-6,charged_profile_floor='20*source state atol; 1e-20 electron /1e-16 ion',reference='first-order extrapolation of independently solved whole-grid bisections; previous extrapolation bounds uncertainty',reference_uncertainty_fraction=.2)
 Path(out).write_text(json.dumps(result,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+'\n')
 return result
if __name__=='__main__':
 if len(sys.argv)!=6:raise SystemExit('usage: ionized_source_accuracy.py MECHANISM NATIVE REFERENCES CASE OUTPUT')
 sys.exit(0 if evaluate(*sys.argv[1:6])['passed'] else 1)
