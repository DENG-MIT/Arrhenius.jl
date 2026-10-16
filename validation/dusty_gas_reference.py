"""Cantera4 DustyGas coefficients and endpoint fluxes; no values enter native initialization."""
import argparse,json
from pathlib import Path
import cantera as ct
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("mechanism",type=Path)
p.add_argument("output",type=Path)
a=p.parse_args()
g=ct.DustyGas(str(a.mechanism))
source={"OH":1.,"H":2.,"O2":3.,"O":1e-8,"H2":1e-8,"H2O":1e-8,
        "H2O2":1e-8,"HO2":1e-8,"AR":1e-8}
mixture="CH4:1,O2:2,N2:7.52" if "CH4" in g.species_names else "H2:2,O2:1,N2:3.76"
compositions=[source,"AR:1",mixture,np.arange(1,g.n_species+1,dtype=float)]
media=[(.2,4.,1.5e-7,1.5e-6,-1.),(.6,2.,5e-7,2e-6,0.),(.35,3.,2e-8,0.,1e-12)]
cases=[]
for T in (300.,500.,1500.,2500.):
    for P in (.1*ct.one_atm,ct.one_atm,10*ct.one_atm):
        for composition in compositions:
            for medium in media:
                # Isolate each property state from the pinned CT4 thermal-cache
                # history issue demonstrated in dusty_gas_cache_probe.py.
                g=ct.DustyGas(str(a.mechanism))
                por,tort,radius,diameter,permeability=medium
                g.porosity=por; g.tortuosity=tort
                g.mean_pore_radius=radius; g.mean_particle_diameter=diameter
                g.permeability=permeability
                g.TPX=T,P,composition
                X=g.X.copy(); Y1=g.Y.copy(); rho1=g.density
                diffusion=g.multi_diff_coeffs.copy()
                conductivity=g.thermal_conductivity
                T2=1.17*T
                X2=.92*X+.08*np.roll(X,1)
                g.TPX=T2,1.2*P,X2
                Y2=g.Y.copy(); rho2=g.density
                delta=.001
                flux=g.molar_fluxes(T,T2,rho1,rho2,Y1,Y2,delta)
                cases.append(dict(T=T,P=P,X=X,medium=np.array(medium),diffusion=diffusion,
                    conductivity=conductivity,T2=T2,rho1=rho1,rho2=rho2,Y1=Y1,Y2=Y2,delta=delta,flux=flux))
arrays={key:np.stack([c[key] for c in cases],axis=-1) for key in cases[0]}
arrays["species_utf8"]=np.frombuffer("\n".join(g.species_names).encode(),dtype=np.uint8)
arrays["cantera_version_utf8"]=np.frombuffer(ct.__version__.encode(),dtype=np.uint8)
np.savez(a.output,**arrays)

# Execute the exact source example calculation sequence independently.
g=ct.DustyGas(str(a.mechanism))
g.TPX=500.,ct.one_atm,source
g.porosity=.2; g.tortuosity=4.; g.mean_pore_radius=1.5e-7
g.mean_particle_diameter=1.5e-6; g.permeability=-1
diffusion=g.multi_diff_coeffs.copy()
conductivity=g.thermal_conductivity
T1,rho1,Y1=g.TDY
g.TP=T1,1.2*ct.one_atm
T2,rho2,Y2=g.TDY
zero=g.molar_fluxes(T1,T1,rho1,rho1,Y1,Y1,.001)
pressure=g.molar_fluxes(T1,T2,rho1,rho2,Y1,Y2,.001)
np.savez(str(a.output).replace(".npz","-example.npz"),diffusion=diffusion,
    conductivity=conductivity,zero_flux=zero,pressure_flux=pressure,
    T1=T1,T2=T2,rho1=rho1,rho2=rho2,Y1=Y1,Y2=Y2)
print(json.dumps(dict(species=g.n_species,cases=len(cases),conductivity=conductivity,
    pressure_flux=pressure.tolist(),cantera_version=ct.__version__),indent=2))
