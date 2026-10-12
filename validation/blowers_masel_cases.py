"""Export the exact Blowers–Masel example's calculations as a Cantera 4 oracle."""
from pathlib import Path
import argparse
import importlib.util
import numpy as np
import cantera as ct

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output",type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
assert ct.__version__.startswith("4.0")
parameters = (3.87e1,2.7,6260*1000*4.184)
reactions = [ct.Reaction(equation="O + H2 <=> H + OH",rate=ct.ArrheniusRate(*parameters)),
             ct.Reaction(equation="O + H2 <=> H + OH",rate=ct.BlowersMaselRate(*parameters,1e9)),
             ct.Reaction(equation="H + CH4 <=> CH3 + H2",rate=ct.BlowersMaselRate(*parameters,1e9))]
gas = ct.Solution(thermo="ideal-gas",kinetics="gas",species=ct.Solution("gri30.yaml").species(),reactions=reactions)
gas.write_yaml(args.output/"blowers-masel.yaml")
spec = importlib.util.spec_from_file_location("sidecar",Path(__file__).parents[1]/"mechanism"/"export_sidecar.py")
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)
exporter.export(args.output/"blowers-masel.yaml",args.output/"blowers-masel.yaml.npz")
temperatures = np.arange(300,3500,100,dtype=float)
forward,reverse,enthalpy = [],[],[]
for T in temperatures:
    gas.TPX = T,ct.one_atm,"H2:1,O2:1,CH4:1,AR:1"
    forward.append(gas.forward_rate_constants.copy())
    reverse.append(gas.reverse_rate_constants.copy())
    enthalpy.append(gas.delta_enthalpy.copy())
# In Cantera 4 the standalone rate object's delta_enthalpy remains zero until
# explicitly assigned, so the source uses the intrinsic barrier for these limits.
E0 = gas.reaction(1).rate.activation_energy
deltaH = np.linspace(-5*E0,5*E0,100)
barriers = []
for desired in deltaH:
    index = gas.species_index("H")
    species = gas.species(index)
    coefficients = species.thermo.coeffs.copy()
    change = (desired-gas.delta_enthalpy[1])/ct.gas_constant
    coefficients[6] += change
    coefficients[13] += change
    species.thermo = ct.NasaPoly2(species.thermo.min_temp,species.thermo.max_temp,
                                species.thermo.reference_pressure,coefficients)
    gas.modify_species(index,species)
    rate = gas.reaction(1).rate
    rate.delta_enthalpy = gas.delta_enthalpy[1]
    barriers.append(rate.activation_energy)
np.savez(args.output/"reference.npz",T=temperatures,forward=np.array(forward).T,
         reverse=np.array(reverse).T,enthalpy=np.array(enthalpy).T,deltaH=deltaH,barriers=barriers)
print(ct.__version__,len(temperatures),"temperatures",len(deltaH),"enthalpy shifts")
