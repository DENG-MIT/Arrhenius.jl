"""Cantera 4 RK property references; no property tables enter native runtime.

Usage: python real_gas_cases.py CANTERA_DATA OUTPUT
The source shock-tube temperatures and 40 atm are retained. Additional 1/100 atm,
burned compositions and coefficient fixtures test the EOS independently of ODEs.
"""
from pathlib import Path
import sys
import json
import shutil
import hashlib
import cantera as ct
import numpy as np

source, output = map(Path, sys.argv[1:3])
output.mkdir(parents=True, exist_ok=True)
scalars = ["T", "P", "density", "enthalpy_mole", "int_energy_mole", "entropy_mole",
           "cp_mole", "cv_mole", "isothermal_compressibility", "thermal_expansion_coeff", "sound_speed"]
vectors = ["X", "partial_molar_int_energies_TV", "partial_molar_enthalpies",
           "partial_molar_volumes", "chemical_potentials"]


def save(name, gas, states, path):
    data = {key: [] for key in scalars + vectors}
    roots = []
    for entry in states:
        temperature, pressure, composition = entry[:3]
        root = entry[3] if len(entry) == 4 else 0
        if len(entry) == 4:
            gas.TDX = temperature, 1000 if root else 0.1, composition
        gas.TPX = temperature, pressure, composition
        roots.append(root)
        for key in data:
            data[key].append(getattr(gas, key))
    result = {key: np.array(value).T for key, value in data.items()}
    result["MW"] = gas.molecular_weights
    result["root"] = np.array(roots, dtype=np.int64)
    result["version_utf8"] = np.frombuffer(ct.__version__.encode(), dtype=np.uint8)
    np.savez(output / f"{name}.npz", **result)
    return {"name": name, "states": len(states), "species": gas.n_species,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


path = output / "dodecane.yaml"
shutil.copyfile(source / "nDodecane_Reitz.yaml", path)
gas = ct.Solution(str(path), "nDodecane_RK")
gas.set_equivalence_ratio(1, "c12h26", "o2:1,n2:3.76")
fresh = gas.X
temperatures = [1250, 1170, 1120, 1080, 1040, 1010, 990, 970, 950, 930, 910, 880, 850, 820, 790, 760]
states = [(t, 40*ct.one_atm, fresh) for t in temperatures]
states += [(t, p*ct.one_atm, comp) for t in (700, 1000, 2000, 3000)
           for p in (1, 40, 100) for comp in (fresh, "co2:12,h2o:13,n2:69.56")]
states += [(1200, 40*ct.one_atm, np.arange(1, gas.n_species+1))]
records = [save("dodecane", gas, states, path)]

base = ct.Solution("gri30.yaml")


def ordinary(obj):
    if isinstance(obj, dict):
        return {k: ordinary(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [ordinary(x) for x in obj]
    return obj


for name in ("binary-temperature", "critical-parameters", "ideal-limit"):
    species = [ordinary(base.species(n).input_data) for n in ("O2", "N2")]
    for i, spec in enumerate(species):
        if name == "critical-parameters":
            spec["critical-parameters"] = {"critical-temperature": [154.581, 126.192][i],
                                            "critical-pressure": [5.043e6, 3.3958e6][i]}
        else:
            spec["equation-of-state"] = {"model": "Redlich-Kwong",
                "a": [([1.74102e12, 1.55976e12][i] if name != "ideal-limit" else 0),
                      ([1.2e8, 2.3e8][i] if name != "ideal-limit" else 0)],
                "b": ([22.08100907, 26.81724983][i] if name != "ideal-limit" else 0)}
    if name == "binary-temperature":
        species[0]["equation-of-state"]["binary-a"] = {"N2": [1.3e12, 1.7e8]}
    doc = {"units": {"length": "cm", "quantity": "mol"},
           "phases": [{"name": "gas", "thermo": "Redlich-Kwong", "elements": ["O", "N"],
                       "species": ["O2", "N2"]}], "species": species}
    path = output / f"{name}.yaml"
    path.write_text(json.dumps(doc, indent=2))
    gas = ct.Solution(str(path))
    states = [(t, p, x) for t in (300, 700, 1200, 3000) for p in (1e4, ct.one_atm, 4e6, 1e7)
              for x in ("O2:1", "N2:1", "O2:.21,N2:.79")]
    records.append(save(name, gas, states, path))

# Explicit stable gas/liquid cubic roots, without asserting phase coexistence.
co2 = ordinary(ct.Solution(str(output / "dodecane.yaml"), "nDodecane_RK").species("co2").input_data)
doc = {"phases": [{"name": "gas", "thermo": "Redlich-Kwong", "species": ["co2"]}], "species": [co2]}
path = output / "co2-roots.yaml"
path.write_text(json.dumps(doc, indent=2))
gas = ct.Solution(str(path))
records.append(save("co2-roots", gas, [(280, 4e6, "co2:1", root) for root in (0, 1)], path))
(output / "provenance.json").write_text(json.dumps({"cantera": ct.__version__, "cases": records}, indent=2))
print(json.dumps(records))
