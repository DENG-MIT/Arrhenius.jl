"""Independent Cantera references for NASA9, Shomate, and constant-cp.

Usage: python validation/species_thermo_cases.py CANTERA_DATA OUTPUT_DIRECTORY

NASA9 species come from upstream airNASA9.yaml. The CO Shomate coefficients
are Cantera's documented example; constant-cp fixtures exercise SI and legacy
cal/mol units. No property tables are used by the Julia evaluator.
https://cantera.org/stable/reference/thermo/species-thermo.html
https://cantera.org/stable/yaml/species.html
"""
from pathlib import Path
import json
import sys
import cantera as ct
import numpy as np

source, output = map(Path, sys.argv[1:3])
output.mkdir(parents=True, exist_ok=True)


def ordinary(value):
    if isinstance(value, dict):
        return {key: ordinary(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [ordinary(item) for item in value]
    return value


def document(species, units=None):
    elements = sorted({element for item in species for element in item["composition"]})
    result = {"phases": [{"name": "gas", "thermo": "ideal-gas", "elements": elements,
                          "species": [item["name"] for item in species]}], "species": species}
    if units:
        result["units"] = units
    return result


def reference(name, data, temperatures):
    yaml = json.dumps(data, indent=2)
    (output / f"{name}.yaml").write_text(yaml)
    gas = ct.Solution(yaml=yaml)
    records = {key: [] for key in ("T", "P", "cp_R", "h_RT", "s_R")}
    for temperature in temperatures:
        for pressure in (1e4, ct.one_atm, 1e7):
            gas.TP = temperature, pressure
            for key, value in (("T", gas.T), ("P", gas.P), ("cp_R", gas.standard_cp_R),
                               ("h_RT", gas.standard_enthalpies_RT), ("s_R", gas.standard_entropies_R)):
                records[key].append(value)
    arrays = {key: np.stack(value, axis=1) if key.endswith(("_R", "_RT")) else np.array(value)
              for key, value in records.items()}
    arrays["version_utf8"] = np.frombuffer(ct.__version__.encode(), dtype=np.uint8)
    arrays["MW"] = gas.molecular_weights
    arrays["elements"] = np.array([[gas.n_atoms(k, e) for k in range(gas.n_species)]
                                    for e in range(gas.n_elements)])
    np.savez(output / f"{name}.npz", **arrays)
    print(f"{name}: {gas.n_species} species, {len(records['T'])} states, Cantera {ct.__version__}", flush=True)


air_species = [ordinary(item.input_data) for item in ct.Species.list_from_file(str(source / "airNASA9.yaml"))]
reference("air-nasa9", document(air_species),
          (200., 298.15, 999.999999, 1000., 1000.000001, 3000., 5999.999999,
           6000., 6000.000001, 8000., 15000., 20000.))

n2 = next(item for item in air_species if item["name"] == "N2")
methane = next(ordinary(item.input_data) for item in ct.Species.list_from_file(str(source / "gri30.yaml"))
               if item.name == "CH4")
co = {"name": "CO", "composition": {"C": 1, "O": 1}, "thermo": {
    "model": "Shomate", "temperature-ranges": [298., 1300., 6000.],
    "data": [[25.56759, 6.096130, 4.054656, -2.671301, .131021, -118.0089, 227.3665],
             [35.15070, 1.300095, -.205921, .013550, -3.282780, -127.8375, 231.7120]]}}
argon = {"name": "AR", "composition": {"Ar": 1}, "thermo": {
    "model": "constant-cp", "T0": "298.15 K", "h0": "0 kJ/mol", "s0": "154.845 J/mol/K",
    "cp0": "20.786156545 J/mol/K", "T-min": "200 K", "T-max": "20000 K"}}
for item in (n2, methane, co, argon):
    item["thermo"]["reference-pressure"] = "1 bar"
reference("mixed-thermo", document([n2, methane, co, argon]),
          (298.15, 500., 999.999999, 1000., 1000.000001, 1299.999999,
           1300., 1300.000001, 3000., 5999.999999, 6000.))

legacy = document([
    {"name": "H", "composition": {"H": 1}, "thermo": {
        "model": "constant-cp", "T0": 1000., "h0": 9220., "s0": -3.02, "cp0": 5.95}},
    {"name": "H2", "composition": {"H": 2}, "thermo": {"model": "constant-cp"}},
], {"energy": "cal", "quantity": "mol", "pressure": "bar"})
reference("constant-cp-units", legacy, (200., 298.15, 1000., 2500., 6000.))

single_shomate = ordinary(co)
single_shomate["thermo"] = dict(co["thermo"])
single_shomate["thermo"]["temperature-ranges"] = [298., 1300.]
single_shomate["thermo"]["data"] = [co["thermo"]["data"][0]]
reference("single-shomate", document([single_shomate]), (298.15, 500., 1000., 1300.))

# Neutral air also exercises native equilibrium with the NASA9 data. Charged
# species properties above do not imply support for charge-constrained solving.
neutral = json.loads(json.dumps(air_species[:5]))
for item in neutral:
    item["thermo"]["reference-pressure"] = "1 atm"
neutral_doc = document(neutral)
reference("neutral-nasa9", neutral_doc, (300., 3500., 8000., 15000.))
gas = ct.Solution(yaml=json.dumps(neutral_doc))
temperatures, compositions = [], []
for temperature in (300., 3500., 8000., 15000.):
    gas.TPX = temperature, ct.one_atm, "N2:.79,O2:.21"
    gas.equilibrate("TP", rtol=1e-11, max_steps=2000)
    temperatures.append(gas.T)
    compositions.append(gas.X.copy())
np.savez(output / "neutral-nasa9-equilibrium.npz", T=np.array(temperatures), X=np.stack(compositions, axis=1))
