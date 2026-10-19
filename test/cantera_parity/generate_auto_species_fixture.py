#!/usr/bin/env python3
"""Generate Cantera's species source derivative fixture for the Stagni 9-point state.

Run from the repository root with the project Cantera Python environment:
  python test/cantera_parity/generate_auto_species_fixture.py
The state fixture itself is generated with test/cantera_parity/ct_dump.py and
 test/data/cantera/free_auto_stagni.json.
"""
from pathlib import Path
import cantera as ct
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "test" / "data" / "cantera"
mechanism = ROOT / "mechanism" / "stagni_2023.yaml"
state = np.load(DATA / "free_auto_stagni.npz")
gas = ct.Solution(str(mechanism))
N = len(state["grid"])
K = gas.n_species
jac = np.zeros((K, K, N))
for j in range(1, N - 1):
    Y = np.array([state[f"values/{name}"][j] for name in gas.species_names])
    gas.TPY = float(state["values/T"][j]), 101325.0, Y
    ddCi = gas.net_production_rates_ddCi
    concentrations = gas.concentrations
    rho = gas.density
    wbar = gas.mean_molecular_weight
    jac[:, :, j] = (rho * ddCi - wbar * (ddCi @ concentrations)[:, None]) / gas.molecular_weights[None, :]
np.savez_compressed(DATA / "free_auto_stagni_source.npz", source_jacobian=jac)
print(f"wrote {K}x{K} source derivative at {N-2} interior nodes")
