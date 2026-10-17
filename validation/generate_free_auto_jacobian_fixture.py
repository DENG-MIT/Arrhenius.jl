#!/usr/bin/env python3
"""Map a Cantera full Jacobian into the native (T/1000, Y, m) coordinates."""
import argparse
import json
from pathlib import Path

import cantera as ct
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cantera-dump", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--native-residual", type=Path, required=True)
    p.add_argument("--mechanism", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--anchor", type=int, default=4)
    a = p.parse_args()

    data = np.load(a.cantera_dump)
    meta = json.loads(a.metadata.read_text())
    names, species = meta["component_names"], meta["species_names"]
    nsp, npoints = len(species), len(data["grid"])
    B, NC = nsp + 2, nsp + 6
    gas = ct.Solution(str(a.mechanism))
    mw = np.array([gas.molecular_weights[gas.species_index(s)] for s in species])
    active_names = ["T", *species, "velocity"]
    extras = ["spreadRate", "Lambda", "eField", "Uo"]
    permutation = np.array([
        j * NC + names.index(name)
        for j in range(npoints)
        for name in [*active_names, *extras]
    ])
    Jct = data["jacobian"][np.ix_(permutation, permutation)]

    # Match residual row units. Nonzero rows use measured same-state ratios;
    # zero rows use equation scales recorded in the source-level audit.
    native_residual = np.asarray(np.load(a.native_residual)["residual"]).reshape(npoints, B)
    Rct_original = np.zeros(npoints * NC)
    for j in range(npoints):
        for k, name in enumerate(active_names):
            Rct_original[j * NC + names.index(name)] = data[f"residual/{name}"][j]
    Rct = Rct_original[permutation]
    Rnative_original = np.zeros(npoints * NC)
    for j in range(npoints):
        for k, name in enumerate(active_names):
            Rnative_original[j * NC + names.index(name)] = native_residual[j, k]
    Rnative = Rnative_original[permutation]

    scales = np.ones(npoints * NC)
    anchor = a.anchor - 1
    for j in range(npoints):
        for k in range(B):
            i = j * NC + k
            if abs(Rct[i]) > 1e-24:
                scales[i] = Rnative[i] / Rct[i]
            if abs(Rct[i]) <= 1e-20 or abs(Rnative[i]) <= 1e-20:
                if k == 0:
                    scales[i] = 1e-3 if j in (0, npoints - 1) else 1e-7
                elif 1 <= k < B - 1:
                    scales[i] = 1.0 if j in (0, npoints - 1) else 1e-4
                elif j == anchor:
                    scales[i] = 1e-3
                elif j == npoints - 1:
                    scales[i] = 1.0
                else:
                    dz = (data["grid"][j + 1] - data["grid"][j] if j < anchor
                          else data["grid"][j] - data["grid"][j - 1])
                    scales[i] = -dz
        scales[j * NC + B:j * NC + NC] = 1.0
        scales[j * NC] = 1e-3 if j in (0, npoints - 1) else 1e-7
        if j not in (0, npoints - 1):
            scales[j * NC + 1:j * NC + B - 1] = 1e-4
        if j == anchor:
            scales[j * NC + B - 1] = 1e-3
        elif j == npoints - 1:
            scales[j * NC + B - 1] = 1.0
        else:
            dz = (data["grid"][j + 1] - data["grid"][j] if j < anchor
                  else data["grid"][j] - data["grid"][j - 1])
            scales[j * NC + B - 1] = -dz

    G = np.zeros((npoints * NC, npoints * NC))
    pressure, gas_constant = 101325.0, 8314.46261815324
    for j in range(npoints):
        y = np.array([data[f"values/{s}"][j] for s in species])
        temperature = float(data["values/T"][j])
        velocity = float(data["values/velocity"][j])
        invmw = np.sum(y / mw)
        rho = pressure / (invmw * gas_constant * temperature)
        q = j * NC
        G[q + names.index("T"), q] = 1000.0
        for k, name in enumerate(species):
            G[q + names.index(name), q + k + 1] = 1.0
        urow = q + names.index("velocity")
        G[urow, q] = velocity * 1000.0 / temperature
        G[urow, q + 1:q + nsp + 1] = velocity * (1.0 / mw) / invmw
        G[urow, q + B - 1] = 1.0 / rho
        for k, name in enumerate(extras):
            G[q + names.index(name), q + B + k] = 1.0

    G = G[permutation, :]
    transformed = (scales[:, None] * Jct) @ G
    active = np.array([j * NC + k for j in range(npoints) for k in range(B)])
    a.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.output, jacobian=transformed[np.ix_(active, active)])
    print(f"wrote {a.output}: {npoints} points, {B} variables/point")


if __name__ == "__main__":
    main()
