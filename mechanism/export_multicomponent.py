"""Export temperature-independent fits for native Julia multicomponent transport.

Run with Cantera 4: python export_multicomponent.py mechanism.yaml output.npz
No transport properties at requested simulation states are stored in this sidecar.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import cantera as ct
import numpy as np


def export_multicomponent(mechanism: str, output: str | Path, phase: str = "") -> None:
    gas = ct.Solution(mechanism, phase) if phase else ct.Solution(mechanism)
    if gas.thermo_model != "ideal-gas" or np.any(gas.charges):
        raise ValueError("multicomponent export requires a neutral ideal gas")
    gas.transport_model = "multicomponent"
    species = gas.species()
    n = gas.n_species
    transports = [s.transport for s in species]
    eps = np.array([t.well_depth for t in transports])
    diameter = np.array([t.diameter for t in transports])
    dipole = np.array([t.dipole for t in transports])
    alpha = np.array([t.polarizability for t in transports])
    crot = np.array([{"atom": 0., "linear": 1., "nonlinear": 1.5}[t.geometry]
                     for t in transports])
    epsilon_ij = np.sqrt(eps[:, None] * eps[None, :])
    # GasTransport::makePolarCorrections, including the induced dipole correction.
    for j in range(n):
        for i in range(n):
            if (dipole[i] > 0) != (dipole[j] > 0):
                polar, nonpolar = (i, j) if dipole[i] > 0 else (j, i)
                alpha_star = alpha[nonpolar] / diameter[nonpolar] ** 3
                mu_star = dipole[polar] / np.sqrt(
                    4 * np.pi * ct.epsilon_0 * diameter[polar] ** 3 * eps[polar])
                xi = 1 + .25 * alpha_star * mu_star**2 * np.sqrt(eps[polar]/eps[nonpolar])
                epsilon_ij[i, j] *= xi**2
    astar = np.zeros((9, n, n))
    bstar = np.zeros_like(astar)
    cstar = np.zeros_like(astar)
    binary = np.zeros((5, n, n))
    for j in range(n):
        for i in range(n):
            astar[:, i, j], bstar[:, i, j], cstar[:, i, j] = (
                gas.get_collision_integral_polynomials(i, j))
            binary[:, i, j] = gas.get_binary_diff_coeffs_polynomial(i, j)
    utf8 = lambda text: np.frombuffer(text.encode("utf-8"), dtype=np.uint8)
    arrays = dict(
        format_utf8=utf8("arrhenius-multicomponent-v1"),
        species_names_utf8=utf8("\n".join(gas.species_names)),
        source_utf8=utf8(mechanism), cantera_version_utf8=utf8(ct.__version__),
        molecular_weights=gas.molecular_weights, epsilon_over_k=eps/ct.boltzmann,
        log_epsilon_ij_over_k=np.log(epsilon_ij/ct.boltzmann),
        rotational_heat_capacity=crot,
        rotational_relaxation=np.array([t.rotational_relaxation for t in transports]),
        viscosity_poly=np.array([gas.get_viscosity_polynomial(i) for i in range(n)]).T,
        binary_poly=binary, astar_poly=astar, bstar_poly=bstar, cstar_poly=cstar,
        gas_constant=np.array([ct.gas_constant]),
        temperature_range=np.array([gas.min_temp, gas.max_temp]),
    )
    if Path(mechanism).is_file():
        arrays["source_sha256_utf8"] = utf8(hashlib.sha256(Path(mechanism).read_bytes()).hexdigest())
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **arrays)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mechanism")
    parser.add_argument("output")
    parser.add_argument("--phase", default="")
    args = parser.parse_args()
    export_multicomponent(args.mechanism, args.output, args.phase)
