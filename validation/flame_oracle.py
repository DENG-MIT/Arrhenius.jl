"""Generate independent transport, equilibrium, and flame reference data."""
from __future__ import annotations
import argparse
import json
import platform
import time
from pathlib import Path

import cantera as ct
import numpy as np


def generate(mechanism: Path, output: Path, flame: bool):
    gas = ct.Solution(str(mechanism))
    output.mkdir(parents=True, exist_ok=True)
    states = []
    for composition in ("H2:1.1,O2:1,AR:5", "H2:2,O2:1,AR:3.76", "AR:1"):
        for temperature in (300., 800., 1500., 2400.):
            for pressure in (ct.one_atm, 5*ct.one_atm):
                gas.TPX = temperature, pressure, composition
                states.append([temperature, pressure, *gas.X, gas.viscosity,
                               gas.thermal_conductivity, *gas.mix_diff_coeffs])
    np.savetxt(output / "transport.csv", states, delimiter=",")
    np.savez(output / "transport.npz", states=np.array(states))
    equilibria = []
    for mode in ("TP", "HP"):
        for temperature in (300., 1000., 2000.):
            gas.TPX = temperature, ct.one_atm, "H2:1.1,O2:1,AR:5"
            x_initial = gas.X
            gas.equilibrate(mode)
            equilibria.append({"mode": mode, "Tin": temperature, "P": ct.one_atm,
                               "Xin": x_initial.tolist(), "T": gas.T, "X": gas.X.tolist()})
    (output / "equilibrium.json").write_text(json.dumps(equilibria, indent=2)+"\n")
    np.savez(output / "equilibrium.npz", Tin=[s["Tin"] for s in equilibria],
             T=[s["T"] for s in equilibria], X=np.array([s["X"] for s in equilibria]).T)
    metadata = {"cantera_version": ct.__version__, "cantera_commit": getattr(ct, "__git_commit__", None),
                "machine": platform.machine(), "platform": platform.platform(),
                "species": gas.species_names, "comparison_target": "Cantera 4.0", "minimum_speed_ratio": 0.95}
    if flame:
        gas.TPX = 300., ct.one_atm, "H2:1.1,O2:1,AR:5"
        f = ct.FreeFlame(gas, width=0.03)
        f.transport_model = "mixture-averaged"
        f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
        start = time.perf_counter()
        f.solve(loglevel=0, auto=True)
        metadata.update(flame_seconds=time.perf_counter()-start, speed=float(f.velocity[0]),
                        grid_points=len(f.grid), temperature_max=float(max(f.T)),
                        flux_gradient_basis="mole", soret=False)
        np.savez(output / "flame.npz", grid=f.grid, T=f.T, Y=f.Y, velocity=f.velocity,
                 density=f.density, mdot=f.density*f.velocity, inlet_Y=f.inlet.Y,
                 fixed_temperature=f.fixed_temperature)
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mechanism", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--flame", action="store_true")
    args = parser.parse_args()
    generate(args.mechanism, args.output, args.flame)
