"""Generate independent Cantera trajectories for native Julia reactor checks.

The H2/O2/Ar conditions adapt Cantera's reactor1.py (1001 K, one atmosphere,
2:1:4 dilution) to the repository's nine-species h2o2 mechanism. The isothermal
1400 K cases and the 1200 K, 2 MPa methane cases are selected regression points.
The exact mechanism files used by Julia are also loaded by this oracle.

References:
https://cantera.org/dev/examples/python/reactors/reactor1.html
https://cantera.org/dev/reference/reactors/ideal-gas-constant-pressure-reactor.html
https://cantera.org/dev/reference/reactors/ideal-gas-reactor.html
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import cantera as ct
import numpy as np


def generate_case(mechanism_dir: Path, output: Path, name: str, mechanism: str,
                  temperature: float, pressure: float, composition: str,
                  end_time: float, constraint: str, energy: str) -> dict:
    gas = ct.Solution(str(mechanism_dir / mechanism))
    gas.TPX = temperature, pressure, composition
    reactor_type = (ct.IdealGasConstPressureReactor if constraint == "constant_pressure"
                    else ct.IdealGasReactor)
    reactor = reactor_type(gas, energy="on" if energy == "adiabatic" else "off", clone=True)
    network = ct.ReactorNet([reactor])
    network.rtol, network.atol = 1e-11, 1e-17
    times = np.linspace(0.0, end_time, 401)
    state = np.empty((gas.n_species + 1, len(times)))
    rhs = np.empty_like(state)
    properties = {key: [] for key in ("pressure", "density", "enthalpy", "internal_energy")}
    for index, time in enumerate(times):
        if time > 0:
            network.advance(float(time))
        phase = reactor.phase if hasattr(reactor, "phase") else reactor.thermo
        state[:, index] = np.r_[phase.Y, phase.T]
        properties["pressure"].append(phase.P)
        properties["density"].append(phase.density)
        properties["enthalpy"].append(phase.enthalpy_mass)
        properties["internal_energy"].append(phase.int_energy_mass)
        production = phase.net_production_rates
        rhs[:-1, index] = production * phase.molecular_weights / phase.density
        species_energy = (phase.partial_molar_enthalpies if constraint == "constant_pressure"
                          else phase.partial_molar_int_energies)
        heat_capacity = phase.cp_mass if constraint == "constant_pressure" else phase.cv_mass
        rhs[-1, index] = (-np.dot(species_energy, production) / (phase.density * heat_capacity)
                          if energy == "adiabatic" else 0.0)
    target = output / f"reference_{name}.npz"
    np.savez(target, time=times, state=state, rhs=rhs,
             mechanism_utf8=np.frombuffer(mechanism.encode(), dtype=np.uint8),
             cantera_version_utf8=np.frombuffer(ct.__version__.encode(), dtype=np.uint8),
             constraint=np.array([constraint == "constant_volume"], dtype=np.int32),
             energy=np.array([energy == "isothermal"], dtype=np.int32),
             **{key: np.asarray(values) for key, values in properties.items()})
    return {"name": name, "mechanism": mechanism,
            "mechanism_sha256": hashlib.sha256((mechanism_dir / mechanism).read_bytes()).hexdigest(),
            "temperature_K": temperature,
            "pressure_Pa": pressure, "composition": composition, "end_time_s": end_time,
            "constraint": constraint, "energy": energy, "reference": target.name,
            "final_temperature_K": float(state[-1, -1])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanism-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "mechanism")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=("regression", "reactor1"), default="regression",
                        help="reactor1 requires the standard Cantera h2o2.yaml including N2")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cases = []
    if args.case == "reactor1":
        cases.append(generate_case(args.mechanism_dir, args.output,
            "reactor1_exact", "h2o2.yaml", 1001.0, ct.one_atm, "H2:2,O2:1,N2:4",
            0.001, "constant_pressure", "adiabatic"))
    else:
        for constraint, suffix in (("constant_pressure", "cp"), ("constant_volume", "cv")):
            for energy in ("adiabatic", "isothermal"):
                cases.append(generate_case(args.mechanism_dir, args.output,
                    f"h2_{suffix}_{energy}", "h2o2.yaml", 1001.0 if energy == "adiabatic" else 1400.0,
                    ct.one_atm, "H2:2,O2:1,AR:4", 0.001 if energy == "adiabatic" else 0.002,
                    constraint, energy))
            cases.append(generate_case(args.mechanism_dir, args.output,
                f"ch4_{suffix}_adiabatic", "gri30.yaml", 1200.0, 2e6,
                "CH4:1,O2:2,N2:7.52", 0.02, constraint, "adiabatic"))
    report = {"cantera_version": ct.__version__, "cases": cases}
    (args.output / "reference_manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
