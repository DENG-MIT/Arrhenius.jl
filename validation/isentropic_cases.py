"""Cantera 4 references for the four nozzle and sound-speed example calculations.

Run: python validation/isentropic_cases.py directory

nozzle_calculation(gas, ...) and acoustics_calculation(gas, ...) are importable
for paired timing harnesses; they run the same source calculations on prepared
stateful gas objects (the caller owns state setup/reset) and write no files.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import cantera as ct
import numpy as np

SOURCE_COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"


def nozzle_calculation(gas, points, units=False):
    """thermo/isentropic.py(/_units.py) workload on a prepared gas.

    The gas must already hold the stagnation state (1200 K, 10 atm,
    H2:1,N2:0.1); the per-pressure SP state assignments below are the source's
    in-loop state changes.
    """
    s0, h0, p0, t0 = gas.entropy_mass, gas.enthalpy_mass, gas.P, gas.T
    rows, states = [], []
    for pressure in np.linspace(.01*p0, .99*p0, points):
        gas.SP = s0, pressure
        speed = np.sqrt(2*(h0-gas.enthalpy_mass))
        area = 1/(gas.density*speed)
        rows.append([area, speed/gas.sound_speed, gas.T if units else gas.T/t0, pressure/p0])
        states.append([gas.T, gas.P, gas.density, gas.enthalpy_mass, gas.entropy_mass])
    rows = np.array(rows)
    throat = rows[:, 0].min()
    rows[:, 0] /= throat
    return dict(data=rows, states=np.array(states), throat_area=np.array([throat]))


def nozzle(directory, name, points, units=False):
    gas = ct.Solution(str(Path(directory)/name))
    gas.TPX = 1200., 10*ct.one_atm, "H2:1,N2:0.1"
    return nozzle_calculation(gas, points, units)


def acoustics_calculation(gas, temperatures, rtol, units=False, refine_entropy=False):
    """thermo/sound_speed.py(/_units.py) workload on a prepared gas.

    The gas must already hold the fresh CH4:1,O2:2,N2:7.52 composition; per
    the published source, each equilibrium composition is carried into the
    next TP state. Temperatures are K, or degF when units=True (converted with
    the explicit SI formula, not pint). refine_entropy adds the independent
    frozen-SP entropy refinement (residual <1e-10 J/kg/K) used for the
    correctness references; it is not part of the published source.
    """
    rows, states = [], []
    for temperature in temperatures:
        kelvin = (temperature-32)*5/9+273.15 if units else temperature
        gas.TP = kelvin, ct.one_atm
        gas.equilibrate("TP", rtol=rtol, max_iter=5000)
        s0, p0, rho0 = gas.entropy_mass, gas.P, gas.density
        p1 = p0*1.0001
        gas.SP = s0, p1
        if refine_entropy:
            # The public SP setter has its own fixed temperature tolerance,
            # independent of equilibrate(rtol). Tighten the same frozen
            # isentrope before differentiating its small density increment.
            for _ in range(12):
                error = gas.entropy_mass-s0
                if abs(error) < 1e-10:
                    break
                gas.TP = gas.T-error*gas.T/gas.cp_mass, p1
            else:
                raise RuntimeError("frozen entropy refinement did not converge")
        frozen = np.sqrt((p1-p0)/(gas.density-rho0))
        gas.equilibrate("SP", rtol=rtol, max_iter=5000)
        equilibrium = np.sqrt((p1-p0)/(gas.density-rho0))
        velocities = np.array([equilibrium, frozen, gas.sound_speed])
        rows.append([temperature, *(velocities/.3048 if units else velocities)])
        states.append([gas.T, gas.P, *gas.X])
    return dict(data=np.array(rows), final_states=np.array(states).T)


def acoustics(directory, name, temperatures, rtol, units=False, refine_entropy=False):
    gas = ct.Solution(str(Path(directory)/name))
    gas.X = "CH4:1,O2:2,N2:7.52"
    return acoustics_calculation(gas, temperatures, rtol, units, refine_entropy)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args(argv)
    args.directory.mkdir(parents=True, exist_ok=True)
    assert ct.__version__.startswith("4.0")
    spec = importlib.util.spec_from_file_location(
        "exporter", Path(__file__).parents[1]/"mechanism"/"export_sidecar.py")
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    mechanisms = {}
    for name in ("h2o2.yaml", "gri30.yaml", "gri30_highT.yaml"):
        source = next(Path(d)/name for d in ct.get_data_directories() if (Path(d)/name).is_file())
        target = args.directory/name
        shutil.copyfile(source, target)
        exporter.export(target, Path(str(target)+".npz"))
        mechanisms[name] = hashlib.sha256(target.read_bytes()).hexdigest()

    outputs = {}
    for case, values in (
            ("isentropic", nozzle(args.directory, "h2o2.yaml", 200)),
            ("isentropic_units", nozzle(args.directory, "gri30.yaml", 10, units=True)),
            ("sound_speed", acoustics(args.directory, "gri30_highT.yaml",
                                      np.arange(300., 5001., 200.), 1e-8)),
            ("sound_speed_units", acoustics(args.directory, "gri30.yaml",
                                            np.linspace(80., 4880., 25), 1e-6, units=True)),
            ("sound_speed_refined", acoustics(args.directory, "gri30_highT.yaml",
                                              np.arange(300., 5001., 200.), 1e-11, refine_entropy=True)),
            ("sound_speed_units_refined", acoustics(args.directory, "gri30.yaml",
                                                    np.linspace(80., 4880., 25), 1e-11,
                                                    units=True, refine_entropy=True)),
            ("sound_speed_refined2", acoustics(args.directory, "gri30_highT.yaml",
                                               np.arange(300., 5001., 200.), 1e-13, refine_entropy=True)),
            ("sound_speed_units_refined2", acoustics(args.directory, "gri30.yaml",
                                                     np.linspace(80., 4880., 25), 1e-13,
                                                     units=True, refine_entropy=True))):
        for key, value in values.items():
            outputs[f"{case}_{key}"] = value
    np.savez(args.directory/"reference.npz", **outputs)
    metadata = {"cantera_version": ct.__version__, "mechanism_sha256": mechanisms,
                "source_commit": SOURCE_COMMIT,
                "scope": "All source calculation states; unit examples use identical SI calculations with explicit conversions; unit-wrapper overhead is not measured",
                "sound_speed_reference": "Published tolerance retained separately; two refinements use frozen entropy error <1e-10 J/kg/K and equilibrium rtol=1e-11, 1e-13"}
    (args.directory/"reference.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print("Generated 200+10 nozzle states and 24+25 acoustic states, plus refined acoustic references")


if __name__ == "__main__":
    main()
