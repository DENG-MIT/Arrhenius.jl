"""Cantera 4 reactor-network trajectories and prepared mechanism fixtures.

mix1 and combustor use the published examples' conditions and devices:
https://cantera.org/dev/examples/python/reactors/mix1.html
https://cantera.org/dev/examples/python/reactors/combustor.html
The closed two-vessel case is a selected conservation regression with a
time-dependent flow controller, a valve, and a fixed conductive wall.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import shutil

import cantera as ct
import numpy as np


def prepare_mechanisms(output: Path) -> Path:
    mechanism_directory = output / "mechanisms"
    mechanism_directory.mkdir(parents=True, exist_ok=True)
    exporter_path = Path(__file__).resolve().parents[1] / "mechanism" / "export_sidecar.py"
    spec = importlib.util.spec_from_file_location("sidecar_exporter", exporter_path)
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    for name in ("air.yaml", "gri30.yaml"):
        source = next(Path(directory) / name for directory in ct.get_data_directories()
                      if (Path(directory) / name).is_file())
        target = mechanism_directory / name
        shutil.copyfile(source, target)
        exporter.export(target, Path(str(target) + ".npz"))
    return mechanism_directory


def build_case(name: str, mechanism_directory: Path):
    def gas(mechanism="gri30.yaml", temperature=300.0, pressure=ct.one_atm,
            composition="O2:0.21,N2:0.78,AR:0.01"):
        phase = ct.Solution(str(mechanism_directory / mechanism))
        phase.TPX = temperature, pressure, composition
        return phase

    if name == "mix1":
        air = gas("air.yaml")
        fuel = gas(composition="CH4:1")
        source_a, source_b = ct.Reservoir(air, clone=True), ct.Reservoir(fuel, clone=True)
        outlet = ct.Reservoir(air, clone=True)
        mixer = ct.IdealGasReactor(gas(), clone=True, volume=1.0)
        inlet_a = ct.MassFlowController(source_a, mixer, mdot=air.density * 2.5 / 0.21)
        inlet_b = ct.MassFlowController(source_b, mixer, mdot=fuel.density)
        valve = ct.Valve(mixer, outlet, K=10.0)
        return [mixer], [inlet_a, inlet_b, valve], [], 2.0
    if name == "combustor":
        phase = gas(composition="CH4:1,O2:4,N2:15.04")
        inlet = ct.Reservoir(phase, clone=True)
        phase.equilibrate("HP")
        combustor = ct.IdealGasReactor(phase, clone=True, volume=1.0)
        exhaust = ct.Reservoir(phase, clone=True)
        controller = ct.MassFlowController(inlet, combustor, mdot=lambda t: combustor.mass / 0.1)
        outlet = ct.PressureController(combustor, exhaust, primary=controller, K=0.01)
        return [combustor], [controller, outlet], [], 2.0
    if name == "closed_pair":
        hot = ct.IdealGasReactor(gas(temperature=800.0, pressure=2*ct.one_atm,
                                    composition="N2:1"), clone=True, volume=0.02)
        cold = ct.IdealGasReactor(gas(temperature=400.0, composition="N2:1"), clone=True, volume=0.03)
        hot.chemistry_enabled = cold.chemistry_enabled = False
        controller = ct.MassFlowController(hot, cold, mdot=lambda t: 0.0005*(1+0.2*np.sin(3*t)))
        valve = ct.Valve(cold, hot, K=1e-8)
        wall = ct.Wall(hot, cold, A=0.1, U=10.0, Q=lambda t: 20*np.sin(t))
        return [hot, cold], [controller, valve], [wall], 2.0
    raise ValueError(name)


def capture_case(name: str, mechanisms: Path, output: Path) -> dict:
    reactors, devices, walls, end_time = build_case(name, mechanisms)
    network = ct.ReactorNet(reactors)
    network.rtol, network.atol = 1e-10, 1e-17
    network.initialize()
    times = np.unique(np.r_[0.0, np.geomspace(1e-8, 1e-2, 40), np.linspace(0.02, end_time, 100)])
    states, pressures, flows, heat_rates, total_energy = [], [], [], [], []

    def state_vector():
        return np.concatenate([np.r_[reactor.mass*reactor.phase.Y, reactor.T] for reactor in reactors])

    for time in times:
        if time > 0:
            network.advance(float(time))
        states.append(state_vector())
        pressures.append([reactor.phase.P for reactor in reactors])
        flows.append([device.mass_flow_rate for device in devices])
        heat_rates.append([wall.heat_rate for wall in walls])
        total_energy.append(sum(reactor.mass*reactor.phase.int_energy_mass for reactor in reactors))
    payload = dict(time=times, state=np.asarray(states).T, pressure=np.asarray(pressures).T,
                   mass_flow=np.asarray(flows).T, total_internal_energy=np.asarray(total_energy),
                   case_utf8=np.frombuffer(name.encode(), dtype=np.uint8),
                   cantera_version_utf8=np.frombuffer(ct.__version__.encode(), dtype=np.uint8))
    if walls:
        payload["wall_heat"] = np.asarray(heat_rates).T
    if name in ("mix1", "combustor"):
        network.solve_steady()
        payload["steady_state"] = state_vector()
    np.savez(output / f"reference_network_{name}.npz", **payload)
    return {"case": name, "cantera_version": ct.__version__, "sampled_points": len(times),
            "final_temperatures_K": [reactor.T for reactor in reactors],
            "reference": f"reference_network_{name}.npz"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    mechanisms = prepare_mechanisms(args.output)
    report = {"cases": [capture_case(name, mechanisms, args.output)
                        for name in ("mix1", "combustor", "closed_pair")]}
    (args.output / "network_reference_manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
