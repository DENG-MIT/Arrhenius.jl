"""Cantera surface-rate, coverage-trajectory, and stationary-state references.

The Pt/H2 point is the catalytic_combustion.py surface initialization at 900 K
and one atmosphere. The diamond point follows diamond_cvd.py at 1200 K, 20 Torr.
Additional positive coverage vectors and unequal gas/surface temperatures are
selected regression points, not published operating ranges.
https://cantera.org/dev/examples/python/onedim/catalytic_combustion.html
https://cantera.org/dev/examples/python/kinetics/diamond_cvd.html
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import cantera as ct
import numpy as np


def main(output):
    output.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("export_surface", Path(__file__).resolve().parents[1]/"mechanism"/"export_surface.py")
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    cases = [
        ("pt_h2", "ptcombust.yaml", "Pt_surf", 900.0, ct.one_atm,
         "H2:0.05,O2:0.21,N2:0.78,AR:0.01", {"PT(S)":0.5,"O(S)":0.5}),
        ("diamond", "diamond.yaml", "diamond_100", 1200.0, 20*ct.one_atm/760,
         "H:0.002,H2:0.988,CH3:0.0002,CH4:0.01", None),
    ]
    manifest = {"cantera_version":ct.__version__, "cases":{}}
    for name, mechanism, phase, T, P, composition, theta0 in cases:
        surface, gas_file = exporter.export_surface(mechanism, phase, output/f"{name}.surface.npz")
        gas = next(p for p in surface.adjacent.values() if p.thermo_model == "ideal-gas")
        surface.TP = T, P
        gas.TPX = T, P, composition
        if theta0 is not None:
            surface.coverages = theta0
        theta0 = surface.coverages.copy()
        names = surface.species_names + gas.species_names + [sp for p in surface.adjacent.values() if p is not gas for sp in p.species_names]
        order = [surface.kinetics_species_index(sp) for sp in names]
        rng = np.random.default_rng(4127)
        theta = [theta0]
        for k in range(8):
            q = np.exp(rng.uniform(-12, 1, surface.n_species))
            theta.append(q/q.sum())
        temperatures = np.array([T, 450, 700, 900, 1000, 1100, 1400, T, T],float)
        gas_temperatures = temperatures.copy()
        gas_temperatures[-2:] = [300, 1300]
        ropf, ropr, kf, kr, wdot = [], [], [], [], []
        for q, Ts, Tg in zip(theta, temperatures, gas_temperatures):
            surface.TP = Ts, P
            gas.TPX = Tg, P, composition
            surface.coverages = q
            for p in surface.adjacent.values():
                if p is not gas:
                    p.TP = Ts, P
            ropf.append(surface.forward_rates_of_progress.copy())
            ropr.append(surface.reverse_rates_of_progress.copy())
            kf.append(surface.forward_rate_constants.copy())
            kr.append(surface.reverse_rate_constants.copy())
            wdot.append(surface.net_production_rates[order].copy())
        surface.TP = T, P
        gas.TPX = T, P, composition
        for p in surface.adjacent.values():
            if p is not gas:
                p.TP = T, P
        surface.coverages = theta0
        times = np.concatenate(([0.0], np.geomspace(1e-12, 1.0, 121)))
        trajectory = [surface.coverages.copy()]
        for left, right in zip(times[:-1], times[1:]):
            surface.advance_coverages(right-left, rtol=1e-11, atol=1e-20, max_steps=100000)
            trajectory.append(surface.coverages.copy())
        surface.advance_coverages_to_steady_state()
        steady = surface.coverages.copy()
        production_steady = surface.net_production_rates[order].copy()
        np.savez(output/f"{name}.reference.npz", temperatures=temperatures,
            gas_temperatures=gas_temperatures, pressure=np.array([P]),
            temperature=np.array([T]), mole_fractions=gas.X.copy(),
            coverages=np.array(theta).T, forward=np.array(ropf).T,
            reverse=np.array(ropr).T, forward_constants=np.array(kf).T,
            reverse_constants=np.array(kr).T, production=np.array(wdot).T,
            times=times, trajectory=np.array(trajectory).T, initial_coverages=theta0,
            steady_coverages=steady, steady_production=production_steady)
        manifest["cases"][name] = {"mechanism":mechanism,"phase":phase,"temperature":T,
            "pressure":P,"composition":composition,"rate_states":len(theta),
            "trajectory_points":len(times),"gas_file":gas_file.name}
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    main(parser.parse_args().output)
