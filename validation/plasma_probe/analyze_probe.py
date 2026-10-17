#!/usr/bin/env python3
"""Compare native, Cantera, and LoKI-B EEDFs on a shared LXCat input.

The script consumes private, user-authorized LXCat-derived runtime files. It
prints only aggregate metrics; do not commit the source or solver EEDF files.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np

E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
K_B = 1.380649e-23
N_A = 6.02214076e23
AR_MOLAR_MASS = 39.948e-3


def read_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    return (np.array([float(x["energy_eV"]) for x in rows]),
            np.array([float(x["eedf"]) for x in rows]))


def read_loki(path: Path) -> tuple[np.ndarray, np.ndarray]:
    energy, eedf = [], []
    for line in path.read_text().splitlines()[1:]:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            energy.append(float(parts[0])); eedf.append(float(parts[1]))
        except ValueError:
            continue
    if not energy:
        raise ValueError(f"no numeric EEDF rows in {path}")
    return np.asarray(energy), np.asarray(eedf)


def interp_cross_section(energy: np.ndarray, table: dict) -> np.ndarray:
    x = np.asarray(table["energy_eV"], dtype=float)
    y = np.asarray(table["cross_section_m2"], dtype=float)
    result = np.interp(energy, x, y, left=0.0, right=y[-1])
    # Excitation and ionization cross sections vanish below their thresholds.
    if table["kind"] in ("excitation", "ionization"):
        result[energy < table["threshold_eV"]] = 0.0
    return result


def common_metrics(energy: np.ndarray, f: np.ndarray, tables: list[dict],
                   mobility: float, neutral_density: float) -> dict:
    norm = float(np.trapezoid(np.sqrt(energy) * f, energy))
    if norm <= 0 or not np.isfinite(norm):
        raise ValueError(f"invalid EEDF normalization {norm}")
    mean = float(np.trapezoid(energy ** 1.5 * f, energy) / norm)
    sigma = {t["kind"]: interp_cross_section(energy, t) for t in tables}
    vscale = math.sqrt(2.0 * E_CHARGE / M_E)
    k_ion = vscale * float(np.trapezoid(energy * sigma["ionization"] * f, energy)) / norm
    sigma_el = sigma["effective"] - sigma["excitation"] - sigma["ionization"]
    min_elastic = float(np.min(sigma_el))
    elastic_loss = vscale * (2.0 * M_E / (AR_MOLAR_MASS / N_A)) * float(
        np.trapezoid(energy ** 2 * sigma_el * f, energy)) / norm
    k_exc = vscale * float(np.trapezoid(energy * sigma["excitation"] * f, energy)) / norm
    k_ion_all = vscale * float(np.trapezoid(energy * sigma["ionization"] * f, energy)) / norm
    p_field = mobility * (1.0e-19 ** 2) * neutral_density
    p_loss_el = elastic_loss
    p_loss_exc = k_exc * 11.5
    p_loss_ion = k_ion_all * 15.8
    # In an ionizing steady EEDF, electron population growth carries its own
    # mean-energy term. This is the density-growth term in LoKI's power balance.
    p_loss_growth = k_ion_all * mean
    p_loss = p_loss_el + p_loss_exc + p_loss_ion + p_loss_growth
    return {"normalization": norm, "mean_energy_eV": mean,
            "mobility_m2_V_s": mobility, "ionization_rate_m3_s": k_ion,
            "energy_balance_field_eV_m3_s": p_field,
            "energy_balance_elastic_loss_eV_m3_s": p_loss_el,
            "energy_balance_excitation_loss_eV_m3_s": p_loss_exc,
            "energy_balance_ionization_loss_eV_m3_s": p_loss_ion,
            "energy_balance_electron_growth_loss_eV_m3_s": p_loss_growth,
            "energy_balance_relative_residual": (p_field-p_loss)/p_field,
            "min_reconstructed_elastic_cross_section_m2": min_elastic}


def relative(a: float, b: float) -> float:
    return abs(a-b) / abs(b) if b != 0 else abs(a-b)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--native", type=Path, required=True)
    ap.add_argument("--cantera", type=Path, required=True)
    ap.add_argument("--loki", type=Path, required=True)
    ap.add_argument("--swarm", type=Path, required=True)
    ap.add_argument("--tables", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--native-metadata", type=Path, required=True)
    ap.add_argument("--cantera-metadata", type=Path, required=True)
    ap.add_argument("--pressure-pa", type=float, default=133.32)
    ap.add_argument("--temperature-k", type=float, default=300.0)
    args = ap.parse_args()
    n = args.pressure_pa / (K_B * args.temperature_k)
    tables = json.loads(args.tables.read_text())
    distributions = {"native": read_csv(args.native), "cantera": read_csv(args.cantera),
                     "loki": read_loki(args.loki)}
    swarm_text = args.swarm.read_text()
    m = re.search(r"Reduced mobility coefficient\s*=\s*([+\-0-9.eE]+)", swarm_text)
    if not m:
        raise ValueError("LoKI reduced mobility not found")
    loki_mobility = float(m.group(1)) / n
    mobilities = {"native": None, "cantera": None, "loki": loki_mobility}
    # Native and Cantera mobilities are emitted in their run log and supplied
    # as the second and third optional fields in the same directory's metadata.
    native_meta = json.loads(args.native_metadata.read_text())
    cantera_meta = json.loads(args.cantera_metadata.read_text())
    mobilities["native"] = float(native_meta["native_mobility_m2_V_s"])
    mobilities["cantera"] = float(cantera_meta["cantera_mobility_m2_V_s"])
    summaries = {}
    for name, (energy, f) in distributions.items():
        summaries[name] = common_metrics(energy, f, tables, mobilities[name], n)
    # Compare normalized distributions at the native grid's points. This is
    # interpolation only; no solver output is rescaled except by its computed norm.
    base_e, base_f = distributions["native"]
    eedf_errors = {}
    for name in ("cantera", "loki"):
        e, f = distributions[name]
        f_on_base = np.interp(base_e, e, f)
        f_on_base /= summaries[name]["normalization"]
        base_normed = base_f / summaries["native"]["normalization"]
        eedf_errors[name] = float(np.max(np.abs(f_on_base-base_normed)) /
                                   np.max(np.abs(base_normed)))
        ion_tail = base_e >= 15.8
        eedf_errors[name + "_ion_tail"] = float(np.max(np.abs(f_on_base[ion_tail]-base_normed[ion_tail])) /
                                                 np.max(np.abs(base_normed[ion_tail])))
    comparisons = {}
    for name in ("cantera", "loki"):
        comparisons[name] = {k: relative(summaries[name][k], summaries["native"][k])
                             for k in ("mean_energy_eV", "mobility_m2_V_s",
                                       "ionization_rate_m3_s")}
        comparisons[name]["eedf_normalized_linf_relative"] = eedf_errors[name]
        comparisons[name]["eedf_ion_tail_linf_relative"] = eedf_errors[name + "_ion_tail"]
    result = {"state": {"gas": "Ar", "pressure_Pa": args.pressure_pa,
                         "temperature_K": args.temperature_k, "reduced_field_Td": 100.0,
                         "neutral_density_m-3": n},
              "summaries": summaries, "relative_to_native": comparisons}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
