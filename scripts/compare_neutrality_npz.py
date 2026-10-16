#!/usr/bin/env python3
"""Exact final-state and deterministic-counter comparison for native case NPZs."""
import sys
import json
from pathlib import Path
import numpy as np

if len(sys.argv) not in (3, 4):
    raise SystemExit("usage: compare_neutrality_npz.py <base-dir> <candidate-dir> [nh3-parity-root]")
base, candidate = map(Path, sys.argv[1:3])
nh3_root = Path(sys.argv[3]) if len(sys.argv) == 4 else None
if nh3_root:
    import cantera as ct
    import csv, json
    specs = {}
    for manifest in base.rglob("manifest.csv"):
        for row in csv.DictReader(manifest.open()):
            specs[row["id"]] = json.loads((nh3_root / row["spec"]).read_text())
    gas_cache = {}
    invariant_max = {"sumY_abs": 0.0, "planar_mass_flux_rel": 0.0}
state_keys = ("grid", "T", "Y", "velocity", "inlet_Y", "SL")
base_files = {p.name: p for p in base.rglob("*.npz")}
candidate_files = {p.name: p for p in candidate.rglob("*.npz")}
files = sorted(base_files)
if len(files) != 26 or set(files) != set(candidate_files):
    raise SystemExit(f"expected matching 26-case output sets; base={len(files)} candidate={len(candidate_files)}")
all_exact = True
counts_summary = {}
for name in files:
    with np.load(base_files[name], allow_pickle=False) as a, np.load(candidate_files[name], allow_pickle=False) as b:
        state_exact = True
        for key in state_keys:
            if key not in a or key not in b or not np.array_equal(a[key], b[key]):
                state_exact = False
                if key in a and key in b and a[key].shape == b[key].shape and a[key].size:
                    aa, bb = a[key].astype(float), b[key].astype(float)
                    scale = max(float(np.max(np.abs(aa))), float(np.max(np.abs(bb))), np.finfo(float).tiny)
                    print(f"  {key}: max_abs={np.max(np.abs(aa-bb)):.17g} max_rel={np.max(np.abs(aa-bb))/scale:.17g}")
                else:
                    print(f"  {key}: missing or shape differs")
        count_keys = sorted(k for k in a.files if k.startswith("stats_"))
        counts_exact = count_keys == sorted(k for k in b.files if k.startswith("stats_")) and all(
            np.array_equal(a[k], b[k]) for k in count_keys if k in b
        )
        if not counts_exact:
            print("  counters differ:", {k: (a[k].tolist() if k in a else None, b[k].tolist() if k in b else None)
                  for k in sorted(set(count_keys) | {k for k in b.files if k.startswith("stats_")})
                  if k not in a or k not in b or not np.array_equal(a[k], b[k])})
        if count_keys and counts_exact:
            print("  counters=" + repr({k: int(a[k][0]) for k in count_keys}))
        print(f"{name[:-4]}: state_exact={state_exact} counters_exact={counts_exact}")
        counts_summary[name[:-4]] = {
            "jacobians": int(a["stats_jacobians"][0]),
            "newton": int(a["stats_newton_iterations"][0]),
            "steps": int(a["stats_time_steps"][0]),
            "failures": int(a["stats_failures"][0]),
        }
        if nh3_root:
            spec = specs[name[:-4]]
            mechanism = spec["mechanism"]
            if mechanism not in gas_cache:
                gas_cache[mechanism] = ct.Solution(str(nh3_root / f"{mechanism}.yaml"))
            gas = gas_cache[mechanism]
            names = bytes(a["species_names_utf8"]).decode("utf-8").splitlines()
            if names != gas.species_names or names != bytes(b["species_names_utf8"]).decode("utf-8").splitlines():
                raise SystemExit(f"species order mismatch for {name}")
            for label, z in (("base", a), ("candidate", b)):
                ys = z["Y"].astype(float)
                sum_error = float(np.max(np.abs(ys.sum(axis=0) - 1.0)))
                invariant_max["sumY_abs"] = max(invariant_max["sumY_abs"], sum_error)
                print(f"  {label}_invariants: max_abs_sumY={sum_error:.17g}")
                if spec["kind"] in ("free", "burner", "fixed"):
                    pressure = float(spec["P_atm"]) * ct.one_atm
                    flux = []
                    for i, temp in enumerate(z["T"]):
                        gas.TPY = float(temp), pressure, ys[:, i]
                        flux.append(gas.density_mass * float(z["velocity"][i]))
                    flux = np.asarray(flux)
                    flux_rel = float(np.max(np.abs(flux - flux[0])) / max(abs(float(flux[0])), np.finfo(float).tiny))
                    invariant_max["planar_mass_flux_rel"] = max(invariant_max["planar_mass_flux_rel"], flux_rel)
                    print(f"  {label}_planar_mass_flux_rel={flux_rel:.17g}")
        all_exact &= state_exact and counts_exact
print("counts_json=" + json.dumps(counts_summary, sort_keys=True))
if nh3_root:
    print("invariant_maxima=" + repr(invariant_max))
print(f"exact={sum(1 for p in files)} cases; verdict={'PASS' if all_exact else 'FAIL'}")
raise SystemExit(0 if all_exact else 1)
