#!/usr/bin/env python3
"""Exact final-state and deterministic-counter comparison for native case NPZs."""
import sys
from pathlib import Path
import numpy as np

if len(sys.argv) not in (3, 4):
    raise SystemExit("usage: compare_neutrality_npz.py <base-dir> <candidate-dir> [nh3-parity-root]")
base, candidate = map(Path, sys.argv[1:3])
nh3_root = Path(sys.argv[3]) if len(sys.argv) == 4 else None
if nh3_root:
    import cantera as ct
    import json
    specs = {p.stem: json.loads(p.read_text()) for d in nh3_root.glob("specs_*") if d.is_dir() for p in d.glob("*.json")}
    gas_cache = {}
    invariant_max = {"sumY_abs": 0.0, "element_mass_fraction_abs": 0.0, "planar_mass_flux_rel": 0.0}
state_keys = ("grid", "T", "Y", "velocity", "inlet_Y", "SL")
files = sorted(p.name for p in base.glob("*.npz"))
if len(files) != 26 or set(files) != {p.name for p in candidate.glob("*.npz")}:
    raise SystemExit(f"expected matching 26-case output sets; base={len(files)} candidate={len(list(candidate.glob('*.npz')))}")
all_exact = True
for name in files:
    with np.load(base/name, allow_pickle=False) as a, np.load(candidate/name, allow_pickle=False) as b:
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
        print(f"{name[:-4]}: state_exact={state_exact} counters_exact={counts_exact}")
        if nh3_root:
            spec = specs[name[:-4]]
            mechanism = spec["mechanism"]
            gas = gas_cache.setdefault(mechanism, ct.Solution(str(nh3_root / f"{mechanism}.yaml")))
            names = bytes(a["species_names_utf8"]).decode("utf-8").splitlines()
            if names != gas.species_names or names != bytes(b["species_names_utf8"]).decode("utf-8").splitlines():
                raise SystemExit(f"species order mismatch for {name}")
            elem_matrix = np.zeros((gas.n_elements, gas.n_species))
            for j, species in enumerate(names):
                comp = gas.species(species).composition
                mw = gas.molecular_weights[j]
                for e, atoms in comp.items():
                    ei = gas.element_index(e)
                    elem_matrix[ei, j] = atoms * gas.atomic_weights[ei] / mw
            for label, z in (("base", a), ("candidate", b)):
                ys = z["Y"].astype(float)
                sum_error = float(np.max(np.abs(ys.sum(axis=0) - 1.0)))
                elem = elem_matrix @ ys
                elem_error = float(np.max(np.abs(elem - elem[:, :1])))
                invariant_max["sumY_abs"] = max(invariant_max["sumY_abs"], sum_error)
                invariant_max["element_mass_fraction_abs"] = max(invariant_max["element_mass_fraction_abs"], elem_error)
                print(f"  {label}_invariants: max_abs_sumY={sum_error:.17g} max_abs_element_mass_fraction_drift={elem_error:.17g}")
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
if nh3_root:
    print("invariant_maxima=" + repr(invariant_max))
print(f"exact={sum(1 for p in files)} cases; verdict={'PASS' if all_exact else 'FAIL'}")
raise SystemExit(0 if all_exact else 1)