#!/usr/bin/env python3
import json
import sys
import numpy as np
import cantera as ct

def read_names(data, gas):
    if "species_names_utf8" not in data.files:
        return list(gas.species_names)
    return bytes(data["species_names_utf8"].tolist()).decode("utf-8").split("\n")

def compare(spec_path, arr_path, ct_path):
    spec = json.load(open(spec_path, encoding="utf-8"))
    arr = np.load(arr_path, allow_pickle=False)
    ref = np.load(ct_path, allow_pickle=False)
    gas = ct.Solution(spec["mechanism"] + ".yaml")
    arr_names = read_names(arr, gas)
    ct_names = read_names(ref, gas)
    if arr_names != ct_names:
        order = [arr_names.index(name) for name in ct_names]
        ya = np.asarray(arr["Y"])[order]
    else:
        ya = np.asarray(arr["Y"])
    za = np.asarray(arr["grid"], dtype=float)
    zr = np.asarray(ref["grid"], dtype=float)
    if not np.array_equal(za, zr):
        fields = {"T": (np.interp(zr, za, arr["T"]), ref["T"]),
                  "u": (np.interp(zr, za, arr["velocity"]), ref["velocity"]),
                  "Y": (np.vstack([np.interp(zr, za, row) for row in ya]), ref["Y"])}
    else:
        fields = {"T": (arr["T"], ref["T"]), "u": (arr["velocity"], ref["velocity"]),
                  "Y": (ya, ref["Y"])}
    errors = {}
    for key, (a, b) in fields.items():
        rel = np.abs(a - b) / np.maximum(np.abs(b), 1e-30)
        errors[key] = {"max_abs": float(np.max(np.abs(a-b))), "max_rel": float(np.max(rel)),
                       "max_rel_location": [int(v) for v in np.unravel_index(int(np.argmax(rel)), rel.shape)]}
        if key == "Y":
            per_species_pointwise = {name: float(np.max(rel[k])) for k, name in enumerate(ct_names)}
            per_species_scaled = {name: float(np.max(np.abs(a[k]-b[k])) / max(float(np.max(np.abs(b[k]))), 1e-12))
                                  for k, name in enumerate(ct_names)}
            errors[key]["per_species_pointwise_max_rel"] = per_species_pointwise
            errors[key]["per_species_scaled_rel"] = per_species_scaled
            errors[key]["max_abs_species"] = float(np.max(np.abs(a-b)))
            errors[key]["max_species_scaled_rel"] = max(per_species_scaled.values())
            errors[key]["worst_scaled_species"] = max(per_species_scaled, key=per_species_scaled.get)
    def mass_flux(data, species_names):
        order = [species_names.index(name) for name in gas.species_names]
        ys = np.asarray(data["Y"])[order]
        ts = np.asarray(data["T"])
        us = np.asarray(data["velocity"])
        vals = np.empty_like(ts, dtype=float)
        for j in range(len(ts)):
            gas.TPY = float(ts[j]), spec["P_atm"]*ct.one_atm, ys[:,j]
            vals[j] = gas.density*us[j]
        return np.interp(zr, np.asarray(data["grid"], dtype=float), vals)
    ma = (np.interp(zr, za, np.asarray(arr["mass_flux"], dtype=float))
          if "mass_flux" in arr.files else mass_flux(arr, arr_names))
    mr = mass_flux(ref, ct_names)
    mrel = np.abs(ma-mr)/np.maximum(np.abs(mr), 1e-30)
    result = {"case": spec["id"], "kind": spec["kind"], "native_grid_n": len(za), "ct_grid_n": len(zr),
              "grid_equal": bool(np.array_equal(za,zr)), "state_errors": errors,
              "mass_flux": {"arr_max_abs": float(np.max(np.abs(ma))), "ct_max_abs": float(np.max(np.abs(mr))),
                            "max_abs_delta": float(np.max(np.abs(ma-mr))), "max_rel_delta": float(np.max(mrel))}}
    print(json.dumps(result, sort_keys=True))

if len(sys.argv) < 4 or (len(sys.argv)-1) % 3:
    raise SystemExit("usage: fixed_t_state_compare.py spec.json arr.npz ct.npz [triples ...]")
for i in range(1, len(sys.argv), 3):
    compare(*sys.argv[i:i+3])
