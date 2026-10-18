#!/usr/bin/env python3
"""Compare Arrhenius and production Cantera BurnerFlame startup arrays."""
import json
import sys

import numpy as np


def names(data):
    return bytes(data["species_names_utf8"].tolist()).decode("utf-8").split("\n")


def compare(spec_path, arr_path, ct_path):
    spec = json.load(open(spec_path, encoding="utf-8"))
    arr = np.load(arr_path, allow_pickle=False)
    ref = np.load(ct_path, allow_pickle=False)
    arr_order = [names(arr).index(name) for name in names(ref)]
    fields = {
        "grid": (arr["grid"], ref["grid"]),
        "T": (arr["T"], ref["T"]),
        "Y": (arr["Y"][arr_order], ref["Y"]),
        "velocity": (arr["velocity"], ref["velocity"]),
    }
    errors = {}
    relative = {}
    for key, (a, b) in fields.items():
        if a.shape != b.shape:
            raise AssertionError(f"{spec['id']} {key} shape differs: {a.shape} vs {b.shape}")
        errors[key] = float(np.max(np.abs(a - b)))
        relative[key] = errors[key] / max(np.finfo(float).tiny, float(np.max(np.abs(b))))
    limits = {"grid": 5e-15, "T": 3e-11, "Y": 4e-12, "velocity": 5e-11}
    failed = {key: (relative[key] if key != "grid" else errors[key])
              for key in fields
              if (relative[key] if key != "grid" else errors[key]) > limits[key]}
    print("INITIAL_COMPARE case={} points={} max_abs={} rel_inf={} species_order_equal={}".format(
        spec["id"], len(ref["grid"]), errors, relative, names(arr) == names(ref)), flush=True)
    if failed:
        raise AssertionError(f"{spec['id']} startup exceeds round-off limits: {failed}")
    return errors


if len(sys.argv) < 4 or (len(sys.argv) - 1) % 3:
    raise SystemExit("usage: fixed_t_initial_compare.py spec.json arr.npz ct.npz [triple ...]")
for i in range(1, len(sys.argv), 3):
    compare(*sys.argv[i:i + 3])
