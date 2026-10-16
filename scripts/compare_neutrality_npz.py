#!/usr/bin/env python3
"""Exact final-state and deterministic-counter comparison for native case NPZs."""
import sys
from pathlib import Path
import numpy as np

if len(sys.argv) != 3:
    raise SystemExit("usage: compare_neutrality_npz.py <base-dir> <candidate-dir>")
base, candidate = map(Path, sys.argv[1:])
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
        all_exact &= state_exact and counts_exact
print(f"exact={sum(1 for p in files)} cases; verdict={'PASS' if all_exact else 'FAIL'}")
raise SystemExit(0 if all_exact else 1)