#!/usr/bin/env python
"""Compare native vs Cantera counterflow residuals (lockstep layer 1).

Usage: compare_cf.py DIR [ct_state|ct_solved]

DIR must hold the ct dump (ct_state.npz or ct_solved.npz + *_meta.json from
ct_dump.py) and native_cf.npz + native_cf_meta.json from native_dump_cf.jl,
both on the identical grid/state.

Native rows: [T/1000; Y_k; V/100; Lambda/10000; mdot]; raw native residuals
carry row scalings (ts=1e-4, Vscale=100, Lscale=1e4, dz), so for every
(component, region) pair the script estimates the scale s = median(native/ct)
over points with significant |ct| and reports:
  - the estimated scale s,
  - max |native - s*ct| relative to the regional max |ct| (should be ~1e-12
    if the discretizations are identical up to row scaling),
  - the spread of native/ct around s,
  - the argmax location of the deviation.
Rows that are structurally different (native dependent-species row = sum(Y)-1
at every point vs Cantera's full PDE at interior points / excess-species
normalization at boundaries; twin/stagnation wall species rows) are flagged
STRUCTURAL and reported separately.
"""
import json
import os
import sys

import numpy as np

d = sys.argv[1]
tag = sys.argv[2] if len(sys.argv) > 2 else "ct_state"
ct = np.load(os.path.join(d, f"{tag}.npz"))
meta = json.load(open(os.path.join(d, f"{tag}_meta.json")))
nt = np.load(os.path.join(d, f"native_cf_{tag}.npz"))
nmeta = json.load(open(os.path.join(d, f"native_cf_{tag}_meta.json")))

names = list(meta["species_names"])
n = len(names)
dep = int(nt["dependent_species"][0])
grid = ct["grid"]
N = len(grid)
kind = nmeta["kind"]
r_native = nt["residual"]          # (n+4, N)
r_ct = {c: ct[f"residual/{c}"] for c in meta["component_names"]}
v_ct = {c: ct[f"values/{c}"] for c in meta["component_names"]}

row_of = {"T": 0, "velocity": n + 3, "spreadRate": n + 1, "Lambda": n + 2}
for k, s in enumerate(names):
    row_of[s] = k + 1

print(f"kind={kind} N={N} species={n} native dependent = {names[dep]} "
      f"(native row {dep + 1})")
print(f"mirror_error = {nt['mirror_error'][0]:.3g}")
bd = {k[9:]: v for k, v in ct.items()
      if k.startswith("residual/") and k[9:] not in meta["component_names"]}
if bd:
    print("boundary-domain residuals:",
          {k: float(np.max(np.abs(v))) if v.size else 0.0 for k, v in bd.items()})

regions = {"left": [0], "interior": list(range(1, N - 1)), "right": [N - 1]}
if N <= 2:
    regions = {"all": list(range(N))}

print(f"\n{'comp':>10} {'region':>8} {'scale':>12} {'rel dev':>10} "
      f"{'ratio spread':>12} {'max|nat|':>10} {'max|ct|':>10}  worst point")
structural = []
for cname in meta["component_names"]:
    if cname == "eField" or cname == "Uo":
        continue
    k = row_of[cname]
    cn, rn = r_ct[cname], r_native[k]
    for rname, idx in regions.items():
        if not idx:
            continue
        c, r = cn[idx], rn[idx]
        scale = 1.0
        big = np.abs(c) > 1e-3 * max(np.abs(c).max(), 1e-300)
        if big.any():
            ratio = r[big] / c[big]
            scale = float(np.median(ratio))
        dev = np.abs(r - scale * c)
        cmax = max(np.abs(c).max(), 1e-300)
        rel = float(dev.max() / cmax) if cmax else float(dev.max())
        j = int(np.argmax(dev))
        spread = float(np.max(np.abs(r[big] / c[big] - scale))) if big.any() else 0.0
        flag = ""
        is_dep = (k == dep)
        if is_dep:
            flag = "STRUCT(dep row: native=sumY-1)"
        if rname == "right" and k >= 1 and k <= n and kind in ("twin", "stagnation") \
                and k != dep:
            flag = "STRUCT(wall: native=Y_N-Y_N-1, CT=flux BC)"
        print(f"{cname:>10} {rname:>8} {scale:12.4g} {rel:10.2e} "
              f"{spread:12.2e} {np.abs(r).max():10.3g} {cmax:10.3g}  "
              f"z={grid[idx[j]]:.5g} {flag}")
