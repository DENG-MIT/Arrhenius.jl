#!/usr/bin/env python3
import json
import os
import sys
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "1"
import cantera as ct
import numpy as np


def run_case(spec_path, out):
    spec = json.load(open(spec_path, encoding="utf-8"))
    gas = ct.Solution(spec["mechanism"] + ".yaml")
    gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
    inlet_y = gas.Y.copy()
    f = ct.BurnerFlame(gas, width=spec["width"])
    f.burner.mdot = spec["mdot"]
    f.transport_model = spec["transport"]
    f.set_refine_criteria(ratio=spec["ratio"], slope=spec["slope"], curve=spec["curve"],
                          prune=min(spec.get("prune", 0.0), spec["slope"], spec["curve"]))
    if spec["kind"] == "fixed":
        f.set_initial_guess()
        z, T = np.asarray(spec["profile_z"], float), np.asarray(spec["profile_T"], float)
        f.flame.set_fixed_temp_profile(z / z[-1], T)
        f.energy_enabled = False
    elif spec["kind"] != "burner":
        raise ValueError("unsupported case kind " + spec["kind"])
    f.solve(loglevel=0, auto=(spec["kind"] == "burner"))
    ordinary = f.solver_stats
    jacobians = sum(ordinary.get("jacobian_evals", []))
    steps = sum(ordinary.get("steps", []))
    regrids = max(0, len(ordinary.get("grid_points", [])) - 1)
    tight_rtol = None
    try:
        f.max_time_step_count = 5000
    except Exception:
        pass
    for rtol in (1e-6, 1e-7, 1e-8):
        saved = f.to_array()
        for domain in f.domains:
            try:
                domain.set_steady_tolerances(default=(rtol, 1e-14))
            except Exception:
                pass
        try:
            f.solve(loglevel=0, refine_grid=False)
            tight_rtol = rtol
        except Exception:
            f.from_array(saved)
            break
    names = "\n".join(gas.species_names).encode("utf-8")
    np.savez_compressed(out, grid=np.asarray(f.grid), T=np.asarray(f.T), Y=np.asarray(f.Y),
        inlet_Y=inlet_y, species_names_utf8=np.frombuffer(names, dtype=np.uint8),
        jacobians=np.array([jacobians]), time_steps=np.array([steps]),
        regrids=np.array([regrids]), tight_rtol=np.array([tight_rtol or 0.0]))
    print("CT_TIGHT case={} kind={} points={} jacobians={} steps={} regrids={} tight_rtol={}".format(
        spec["id"], spec["kind"], len(f.grid), jacobians, steps, regrids, tight_rtol), flush=True)


args = sys.argv[1:]
if len(args) < 2 or len(args) % 2:
    raise SystemExit("usage: fixed_t_six_ct.py spec.json output.npz [spec.json output.npz ...]")
for i in range(0, len(args), 2):
    run_case(args[i], args[i + 1])
