#!/usr/bin/env python3
"""Capture Cantera's first fixed-temperature burner decisions from its ramp guess."""
import json, os, sys
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "1"
import cantera as ct
import numpy as np
spec = json.load(open(sys.argv[1], encoding="utf-8"))
gas = ct.Solution(spec["mechanism"] + ".yaml")
gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
f = ct.BurnerFlame(gas, width=spec["width"])
f.burner.mdot = spec["mdot"]
f.transport_model = spec["transport"]
f.set_refine_criteria(ratio=spec["ratio"], slope=spec["slope"], curve=spec["curve"], prune=0.0)
f.set_initial_guess()
z, T = np.asarray(spec["profile_z"], float), np.asarray(spec["profile_T"], float)
f.flame.set_fixed_temp_profile(z / z[-1], T)
if len(sys.argv) > 2:
    np.savez_compressed(sys.argv[2], grid=np.asarray(f.grid), T=np.asarray(f.T),
                        Y=np.asarray(f.Y).T, velocity=np.asarray(f.velocity))
f.energy_enabled = False
f.max_time_step_count = 10
print("INITIAL_GUESS", "grid=", np.asarray(f.grid), "T=", np.asarray(f.T), "u=", np.asarray(f.velocity), flush=True)
try:
    f.solve(loglevel=4, refine_grid=False, auto=False)
except ct.CanteraError as exc:
    print("STOP", str(exc).splitlines()[0], flush=True)
print("FINAL", "points=", len(f.grid), "steps=", sum(f.solver_stats.get("steps", [])),
      "jacobians=", sum(f.solver_stats.get("jacobian_evals", [])), flush=True)