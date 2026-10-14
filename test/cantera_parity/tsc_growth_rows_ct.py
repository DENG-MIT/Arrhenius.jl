"""Cantera reference half of the timestep-controls growth-strategy rows: the 8
BurnerFlame rows of flame_solver_timestep_controls.py (gri30 CH4/air phi=0.8,
300 K, 1 atm, mdot 0.09), recording per row the accepted time-step count, the
largest step, Jacobian evaluations, and the converged T_max and burner
velocity. Usage: python3 tsc_growth_rows_ct.py <tag>; writes <tag>_growth_ct.csv.
"""
import sys
import time
from pathlib import Path

import numpy as np
import cantera as ct

TAG = sys.argv[1]
MECH = str(Path(__file__).resolve().parents[2] / "mechanism" / "gri30.yaml")

ROWS = [
    ("fixed-growth", 1.5),
    ("fixed-growth", 2.0),
    ("steady-norm", 2.0),
    ("transient-residual", 2.0),
    ("residual-ratio", 2.0),
    ("newton-iterations", 2.0),
    ("steady-norm", 6.0),
    ("newton-iterations", 6.0),
]

out = open(f"{TAG}_growth_ct.csv", "w")
out.write("row,ok,pts,Tmax,speed,steps,jac,dt_max\n")
print("row,ok,pts,Tmax,speed,steps,jac,dt_max,wall_s")
for strategy, growth in ROWS:
    gas = ct.Solution(MECH)
    gas.set_equivalence_ratio(0.8, "CH4", {"O2": 1.0, "N2": 3.76})
    gas.TP = 300.0, ct.one_atm
    flame = ct.BurnerFlame(gas, width=0.02)
    flame.burner.mdot = 0.09
    flame.set_refine_criteria(ratio=3.0, slope=0.25, curve=0.35, prune=0.05)
    dts = []
    flame.set_time_step_callback(lambda dt: dts.append(float(dt)) or 0)
    flame.time_step_growth_strategy = strategy
    flame.time_step_growth_factor = growth
    tic = time.perf_counter()
    try:
        flame.solve(loglevel=0, refine_grid=True)
        ok = True
    except ct.CanteraError as err:
        ok = False
        print(f"{strategy} failed: {err}", file=sys.stderr)
    wall = time.perf_counter() - tic
    steps = len(dts)
    jac = int(sum(flame.solver_stats["jacobian_evals"]))
    tmax = float(np.max(flame.T)) if ok else float("nan")
    speed = float(flame.velocity[0]) if ok else float("nan")
    dt_max = max(dts) if dts else 0.0
    label = f"{strategy} g={growth}"
    print(f"{label},{ok},{len(flame.grid)},{tmax:.2f},{speed:.6f},"
          f"{steps},{jac},{dt_max:.3e},{wall:.1f}")
    out.write(f"{label},{ok},{len(flame.grid)},{tmax!r},{speed!r},"
              f"{steps},{jac},{dt_max!r}\n")
    out.flush()
out.close()
