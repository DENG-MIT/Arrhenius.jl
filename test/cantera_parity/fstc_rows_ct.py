"""Cantera reference half of the fstc counterflow rows: runs the 7 distinct
counterflow cases of flame_solver_timestep_controls.py (no burner section, no
plots), records the time-step callback dt sequence per case, the per-attempt
solver_stats (grid_points/steps), and prints the same summary row as the
native harness. Usage: python3 fstc_rows_ct.py <tag>; writes <tag>_dts_ct.csv.
"""
import os
import sys
import time
import numpy as np
import cantera as ct

TAG = sys.argv[1]

ROWS = [
    ("7MPa regrid=0",            7e6, 0,  200,  "fixed-growth",     1.5),
    ("7MPa regrid=1",            7e6, 1,  200,  "fixed-growth",     1.5),
    ("7MPa regrid=3",            7e6, 3,  200,  "fixed-growth",     1.5),
    ("15MPa budget=200",        15e6, 3,  200,  "fixed-growth",     1.5),
    ("15MPa budget=1250",       15e6, 3,  1250, "fixed-growth",     1.5),
    ("15MPa residual-ratio",    15e6, 3,  1250, "residual-ratio",   1.5),
    ("15MPa newton-iterations", 15e6, 3,  1250, "newton-iterations",1.5),
]

dt_out = open(f"{TAG}_dts_ct.csv", "w")
dt_out.write("row,step,dt\n")

print("row,ok,pts,Tmax,steps,jac,wall_s")
for label, P, regrid_max, budget, strategy, growth in ROWS:
    if os.environ.get("ROW_FILTER") and os.environ["ROW_FILTER"] not in label: continue
    gas = ct.Solution("h2o2.yaml")
    flame = ct.CounterflowDiffusionFlame(
        gas, grid=np.linspace(0.0, 30e-3, 20))
    flame.flame.jacobian_mode = "finite-difference"
    flame.max_time_step_count = budget
    flame.set_refine_criteria(ratio=2.0, slope=0.06, curve=0.08, prune=0.02)
    flame.time_step_regrid = regrid_max
    flame.time_step_growth_strategy = strategy
    flame.time_step_growth_factor = growth
    flame.P = P
    flame.fuel_inlet.X = "H2:1.0"
    flame.fuel_inlet.T = 800.0
    flame.oxidizer_inlet.X = "O2:1.0"
    flame.oxidizer_inlet.T = 711.0
    rho_fuel = flame.fuel_inlet.phase.density
    rho_oxidizer = flame.oxidizer_inlet.phase.density
    flame.fuel_inlet.mdot = 0.3
    flame.oxidizer_inlet.mdot = 0.3 / rho_fuel * rho_oxidizer * 3.0
    flame.clear_stats()

    dts = []
    flame.set_time_step_callback(lambda dt: dts.append(float(dt)) or 0)
    tic = time.perf_counter()
    try:
        flame.solve(loglevel=0, auto=False)
        ok = True
    except ct.CanteraError as err:
        ok = False
        print(f"{label} failed: {err}", file=sys.stderr)
    if ok:
        safe="".join(c if c.isalnum() else "_" for c in label)
        np.savez(f"{TAG}_{safe}.npz",grid=flame.grid,T=flame.T,Y=flame.Y,velocity=flame.velocity,spreadRate=flame.spread_rate,Lambda=flame.L,species_names_utf8=np.frombuffer("\n".join(gas.species_names).encode(),dtype=np.uint8))
    wall = time.perf_counter() - tic
    steps = int(sum(flame.solver_stats["steps"]))
    jac = int(sum(flame.solver_stats["jacobian_evals"]))
    tmax = float(np.max(flame.T))
    print(f"{label},{ok},{len(flame.grid)},{tmax:.1f},{steps},{jac},{wall:.1f}")
    for i in range(len(flame.solver_stats["grid_points"])):
        print(f"  attempt {i+1}: grid={flame.solver_stats['grid_points'][i]} "
              f"steps={flame.solver_stats['steps'][i]} "
              f"jac={flame.solver_stats['jacobian_evals'][i]}")
    for i, dt in enumerate(dts):
        dt_out.write(f"{label},{i+1},{dt!r}\n")
    dt_out.flush()

dt_out.close()
