#!/usr/bin/env python3
"""Generate one gate-ready Cantera tight archive from a case spec.

The archive records the effective steady tolerances and the successful return
from Cantera 4.0.0a2. Run from the NH3 harness directory so mechanism files
resolve; output belongs in a task-owned directory, never a shared refs folder.
"""
import json
import os
import sys
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import cantera as ct
import numpy as np

RTOL = 1.0e-8
ATOL = 1.0e-14


def build(spec):
    gas = ct.Solution(spec["mechanism"] + ".yaml")
    kind = spec["kind"]
    P = float(spec["P_atm"]) * ct.one_atm
    if kind == "diffusion":
        gas.TP = spec["T"], P
        f = ct.CounterflowDiffusionFlame(gas, width=spec["width"])
        f.fuel_inlet.mdot, f.fuel_inlet.X, f.fuel_inlet.T = (
            spec["mdot_fuel"], spec["fuel"], spec["T"])
        f.oxidizer_inlet.mdot, f.oxidizer_inlet.X, f.oxidizer_inlet.T = (
            spec["mdot_oxidizer"], spec["oxidizer"], spec["T"])
    else:
        gas.TPX = spec["T"], P, spec.get("X", "N2:1")
        if kind == "premixed":
            f = ct.CounterflowPremixedFlame(gas, width=spec["width"])
            f.reactants.mdot = spec["mdot"]
            gas.equilibrate("HP")
            f.products.X, f.products.T = gas.X, gas.T
            f.products.mdot = spec["mdot_products"]
        elif kind == "twin":
            f = ct.CounterflowTwinPremixedFlame(gas, width=spec["width"])
            f.reactants.mdot = spec["mdot"]
        elif kind == "stagnation":
            f = ct.ImpingingJet(gas, width=spec["width"])
            f.inlet.mdot = spec["mdot"]
            f.surface.T = spec["T_surface"]
            f.set_initial_guess("equil")
        elif kind == "free":
            f = ct.FreeFlame(gas, width=spec["width"])
        elif kind in ("burner", "fixed"):
            f = ct.BurnerFlame(gas, width=spec["width"])
            f.burner.mdot = spec["mdot"]
            if kind == "fixed":
                z, T = np.asarray(spec["profile_z"]), np.asarray(spec["profile_T"])
                f.flame.set_fixed_temp_profile(z / z[-1], T)
                f.energy_enabled = False
        else:
            raise ValueError(f"unsupported case kind: {kind}")
    f.transport_model = spec["transport"]
    return f


def main(spec_path, output_path):
    spec = json.loads(Path(spec_path).read_text())
    flame = build(spec)
    for domain in flame.domains:
        domain.set_steady_tolerances(default=(RTOL, ATOL))
    effective_rtol = [np.asarray(d.steady_reltol()) for d in flame.domains]
    effective_atol = [np.asarray(d.steady_abstol()) for d in flame.domains]
    if not all(np.all(x == RTOL) for x in effective_rtol):
        raise RuntimeError("effective Cantera steady relative tolerances differ from requested value")
    if not all(np.all(x == ATOL) for x in effective_atol):
        raise RuntimeError("effective Cantera steady absolute tolerances differ from requested value")
    prune = spec.get("prune", 0.0)
    slope, curve = float(spec["slope"]), float(spec["curve"])
    flame.set_refine_criteria(ratio=spec["ratio"], slope=slope, curve=curve,
                              prune=min(prune, slope, curve))
    # Fixed-temperature flames must retain their prescribed profile.
    flame.solve(loglevel=0, auto=(spec["kind"] != "fixed"))
    stats = flame.solver_stats
    data = dict(
        grid=flame.grid.copy(), T=flame.T.copy(), Y=flame.Y.copy(),
        velocity=flame.velocity.copy(),
        SL=np.asarray(float(flame.velocity[0]) if spec["kind"] == "free"
                       else float(np.max(flame.T))),
        mdot=np.asarray(float(flame.density[0] * flame.velocity[0])),
        species_names=np.asarray(flame.gas.species_names),
        newton_ok=np.asarray([1], dtype=np.int8),
        solver_rtol=np.asarray([RTOL]), solver_atol=np.asarray([ATOL]),
    )
    for name, key in (("jacobians", "jacobian_evals"),
                      ("newton_iterations", "linear_solves"),
                      ("time_steps", "steps")):
        if key not in stats:
            raise RuntimeError(f"Cantera solver_stats omitted required count {key}")
        data[f"stats_{name}"] = np.asarray([sum(stats[key])], dtype=np.int64)
    if "grid_points" not in stats:
        raise RuntimeError("Cantera solver_stats omitted required count grid_points")
    data["stats_regrids"] = np.asarray([max(0, len(stats["grid_points"]) - 1)], dtype=np.int64)
    data["stats_final_points"] = np.asarray([len(flame.grid)], dtype=np.int64)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **data)
    print(f"tight {spec['id']} newton_ok=1 rtol={RTOL:.0e} atol={ATOL:.0e} "
          f"Tmax={max(flame.T):.6g} points={len(flame.grid)} "
          f"jac={int(data['stats_jacobians'][0])} steps={int(data['stats_time_steps'][0])} "
          f"out={output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: ct_tight_reference.py SPEC.json OUT.ct.npz")
    main(sys.argv[1], sys.argv[2])



