#!/usr/bin/env python3
"""Capture Cantera initial, successful time-step, steady, and grid stages."""
import json, os, sys, time
from pathlib import Path
import numpy as np
import cantera as ct

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ct_dump import build
from jacobian_trace import cantera_attempt_trace

def snapshot(f, root, idx, event, dt=None, kind=None):
    # Capture every domain, including boundary state, so replay inputs preserve
    # the full Sim1D vector and domain layout.
    arrays = {}
    layout = {}
    for d in f.domains:
        layout[d.name] = {"components": [d.component_name(i) for i in range(d.n_components)],
                          "n_points": d.n_points}
        for i in range(d.n_components):
            name = d.component_name(i)
            arrays[f"{d.name}/{name}"] = np.asarray(d.values(name)).copy()
        arrays[f"{d.name}/grid"] = np.asarray(d.grid).copy()
    if hasattr(f, "density"):
        arrays["flame/density"] = np.asarray(f.density).copy()
    # Canonical aliases consumed by the shared xresid harness.
    arrays["grid"] = arrays["flame/grid"]
    for name in layout["flame"]["components"]:
        if name in ("eField", "Uo"):
            continue
        if name in ("spreadRate", "Lambda") and kind not in ("diffusion", "premixed", "twin", "stagnation"):
            continue
        arrays[f"values/{name}"] = arrays[f"flame/{name}"]
    tag = f"{idx:04d}" if isinstance(idx, int) else str(idx)
    path = root / f"stage_{tag}.npz"
    np.savez_compressed(path, **arrays)
    return {"index": idx, "event": event, "dt": dt, "file": path.name,
            "wall_seconds": time.perf_counter(), "grid": arrays["flame/grid"].tolist(),
            "layout": layout, "solver_stats": f.solver_stats}

def decision_trace(spec, out):
    """Record the fixed initial-grid decision trace from the captured guess."""
    f = build(spec)
    if spec["kind"] != "stagnation": f.set_initial_guess()
    f.energy_enabled = True
    trace_initial = snapshot(f, out, "decision_trace_initial", "decision_trace_initial", kind=spec["kind"])
    captured = np.load(out / "stage_0000.npz")
    traced = np.load(out / trace_initial["file"])
    common = sorted(set(captured.files) & set(traced.files))
    exact_initial = common == sorted(set(captured.files)) == sorted(set(traced.files)) and all(
        np.array_equal(captured[key], traced[key]) for key in common)
    initial_max_abs_delta = max((float(np.max(np.abs(captured[key] - traced[key])))
                                 for key in common if captured[key].dtype.kind in "fiu"), default=0.0)
    logpath = out / "ct_decisions.log"
    logfd = os.open(logpath, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    oldfd = os.dup(1)
    error = None
    os.dup2(logfd, 1)
    try:
        f.solve(loglevel=4, refine_grid=False, auto=False)
    except Exception as exc:
        error = str(exc).splitlines()[0]
    finally:
        # Leave stdout redirected until process exit: C++ streams may flush
        # after solve() returns. Write the concise status to the saved fd.
        os.close(logfd)
    (out / "ct_decisions.json").write_text(json.dumps({
        "error": error, "solver_stats": f.solver_stats,
        "grid": np.asarray(f.flame.grid).tolist()}, indent=2))
    os.write(oldfd, f"Cantera fixed-grid trace error={error or 'none'} "
             f"log={os.path.getsize(logpath)} bytes\n".encode())
    os.close(oldfd)
    # Parse rejected and accepted transient attempts. Refresh markers can be
    # emitted before the cumulative N_jac increment appears in the next row;
    # the shared parser assigns that evaluation to the triggering retry.
    lines = logpath.read_text(errors="replace").splitlines()
    attempts = cantera_attempt_trace(logpath)
    (out / "ct_decisions.json").write_text(json.dumps({"error": error,
        "initial_snapshot_equal": exact_initial,
        "initial_snapshot_max_abs_delta": initial_max_abs_delta,
        "initial_steady_failed": "Newton steady-state solve failed" in "\n".join(lines),
        "attempts": attempts, "solver_stats": f.solver_stats}, indent=2))

def main(spec_path, out_path):
    spec = json.loads(Path(spec_path).read_text())
    out = Path(out_path); out.mkdir(parents=True, exist_ok=True)
    f = build(spec)
    if spec["kind"] != "stagnation": f.set_initial_guess()
    f.energy_enabled = True
    stages = []; tic = time.perf_counter()
    stages.append(snapshot(f, out, len(stages), "initial_guess", kind=spec["kind"]))
    previous_grid = [np.asarray(f.flame.grid).tolist()]
    prior_jacobians = [0]
    def on_step(dt):
        jacobians = sum(f.solver_stats.get("jacobian_evals", []))
        # Cantera calls back after adapting dt: it grows dt only when this
        # successful step evaluated no new Jacobian (default growth=1.5).
        used_dt = float(dt) / (1.5 if jacobians == prior_jacobians[0] else 1.0)
        stage = snapshot(f, out, len(stages), "time_step", used_dt, spec["kind"])
        stage["dt_next"] = float(dt)
        stage["jacobians_total"] = jacobians
        stages.append(stage)
        prior_jacobians[0] = jacobians
        return 0.0
    def on_steady(_):
        current = np.asarray(f.flame.grid).tolist()
        stage = snapshot(f, out, len(stages), "steady_solve", kind=spec["kind"])
        stage["grid_before"] = previous_grid[0]
        stage["refined_since_previous"] = current != previous_grid[0]
        stages.append(stage)
        previous_grid[0] = current
        return 0.0
    f.set_time_step_callback(on_step)
    f.set_steady_callback(on_steady)
    error = None
    try: f.solve(loglevel=0, auto=True)
    except Exception as e: error = str(e).splitlines()[0]
    # A changed grid after the last pre-refine steady callback marks a refine.
    if stages[-1]["grid"] != np.asarray(f.flame.grid).tolist():
        stages.append(snapshot(f, out, len(stages), "final_grid", kind=spec["kind"]))
    previous_time = tic
    for s in stages:
        captured_at = s.pop("wall_seconds")
        s["elapsed_seconds"] = captured_at - tic
        s["stage_wall_seconds"] = captured_at - previous_time
        previous_time = captured_at
    (out / "manifest.json").write_text(json.dumps({"case": spec_path, "kind": spec["kind"],
        "converged": error is None, "error": error, "total_seconds": time.perf_counter()-tic,
        "stages": stages, "solver_stats": f.solver_stats}, indent=2))
    print(f"captured stages={len(stages)} converged={error is None} seconds={time.perf_counter()-tic:.3f}")
    decision_trace(spec, out)

if __name__ == "__main__": main(sys.argv[1], sys.argv[2])
