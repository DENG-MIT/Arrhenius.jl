#!/usr/bin/env python3
"""Bounded EEDF timing harness: 10 complete fresh ct.Solution calculations.

Diagnostic only. Emits timings.npz + receipt.json and applies NO speed or
physics qualification; the parent applies speed>=0.95 and independent
physics gates separately.
Provenance: cantera source commit 726522be4e2a13454d8415b7ef799d621f665cf3.
Usage: ct_eedf_timing.py CONFIG_JSON OUTPUT_DIRECTORY
"""
import hashlib, json, os, platform, sys, time

import numpy as np

COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
T, P = 300., 101325.
XSTR = "N2:0.79,O2:0.21,N2+:1E-10,Electron:1E-10"
SPECIES = {"N2": 0.79, "O2": 0.21, "N2+": 1e-10, "Electron": 1e-10}
XSUM = sum(SPECIES.values())
EN = 200e-21  # reduced electric field, V*m^2
NREP, NCELL = 10, 41
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fresh_run(ct, model):
    """One COMPLETE FRESH calculation; nothing cached across calls."""
    gas = ct.Solution(str(model))
    gas.TPX = T, P, XSTR
    gas.reduced_electric_field = EN
    gas.update_electron_energy_distribution()
    grid = gas.electron_energy_levels.copy()
    eedf = gas.electron_energy_distribution.copy()
    return gas, grid, eedf


def check(gas, grid, eedf):
    if (grid.shape != (NCELL,) or not np.all(np.isfinite(grid))
            or not np.array_equal(grid, np.arange(NCELL, dtype=float))):
        raise ValueError("energy grid must be 41 finite edges exactly np.arange(41.)")
    if (eedf.shape != (NCELL,) or not np.all(np.isfinite(eedf))
            or (eedf < 0).any() or not (eedf > 0).any()):
        raise ValueError("EEDF must be 41 finite, nonnegative, with positive values")
    if gas.T != T or gas.P != P or gas.reduced_electric_field != EN:
        raise ValueError("T/P/reduced_electric_field not as prescribed")
    for k in gas.species_names:
        v = SPECIES.get(k, 0.0)
        if not np.isclose(gas[k].X[0], v / XSUM, rtol=1e-9, atol=0.):
            raise ValueError(f"mole fraction mismatch for {k}")


def write_artifacts(outdir, receipt, times, grids, eedfs, xs):
    np.savez(os.path.join(outdir, "timings.npz"),
             times=np.asarray(times, float),
             grids=np.asarray(grids, float).reshape(-1, NCELL),
             eedfs=np.asarray(eedfs, float).reshape(-1, NCELL),
             mole_fractions=np.asarray(xs, float), n_completed=len(grids))
    with open(os.path.join(outdir, "receipt.json"), "w") as fh:
        json.dump(receipt, fh, indent=2, sort_keys=True)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: ct_eedf_timing.py CONFIG_JSON OUTPUT_DIRECTORY")
    cfg_path, outdir = sys.argv[1], sys.argv[2]
    with open(cfg_path) as fh:
        cfg = json.load(fh)
    model = cfg.get("model_path") or cfg.get("model")
    receipt = {
        "schema": "eedf-timing-receipt/v1",
        "provenance": {"cantera_source_commit": COMMIT},
        "argv": sys.argv, "platform": platform.platform(),
        "python_version": sys.version.split()[0],
        "thread_env": {k: os.environ.get(k) for k in THREAD_VARS},
        "qualification": {"diagnostic": True, "final": False},
        "conditions": {"T":T,"P":P,"EN":EN,"X":{k:v/XSUM for k,v in SPECIES.items()}},
        "script_sha256": sha256(os.path.abspath(__file__)),
    }
    bad = [k for k in THREAD_VARS if os.environ.get(k) != "1"]
    if bad:
        raise SystemExit(f"thread env vars must all be '1': {bad}")
    if os.path.exists(outdir):
        raise SystemExit(f"output directory must be fresh/nonexistent: {outdir}")
    os.makedirs(outdir)
    times, grids, eedfs, xs = [], [], [], []
    try:
        if sha256(model) != cfg["expected_model_sha256"]:
            raise ValueError("model sha256 mismatch")
        receipt["model_path"] = model
        receipt["model_sha256"] = cfg["expected_model_sha256"]
        import cantera as ct
        if ct.__version__ != cfg["expected_cantera_version"]:
            raise ValueError(f"cantera version mismatch: {ct.__version__}")
        receipt["cantera_version"] = ct.__version__
        ext = ct._cantera.__file__
        if sha256(ext) != cfg["expected_extension_sha256"]:
            raise ValueError("extension sha256 mismatch")
        receipt["extension_path"] = ext
        receipt["extension_sha256"] = cfg["expected_extension_sha256"]
        for rep in range(NREP):  # rep 0 = first call; reps 1..9 warm
            t0 = time.perf_counter()
            gas, grid, eedf = fresh_run(ct, model)
            times.append(time.perf_counter() - t0)
            grids.append(grid); eedfs.append(eedf); xs.append(gas.X.copy())
            check(gas, grid, eedf)  # validation/logging outside the timer
            del gas # object destruction is outside the complete-calculation timer
        if not all(np.array_equal(grids[0], g) for g in grids[1:]):
            raise ValueError("energy grid differs across repetitions")
        if not all(np.array_equal(xs[0], x) for x in xs[1:]):
            raise ValueError("composition differs across repetitions")
        if (sha256(model) != receipt["model_sha256"]
                or sha256(ext) != receipt["extension_sha256"]):
            raise ValueError("model/extension hash changed during run")
        receipt.update(raw_times=times, n_repetitions=NREP,
                       first_call_time=times[0],
                       median_warm_time=float(np.median(times[1:])),
                       hashes_unchanged_after_run=True, status="completed")
        write_artifacts(outdir, receipt, times, grids, eedfs, xs)
    except BaseException as exc:
        receipt.update(status="failed", raw_times=times,
                       failure={"type": type(exc).__name__, "message": str(exc),
                                "n_completed": len(times)})
        try:
            write_artifacts(outdir, receipt, times, grids, eedfs, xs)
        finally:
            raise


if __name__ == "__main__":
    main()
