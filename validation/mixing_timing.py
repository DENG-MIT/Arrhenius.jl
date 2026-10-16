"""Check and time the complete constant-pressure mixing example with Cantera 4.

Run mixing_timing.jl first. Both sides use prepared mechanisms, create the two
streams, mix at HP, equilibrate at TP, and collect both state/property reports.
Imports, mechanism loading, text formatting, resets and checks are not timed.
Nine batches of complete warm calculations give normalized per-call samples.
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import math
import statistics
import time
import os
from benchmark_environment import host_metadata, matches_target, cantera_library_hashes, verify_numerical_threads

started = time.perf_counter()
import numpy as np
import cantera as ct
import cantera._cantera as compiled
IMPORT_SECONDS = time.perf_counter()-started
COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
SOURCE_SHA256 = "cfbdefb42f135d7d17915bf5bd5aa0ef870588dc8c2f8f678d0c6145a854c362"
SCALARS = ("T", "P", "density", "mean_molecular_weight", "enthalpy_mole",
           "int_energy_mole", "entropy_mole", "gibbs_mole", "cp_mole", "cv_mole",
           "moles", "mass", "enthalpy")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def report_state(quantity):
    return {"scalars": np.array([getattr(quantity, key) for key in SCALARS]),
            "X": quantity.X.copy(), "Y": quantity.Y.copy(),
            "mu_RT": quantity.chemical_potentials/(ct.gas_constant*quantity.T)}


def calculation(gas):
    air = ct.Quantity(gas, constant="HP")
    air.TPX = 300., ct.one_atm, "O2:0.21,N2:0.78,AR:0.01"
    methane = ct.Quantity(gas, constant="HP")
    methane.TPX = 300., ct.one_atm, "CH4:1"
    air.moles = 1.
    methane.moles = air.X[air.species_index("O2")]*.5
    mixed = air+ methane
    before = report_state(mixed)
    mixed.equilibrate("TP")
    after = report_state(mixed)
    return {key: np.column_stack([before[key], after[key]]) for key in before}


def identical(left, right):
    return left.keys() == right.keys() and all(np.array_equal(left[k], right[k]) for k in left)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--source-example", type=Path, required=True)
    parser.add_argument("--julia-result", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cantera-build-record", type=Path, required=True)
    parser.add_argument("--target", choices=("wsl", "apple-m4"), required=True)
    parser.add_argument("--qualification", choices=("controlled", "informational", "validate-only"), default="informational")
    parser.add_argument("--repetitions", type=int, default=9)
    args = parser.parse_args()
    if args.repetitions < 9:
        parser.error("at least nine warm batches required")
    if sha(args.source_example) != SOURCE_SHA256:
        parser.error("source example differs from the pinned Cantera calculation")
    native = np.load(args.julia_result)
    meta = {k[:-5]: bytes(native[k]).decode() for k in native.files if k.endswith("_utf8")}
    if meta["mechanism_sha256"] != sha(args.mechanism):
        parser.error("native and reference mechanisms differ")
    if meta["scalar_names"].splitlines() != list(SCALARS):
        parser.error("native property ordering differs")
    gas = ct.Solution(str(args.mechanism))
    if gas.n_species != 53 or meta["species_names"].splitlines() != gas.species_names:
        parser.error("expected matching 53-species GRI-Mech phase")
    sidecar = Path(str(args.mechanism)+".npz")
    if meta["sidecar_sha256"] != sha(sidecar):
        parser.error("native mechanism sidecar differs")
    sidecar_values = np.load(sidecar)
    if bytes(sidecar_values["source_sha256_utf8"]).decode() != sha(args.mechanism):
        parser.error("native sidecar was exported from a different mechanism")
    build = json.loads(args.cantera_build_record.read_text())
    libraries = cantera_library_hashes(ct.__file__)
    recorded = {Path(k).name: v for k, v in build["library_hashes"].items()}
    build_matches = bool(libraries and all(recorded.get(k) == v for k, v in libraries.items())
                         and recorded.get(Path(compiled.__file__).name) == sha(compiled.__file__))
    if not build_matches or build["source"]["commit"] != COMMIT or not ct.__version__.startswith("4.0"):
        parser.error("loaded Cantera libraries do not match the pinned development build")
    threads_before = verify_numerical_threads(set_accelerate=True)
    hardware = host_metadata()
    reset = lambda: setattr(gas, "TPX", (300., ct.one_atm, "CH4:1"))
    reset()
    start = time.perf_counter()
    output = calculation(gas)
    first_seconds = time.perf_counter()-start
    samples, batch_size = [], 0
    if args.qualification != "validate-only":
        reset()
        start = time.perf_counter()
        warmup = calculation(gas)
        warmup_seconds = time.perf_counter()-start
        if not identical(output, warmup):
            raise RuntimeError("Cantera warmup changed output")
        batch_size = max(1, min(1000, math.ceil(.020/max(warmup_seconds, 1e-9))))
        gc.collect()
        for _ in range(args.repetitions):
            elapsed = 0.
            for _ in range(batch_size):
                reset()
                start = time.perf_counter()
                repeated = calculation(gas)
                elapsed += time.perf_counter()-start
                if not identical(output, repeated):
                    raise RuntimeError("Cantera repeated calculation changed output")
            samples.append(elapsed/batch_size)
    threads_after = verify_numerical_threads()
    checks = {}
    for key, reference in output.items():
        candidate = native[key]
        if candidate.shape != reference.shape:
            raise RuntimeError(f"{key}: native output dimensions differ")
        rtol, atol = (2e-8, 2e-7) if key in ("scalars", "mu_RT") else (2e-8, 2e-12)
        # Quantity.report displays chemical potentials only for non-minor
        # species. Zero/trace fractions use different log regularizations in
        # the two libraries; all mole/mass fractions are checked above.
        mask = output["X"] >= 1e-14 if key == "mu_RT" else np.ones(reference.shape, dtype=bool)
        checks[key] = {"pass": bool(np.all(np.isfinite(candidate)) and np.allclose(candidate[mask], reference[mask], rtol=rtol, atol=atol)),
                       "maximum_absolute_error": float(np.max(np.abs(candidate[mask]-reference[mask]))),
                       "relative_tolerance": rtol, "absolute_tolerance": atol,
                       "native": candidate.tolist(), "cantera": reference.tolist()}
        if key == "mu_RT":
            checks[key]["scope"] = "chemical potentials of every species displayed individually by the source report (X >= 1e-14); full X/Y arrays checked separately"
    # Independent frozen-mixing enthalpy and mass closure; TP chemistry changes
    # enthalpy but conserves atoms and total mass.
    air = ct.Solution(str(args.mechanism))
    fuel = ct.Solution(str(args.mechanism))
    air.TPX, fuel.TPX = (300., ct.one_atm, "O2:.21,N2:.78,AR:.01"), (300., ct.one_atm, "CH4:1")
    expected_mass = air.mean_molecular_weight+.105*fuel.mean_molecular_weight
    expected_enthalpy = air.enthalpy_mole+.105*fuel.enthalpy_mole
    s = native["scalars"]
    atoms = np.array([[gas.n_atoms(k, e) for k in range(gas.n_species)] for e in range(gas.n_elements)])
    element_amounts = atoms @ (native["X"]*s[10])
    checks["conservation"] = {
        "pass": bool(np.allclose(s[11], expected_mass, rtol=1e-12, atol=1e-12)
                     and abs(s[12,0]-expected_enthalpy) <= 1e-10*max(abs(expected_enthalpy), 1.)
                     and np.allclose(element_amounts[:,0], element_amounts[:,1], rtol=1e-10, atol=1e-12)),
        "mass_error_kg": float(np.max(np.abs(s[11]-expected_mass))),
        "mixing_enthalpy_error_J": float(abs(s[12,0]-expected_enthalpy)),
        "element_amount_error_kmol": float(np.max(np.abs(element_amounts[:,0]-element_amounts[:,1])))}
    native_samples = native["seconds"].tolist()
    native_batch = int(native["batch_size"][0])
    native_checked = int(native["warm_outputs_checked"][0])
    if native_checked != len(native_samples)*native_batch:
        raise RuntimeError("inconsistent native replay count")
    if any(not math.isfinite(v) or v <= 0 for v in samples+native_samples):
        raise RuntimeError("invalid elapsed time")
    ratio = statistics.median(samples)/statistics.median(native_samples) if samples and native_samples else None
    correct = all(check["pass"] for check in checks.values())
    controlled = bool(args.qualification == meta["qualification"] == "controlled"
                      and matches_target(hardware, args.target) and matches_target(meta, args.target)
                      and all(hardware[k] == meta[k] for k in ("cpu", "kernel_release"))
                      and int(native["julia_threads"][0]) == int(native["blas_threads"][0]) == 1
                      and len(samples) >= 9 and len(native_samples) >= 9 and 1 <= native_batch <= 1000)
    report = {"example": "thermo/mixing", "cantera_version": ct.__version__, "cantera_source_sha": COMMIT,
              "published_source_sha256": SOURCE_SHA256, "mechanism_sha256": sha(args.mechanism),
              "cantera_build_record_sha256": sha(args.cantera_build_record), "loaded_libraries_match_build_record": build_matches,
              "cantera_shared_libraries_sha256": libraries, "cantera_extension_sha256": sha(compiled.__file__),
              "harness_sha256": sha(__file__), "environment_helper_sha256": sha(Path(__file__).with_name("benchmark_environment.py")),
              "native_metadata": meta, "native_artifact_sha256": sha(args.julia_result), "hardware": hardware,
              "benchmark_target": args.target, "numerical_threads": threads_before, "numerical_threads_after": threads_after,
              "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "scope": "prepared phase; fresh air/methane streams, HP mixing, TP equilibrium and both numerical state reports; no formatting, resets, imports, model loading or I/O",
              "checks": checks, "correctness_pass": correct, "scalar_names": SCALARS,
              "cantera_import_seconds": IMPORT_SECONDS, "julia_import_seconds": float(native["import_seconds"][0]),
              "cantera_first_seconds": first_seconds, "julia_first_seconds": float(native["first_seconds"][0]),
              "cantera_warm_seconds": samples, "julia_warm_seconds": native_samples,
              "sample_definition": "mean timed seconds per complete calculation within each warm batch; resets and checks excluded",
              "batch_policy": "ceil(0.020 / warmup seconds), clamped to 1..1000 complete calculations",
              "cantera_batch_size": batch_size, "julia_batch_size": native_batch,
              "warm_outputs_checked": len(samples)*batch_size, "native_warm_outputs_checked": native_checked,
              "speed_ratio": ratio, "minimum_speed_ratio": .95,
              "qualification": "controlled" if controlled else "not_qualified",
              "performance_pass": bool(controlled and correct and ratio >= .95)}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
    print("mixing", "correctness", correct, "speed ratio", ratio, "controlled", controlled, flush=True)
    if not correct:
        raise SystemExit("native mixing outputs failed independent checks")


if __name__ == "__main__":
    main()
