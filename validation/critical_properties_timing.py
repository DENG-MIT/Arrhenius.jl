"""Check every property in the complete eight-fluid critical-properties example.

Prepared pure-fluid models; reading Tc/Pc/rhoc/MW, calculating Zc and collecting
the numerical table are timed. Imports, model loading, plots and printing are
excluded. All 40 output values are checked after every timed calculation.
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import math
import statistics
import time
from benchmark_environment import host_metadata, matches_target, cantera_library_hashes, verify_numerical_threads

started = time.perf_counter()
import numpy as np
import cantera as ct
import cantera._cantera as compiled
IMPORT_SECONDS = time.perf_counter()-started
COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
SOURCE_SHA256 = "fa710297d699c5de3af42b85799518fe752694871020a09d383441e86f8219b6"
NAMES = ["water", "nitrogen", "methane", "hydrogen", "oxygen", "carbon dioxide", "heptane", "HFC-134a"]
PROPERTIES = ["critical_temperature", "critical_pressure", "critical_density", "mean_molecular_weight", "critical_compressibility"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculation(fluids):
    output = np.empty((5, len(fluids)))
    for column, name in enumerate(fluids):
        fluid = fluids[name]
        tc, pc = fluid.critical_temperature, fluid.critical_pressure
        rc, mw = fluid.critical_density, fluid.mean_molecular_weight
        zc = pc*mw/(rc*ct.gas_constant*tc)
        output[:, column] = tc, pc, rc, mw, zc
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
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
        parser.error("source example differs from the pinned calculation")
    native = np.load(args.julia_result)
    meta = {k[:-5]: bytes(native[k]).decode() for k in native.files if k.endswith("_utf8")}
    if meta["fluid_names"].splitlines() != NAMES or meta["property_names"].splitlines() != PROPERTIES:
        parser.error("native property or fluid ordering differs")
    for key, digest in meta.items():
        if key.endswith(".jl_sha256"):
            path = Path(__file__).resolve().parent.parent/key.removesuffix("_sha256")
            if sha(path) != digest:
                parser.error(f"native source changed after measurement: {path}")
    if meta["harness_sha256"] != sha(Path(__file__).with_suffix(".jl")):
        parser.error("native timing driver differs")
    build = json.loads(args.cantera_build_record.read_text())
    libraries = cantera_library_hashes(ct.__file__)
    recorded = {Path(k).name: v for k, v in build["library_hashes"].items()}
    build_matches = bool(libraries and all(recorded.get(k) == v for k, v in libraries.items())
                         and recorded.get(Path(compiled.__file__).name) == sha(compiled.__file__))
    if not build_matches or build["source"]["commit"] != COMMIT or not ct.__version__.startswith("4.0"):
        parser.error("loaded Cantera libraries differ from the pinned development build")
    threads_before = verify_numerical_threads(set_accelerate=True)
    hardware = host_metadata()
    started = time.perf_counter()
    fluids = dict(zip(NAMES, [ct.Water(), ct.Nitrogen(), ct.Methane(), ct.Hydrogen(),
                             ct.Oxygen(), ct.CarbonDioxide(), ct.Heptane(), ct.Hfc134a()]))
    construction_seconds = time.perf_counter()-started
    started = time.perf_counter()
    output = calculation(fluids)
    first_seconds = time.perf_counter()-started
    samples, batch_size = [], 0
    if args.qualification != "validate-only":
        started = time.perf_counter()
        warmup = calculation(fluids)
        warmup_seconds = time.perf_counter()-started
        if not np.array_equal(warmup, output):
            raise RuntimeError("Cantera warmup changed output")
        batch_size = max(1, min(1000, math.ceil(.020/max(warmup_seconds, 1e-9))))
        gc.collect()
        for _ in range(args.repetitions):
            elapsed = 0.
            for _ in range(batch_size):
                started = time.perf_counter()
                repeated = calculation(fluids)
                elapsed += time.perf_counter()-started
                if not np.array_equal(output, repeated):
                    raise RuntimeError("Cantera repetition changed output")
            samples.append(elapsed/batch_size)
            verify_numerical_threads()
    threads_after = verify_numerical_threads()
    candidate = native["properties"]
    if candidate.shape != output.shape or candidate.shape != (5, 8):
        raise RuntimeError("native result must contain all five properties of all eight fluids")
    checks = {key: {"pass": bool(np.all(np.isfinite(candidate[i])) and np.allclose(candidate[i], output[i], rtol=2e-14, atol=0)),
                    "maximum_relative_error": float(np.max(np.abs(candidate[i]/output[i]-1))),
                    "relative_tolerance": 2e-14, "absolute_tolerance": 0,
                    "native": candidate[i].tolist(), "cantera": output[i].tolist()}
              for i, key in enumerate(PROPERTIES)}
    native_samples = native["seconds"].tolist()
    native_batch = int(native["batch_size"][0])
    native_checked = int(native["warm_outputs_checked"][0])
    if native_checked != len(native_samples)*native_batch:
        raise RuntimeError("inconsistent native replay count")
    if any(not math.isfinite(v) or v <= 0 for v in samples+native_samples):
        raise RuntimeError("invalid elapsed time")
    thread_names = meta["numerical_thread_names"].splitlines()
    native_threads = dict(zip(thread_names, native["numerical_threads"].tolist(), strict=True))
    threads_match = (native_threads.get("julia_threads") == native_threads.get("blas_threads") == 1
                     and all(v == 1 for v in native_threads.values())
                     and (args.target != "apple-m4" or native_threads.get("accelerate_threading_mode") == 1))
    ratio = statistics.median(samples)/statistics.median(native_samples) if samples and native_samples else None
    correct = all(check["pass"] for check in checks.values())
    controlled = bool(args.qualification == meta["qualification"] == "controlled"
                      and matches_target(hardware, args.target) and matches_target(meta, args.target)
                      and all(hardware[k] == meta[k] for k in ("cpu", "kernel_release"))
                      and threads_match and len(samples) >= 9 and len(native_samples) >= 9
                      and 1 <= native_batch <= 1000)
    report = {"example": "thermo/critical_properties", "cantera_version": ct.__version__, "cantera_source_sha": COMMIT,
              "published_source_sha256": SOURCE_SHA256, "cantera_build_record_sha256": sha(args.cantera_build_record),
              "loaded_libraries_match_build_record": build_matches, "cantera_shared_libraries_sha256": libraries,
              "cantera_extension_sha256": sha(compiled.__file__), "harness_sha256": sha(__file__),
              "environment_helper_sha256": sha(Path(__file__).with_name("benchmark_environment.py")),
              "native_metadata": meta, "native_artifact_sha256": sha(args.julia_result), "hardware": hardware,
              "benchmark_target": args.target, "numerical_threads": threads_before, "numerical_threads_after": threads_after,
              "native_numerical_threads": native_threads, "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "scope": "complete eight-fluid table on prepared models: Tc, Pc, critical density, molecular weight and calculated Zc; plots, formatting, imports and model construction excluded",
              "checks": checks, "correctness_pass": correct, "fluid_names": NAMES,
              "cantera_import_seconds": IMPORT_SECONDS, "cantera_model_construction_seconds": construction_seconds,
              "julia_import_seconds": float(native["import_seconds"][0]),
              "cantera_first_seconds": first_seconds, "julia_first_seconds": float(native["first_seconds"][0]),
              "cantera_warm_seconds": samples, "julia_warm_seconds": native_samples,
              "sample_definition": "mean timed seconds per complete calculation within each warm batch; checks excluded",
              "batch_policy": "ceil(0.020 / warmup seconds), clamped to 1..1000 complete calculations",
              "cantera_batch_size": batch_size, "julia_batch_size": native_batch,
              "warm_outputs_checked": len(samples)*batch_size, "native_warm_outputs_checked": native_checked,
              "speed_ratio": ratio, "minimum_speed_ratio": .95, "qualification": "controlled" if controlled else "not_qualified",
              "performance_pass": bool(controlled and correct and ratio >= .95)}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
    print("critical properties", "correctness", correct, "speed ratio", ratio, "controlled", controlled, flush=True)
    if not correct:
        raise SystemExit("native critical properties failed reference checks")


if __name__ == "__main__":
    main()
