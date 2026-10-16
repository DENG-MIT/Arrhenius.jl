"""Check and time the complete coverage-dependent surface sweep with Cantera 4.

Run coverage_thermo_timing.jl first. Both sides sweep five prepared interface
phases at 300 K/1 atm: four 101-state CO* enthalpy/entropy curves (linear,
piecewise-linear, polynomial, interpolative) plus the 5151-state triangular
CO*/O* cross-interaction enthalpy map, 5555 states per complete call. Phase
construction, imports, plotting/printing and all checks stay outside every
timer. Nine batches of complete warm sweeps give normalized per-call samples.
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
SOURCE_SHA256 = "b92ddf6180e3798bb5ef3d7475385f94f76f2ff7491bc9ac00eb6521bbef931d"
PHASE_NAMES = ["covdep_lin", "covdep_pwlin", "covdep_poly", "covdep_int", "covdep_cross"]
SHAPES = {"coverages": (101,), "curves": (2, 101, 4), "cross": (101, 101)}
STATES_PER_SWEEP = 4*101+101*102//2


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculation(phases):
    coverages = np.linspace(0., 1., 101)
    curves = np.empty((2, 101, 4))
    for k in range(4):
        phase = phases[k]
        for j, co in enumerate(coverages):
            phase.coverages = [1-co, co]
            curves[0, j, k] = phase.standard_enthalpies_RT[1]
            curves[1, j, k] = phase.standard_entropies_R[1]
    phase = phases[4]
    cross = np.zeros((101, 101))
    for i, co in enumerate(coverages):
        for j in range(101-i):
            oxygen = coverages[j]
            phase.coverages = [1.-co-oxygen, co, oxygen]
            cross[i, j] = phase.standard_enthalpies_RT[1]
    return {"coverages": coverages, "curves": curves, "cross": cross}


def identical(left, right):
    return left.keys() == right.keys() and all(np.array_equal(left[k], right[k]) for k in left)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, required=True)
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
    mechanism_sha256 = sha(args.mechanism)
    archives = {}
    for name in PHASE_NAMES:
        path = args.parameters/(name+".coverage.json")
        if not path.is_file():
            parser.error(f"missing coverage parameter archive: {path}")
        archive = json.loads(path.read_text())
        if archive.get("format") != "arrhenius-coverage-thermo-v1":
            parser.error(f"unsupported coverage archive format: {path}")
        if archive.get("phase_name") != name:
            parser.error(f"archive {path} holds phase {archive.get('phase_name')}, expected {name}")
        if archive.get("source_sha256") != mechanism_sha256:
            parser.error(f"archive {path} was not exported from the supplied mechanism")
        archives[name] = sha(path)
    exporter = Path(__file__).resolve().parent.parent/"mechanism/export_coverage_thermo.py"
    native = np.load(args.julia_result)
    meta = {k[:-5]: bytes(native[k]).decode() for k in native.files if k.endswith("_utf8")}
    if meta["phase_names"].splitlines() != PHASE_NAMES:
        parser.error("native phase ordering differs")
    for name in PHASE_NAMES:
        if meta["archive_"+name+"_sha256"] != archives[name]:
            parser.error(f"native parameter archive changed after measurement: {name}.coverage.json")
    if meta["mechanism_source_sha256"] != mechanism_sha256:
        parser.error("native archives were exported from a different mechanism")
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
    phases = [ct.Interface(str(args.mechanism), name) for name in PHASE_NAMES]
    for phase in phases:
        phase.TP = 300., ct.one_atm
    prepare_seconds = time.perf_counter()-started
    started = time.perf_counter()
    output = calculation(phases)
    first_seconds = time.perf_counter()-started
    samples, batch_size = [], 0
    if args.qualification != "validate-only":
        started = time.perf_counter()
        warmup = calculation(phases)
        warmup_seconds = time.perf_counter()-started
        if not identical(output, warmup):
            raise RuntimeError("Cantera warmup changed output")
        batch_size = max(1, min(1000, math.ceil(.020/max(warmup_seconds, 1e-9))))
        gc.collect()
        for _ in range(args.repetitions):
            elapsed = 0.
            for _ in range(batch_size):
                started = time.perf_counter()
                repeated = calculation(phases)
                elapsed += time.perf_counter()-started
                if not identical(output, repeated):
                    raise RuntimeError("Cantera repetition changed output")
            samples.append(elapsed/batch_size)
            verify_numerical_threads()
    threads_after = verify_numerical_threads()
    reference_path = args.output.with_name(args.output.stem+".reference.npz")
    np.savez(reference_path, **output)
    checks = {}
    for key, expected_shape in SHAPES.items():
        candidate, reference = native[key], output[key]
        difference = np.abs(candidate-reference) if candidate.shape == reference.shape else None
        worst = np.unravel_index(np.argmax(difference), difference.shape) if difference is not None else None
        checks[key] = {
            "pass": bool(candidate.shape == reference.shape == expected_shape
                         and np.all(np.isfinite(candidate)) and np.all(np.isfinite(reference))
                         and np.allclose(candidate, reference, rtol=1e-11, atol=1e-12)),
            "native_shape": list(candidate.shape), "reference_shape": list(reference.shape),
            "expected_shape": list(expected_shape),
            "relative_tolerance": 1e-11, "absolute_tolerance": 1e-12,
            "maximum_absolute_error": float(difference[worst]) if worst is not None else None,
            "worst_index": [int(i) for i in worst] if worst is not None else None,
            "native_at_worst": float(candidate[worst]) if worst is not None else None,
            "cantera_at_worst": float(reference[worst]) if worst is not None else None}
    ii, jj = np.indices(SHAPES["cross"])
    zero_mask = ii+jj > 100
    checks["cross_upper_triangle"] = {
        "pass": bool(int(np.count_nonzero(zero_mask)) == 5050
                     and np.all(native["cross"][zero_mask] == 0.) and np.all(output["cross"][zero_mask] == 0.)),
        "unevaluated_states": int(np.count_nonzero(zero_mask)),
        "evaluated_states": STATES_PER_SWEEP,
        "scope": "exactly zero entries above the occupied-site boundary in both implementations; all 5151 triangular states plus the 404 curve states are evaluated, endpoints included"}
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
    report = {"example": "thermo/coverage_dependent_surf", "cantera_version": ct.__version__,
              "cantera_source_sha": COMMIT, "published_source_sha256": SOURCE_SHA256,
              "mechanism_sha256": mechanism_sha256, "parameter_archives_sha256": archives,
              "parameter_exporter_sha256": sha(exporter),
              "cantera_build_record_sha256": sha(args.cantera_build_record),
              "loaded_libraries_match_build_record": build_matches,
              "cantera_shared_libraries_sha256": libraries, "cantera_extension_sha256": sha(compiled.__file__),
              "harness_sha256": sha(__file__), "environment_helper_sha256": sha(Path(__file__).with_name("benchmark_environment.py")),
              "native_metadata": meta, "native_artifact_sha256": sha(args.julia_result),
              "reference_artifact": reference_path.name, "reference_artifact_sha256": sha(reference_path),
              "hardware": hardware, "benchmark_target": args.target,
              "numerical_threads": threads_before, "numerical_threads_after": threads_after,
              "native_numerical_threads": native_threads,
              "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "scope": "five prepared interface phases at 300 K/1 atm; four 101-state h/RT and s/R curves plus the 5151-state triangular cross-interaction h/RT map per complete call; phase construction, imports, plots and checks excluded",
              "checks": checks, "correctness_pass": correct, "phase_names": PHASE_NAMES,
              "cantera_import_seconds": IMPORT_SECONDS, "cantera_prepare_seconds": prepare_seconds,
              "julia_import_seconds": float(native["import_seconds"][0]),
              "julia_prepare_seconds": float(native["prepare_seconds"][0]),
              "cantera_first_seconds": first_seconds, "julia_first_seconds": float(native["first_seconds"][0]),
              "cantera_warm_seconds": samples, "julia_warm_seconds": native_samples,
              "sample_definition": "mean timed seconds per complete 5555-state sweep within each warm batch; checks excluded",
              "batch_policy": "ceil(0.020 / warmup seconds), clamped to 1..1000 complete calculations",
              "cantera_batch_size": batch_size, "julia_batch_size": native_batch,
              "warm_outputs_checked": len(samples)*batch_size, "native_warm_outputs_checked": native_checked,
              "speed_ratio": ratio, "minimum_speed_ratio": .95,
              "qualification": "controlled" if controlled else "not_qualified",
              "performance_pass": bool(controlled and correct and ratio >= .95)}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
    print("coverage sweep", "correctness", correct, "speed ratio", ratio, "controlled", controlled, flush=True)
    if not correct:
        raise SystemExit("native coverage sweep failed reference checks")


if __name__ == "__main__":
    main()
