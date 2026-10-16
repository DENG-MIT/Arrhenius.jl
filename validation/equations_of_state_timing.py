"""Time the complete pinned CO2 equations-of-state example against native Julia output.

The timed call is the original pinned source calculation: ``get_thermo_Cantera``
on the CO2-Ideal and CO2-RK phases and ``get_thermo_CoolProp``, over all 1000
original pressures at 300 K, returning the original (h, u, s, cp, cv) tuples per
equation of state. Imports, model loading, pristine-state restore, plotting and
every comparison are excluded from the timers. The pristine phase state is
restored before every complete call, and every warm output is checked bitwise
against the first call. Density, phase and saturation references are computed
independently outside the measured source call.
"""
from pathlib import Path
import argparse
import ast
import gc
import hashlib
import json
import math
import statistics
import time

from benchmark_environment import (host_metadata, matches_target, cantera_library_hashes,
                                   loaded_library_paths, verify_numerical_threads)

started = time.perf_counter()
from equations_of_state_cases import SOURCE_SHA256
import numpy as np
import cantera as ct
import cantera._cantera as compiled
import CoolProp
import CoolProp.CoolProp as coolprop_compiled
from CoolProp.CoolProp import PropsSI, PhaseSI, get_fluid_param_string
IMPORT_SECONDS = time.perf_counter()-started

COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
EXPECTED_CLAPEYRON = "0.6.28"
ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = Path("example/thermodynamics/equations_of_state.jl")
MODELS = ("ideal", "rk", "helmholtz")
PHASES = {"ideal": "CO2-Ideal", "rk": "CO2-RK"}
PROPERTIES = ["relative_enthalpy", "relative_internal_energy", "relative_entropy", "cp", "cv"]
PROPERTY_UNITS = ["kJ/kg", "kJ/kg", "kJ/kg/K", "kJ/kg/K", "kJ/kg/K"]
PROPERTY_RTOL = PROPERTY_ATOL = 1e-7
STATE_RTOL = STATE_ATOL = 1e-8
IDENTITY_RTOL = IDENTITY_ATOL = 1e-9
COOLPROP_PHASES = {"gas", "liquid", "supercritical_liquid"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_functions(path):
    selected = [n for n in ast.parse(Path(path).read_text()).body
                if isinstance(n, ast.FunctionDef)
                and n.name in ("get_thermo_Cantera", "get_thermo_CoolProp")]
    if len(selected) != 2:
        raise RuntimeError("pinned source must define exactly the two computational functions")
    namespace = dict(ct=ct, np=np, PropsSI=PropsSI)
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), "exec"), namespace)
    return namespace["get_thermo_Cantera"], namespace["get_thermo_CoolProp"]


def capture_state(phase):
    return np.array(phase.state, copy=True)


def restore_state(phase, pristine):
    phase.state = pristine
    if not np.array_equal(phase.state, pristine):
        raise RuntimeError("source phase state restore failed")


def calculation(get_thermo_cantera, get_thermo_coolprop, ideal, rk, T, p):
    ideal_output = get_thermo_cantera(ideal, T, p)
    rk_output = get_thermo_cantera(rk, T, p)
    helmholtz_output = get_thermo_coolprop(T, p)
    return ideal_output, rk_output, helmholtz_output


def stacked(outputs):
    return tuple(np.asarray(row) for model in outputs for row in model)


def check_bitwise(reference, candidate):
    for first, repeated in zip(reference, stacked(candidate), strict=True):
        if first.shape != repeated.shape or not np.array_equal(first, repeated):
            raise RuntimeError("complete source calculation changed output between calls")


def independent_references(mechanism, T, p):
    references = {}
    for key, phase in PHASES.items():
        states = ct.SolutionArray(ct.Solution(str(mechanism), phase), len(p))
        states.TPX = T, p, "CO2:1"
        references[key + "_density"] = np.array(states.density)
    references["helmholtz_density"] = np.array([PropsSI("D", "P", pi, "T", T, "CO2") for pi in p])
    labels = [PhaseSI("P", pi, "T", T, "CO2") for pi in p]
    if not set(labels) <= COOLPROP_PHASES:
        raise RuntimeError(f"unexpected CoolProp phase labels at 300 K: {sorted(set(labels))}")
    # CoolProp labels the compressed liquid above Pc "supercritical_liquid" even
    # though this 300 K sweep remains below Tc; both liquid labels map to 2.
    references["helmholtz_phase"] = np.array([1 if label == "gas" else 2 for label in labels])
    references["coolprop_phase_labels"] = sorted(set(labels))
    references["saturation_pressure"] = np.array([PropsSI("P", "T", T, "Q", 0, "CO2")])
    return references


def property_check(actual, expected):
    same_shape = actual.shape == expected.shape
    difference = np.abs(actual - expected) if same_shape else None
    passed = (same_shape and bool(np.all(np.isfinite(actual)))
              and bool(np.allclose(actual, expected, rtol=PROPERTY_RTOL, atol=PROPERTY_ATOL)))
    return {"pass": bool(passed), "relative_tolerance": PROPERTY_RTOL,
            "absolute_tolerance": PROPERTY_ATOL, "property_names": PROPERTIES,
            "property_units": PROPERTY_UNITS,
            "maximum_absolute_error": float(difference.max()) if same_shape else None,
            "maximum_absolute_error_by_property": difference.max(axis=1).tolist() if same_shape else None}


def state_check(actual, expected, rtol=STATE_RTOL, atol=STATE_ATOL):
    same_shape = actual.shape == expected.shape
    passed = (same_shape and bool(np.all(np.isfinite(actual)))
              and bool(np.allclose(actual, expected, rtol=rtol, atol=atol)))
    return {"pass": bool(passed), "relative_tolerance": rtol, "absolute_tolerance": atol,
            "maximum_absolute_error": float(np.max(np.abs(actual - expected))) if same_shape else None}


def native_diagnostics(native, p):
    pressures = {"ideal": p, "rk": np.asarray(native["rk_pressure"], dtype=float),
                 "helmholtz": np.asarray(native["helmholtz_pressure"], dtype=float)}
    diagnostics = {}
    for model in MODELS:
        matrix = np.asarray(native[model], dtype=float)
        density = np.asarray(native[model + "_density"], dtype=float)
        pv = pressures[model]/density/1000
        diagnostics[model] = {
            "reference_point_exact_zero": bool(np.all(matrix[:3, 0] == 0.0)),
            "enthalpy_internal_energy_identity": {
                "pass": bool(np.allclose(matrix[0] - matrix[1], pv - pv[0],
                                         rtol=IDENTITY_RTOL, atol=IDENTITY_ATOL)),
                "relative_tolerance": IDENTITY_RTOL, "absolute_tolerance": IDENTITY_ATOL,
                "maximum_absolute_error": float(np.max(np.abs((matrix[0] - matrix[1]) - (pv - pv[0])))),
                "definition": "h - u equals p/density, first pressure point referenced, kJ/kg"}}
    native_pressure = np.asarray(native["pressure"], dtype=float)
    native_saturation = float(native["saturation_pressure"][0])
    derived_phase = np.where(native_pressure < native_saturation, 1, 2)
    diagnostics["phase_saturation_consistency"] = {
        "pass": bool(np.array_equal(np.asarray(native["helmholtz_phase"]), derived_phase)),
        "definition": "native phase flags reproduce the native p < psat(300 K) split, 1 gas / 2 liquid"}
    return diagnostics


def verify_core_sources(meta, root):
    src = root/"src"
    recorded_paths = meta["core_source_paths"].splitlines()
    recorded_hashes = meta["core_source_sha256"].splitlines()
    if len(recorded_paths) != len(recorded_hashes) or recorded_paths != sorted(recorded_paths):
        raise RuntimeError("native core source record is malformed")
    actual = sorted(path.relative_to(src).as_posix() for path in src.rglob("*") if path.is_file())
    if recorded_paths != actual:
        missing = sorted(set(recorded_paths) - set(actual))
        extra = sorted(set(actual) - set(recorded_paths))
        raise RuntimeError(f"native core source tree differs: missing={missing[:5]} extra={extra[:5]}")
    for relative, digest in zip(recorded_paths, recorded_hashes, strict=True):
        if sha(src/relative) != digest:
            raise RuntimeError(f"native source changed after measurement: src/{relative}")
    return len(recorded_paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-example", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--julia-result", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cantera-build-record", type=Path, required=True)
    parser.add_argument("--native-package-root", type=Path)
    parser.add_argument("--target", choices=("wsl", "apple-m4"), required=True)
    parser.add_argument("--qualification", choices=("controlled", "informational", "validate-only"),
                        default="informational")
    parser.add_argument("--repetitions", type=int, default=9)
    args = parser.parse_args()
    if args.repetitions < 9:
        parser.error("at least nine warm batches required")
    if sha(args.source_example) != SOURCE_SHA256:
        parser.error("source example differs from the pinned calculation")
    mechanism = args.input_dir/"co2-thermo.yaml"
    sidecar = args.input_dir/"co2-thermo.yaml.npz"
    parameters = args.input_dir/"carbon-dioxide.json"
    for path in (mechanism, sidecar, parameters):
        if not path.is_file():
            parser.error(f"missing input file: {path.name}")
    root = args.native_package_root or ROOT
    native = np.load(args.julia_result)
    meta = {k[:-5]: bytes(native[k]).decode() for k in native.files if k.endswith("_utf8")}
    required_meta = ("julia_version", "clapeyron_version", "qualification", "cpu", "system",
                     "kernel_release", "timestamp_utc", "harness_sha256", "example_sha256",
                     "parameter_sha256", "mechanism_sha256", "sidecar_sha256", "thread_helper_sha256", "core_source_paths",
                     "core_source_sha256", "first_call_scope", "numerical_thread_names")
    missing_meta = [key for key in required_meta if key not in meta]
    if missing_meta:
        raise RuntimeError(f"native metadata incomplete: {missing_meta}")
    if meta["parameter_sha256"] != sha(parameters):
        raise RuntimeError("native Helmholtz parameters differ from the staged file")
    if meta["mechanism_sha256"] != sha(mechanism):
        raise RuntimeError("native mechanism differs from the staged file")
    if meta["sidecar_sha256"] != sha(sidecar):
        raise RuntimeError("native mechanism sidecar differs from the staged file")
    if meta["thread_helper_sha256"] != sha(Path(__file__).with_name("numerical_threads.jl")):
        raise RuntimeError("native numerical thread helper changed after measurement")
    if meta["example_sha256"] != sha(root/EXAMPLE):
        raise RuntimeError("native example source changed after measurement")
    core_source_count = verify_core_sources(meta, root)
    harness = Path(__file__).with_suffix(".jl")
    harness_verified = harness.is_file() and sha(harness) == meta["harness_sha256"]
    if not harness_verified:
        raise RuntimeError("native timing driver differs from the measured harness")
    parsed = json.loads(parameters.read_text(encoding="utf8"))
    exported = json.loads(get_fluid_param_string("CO2", "JSON"))
    parameters_match = len(exported) == 1 and exported[0] == parsed
    if not parameters_match:
        raise RuntimeError("staged Helmholtz parameters differ from the current CoolProp full export")
    build = json.loads(args.cantera_build_record.read_text())
    libraries = cantera_library_hashes(ct.__file__)
    recorded = {Path(k).name: v for k, v in build["library_hashes"].items()}
    build_matches = bool(libraries and all(recorded.get(k) == v for k, v in libraries.items())
                         and recorded.get(Path(compiled.__file__).name) == sha(compiled.__file__))
    if not build_matches or build["source"]["commit"] != COMMIT or not ct.__version__.startswith("4.0"):
        parser.error("loaded Cantera libraries differ from the pinned development build")
    coolprop_libraries = {p.name: sha(p) for p in loaded_library_paths()
                          if "coolprop" in p.name.lower() and p.is_file()}
    cantera_extension_sha = sha(compiled.__file__)
    coolprop_extension_sha = sha(coolprop_compiled.__file__)
    helper = Path(__file__).with_name("benchmark_environment.py")
    input_paths = (args.source_example, mechanism, sidecar, parameters, args.cantera_build_record,
                   Path(__file__), harness, helper, root/EXAMPLE)
    input_hashes = tuple(sha(path) for path in input_paths)
    threads_before = verify_numerical_threads(set_accelerate=True)
    hardware = host_metadata()
    get_thermo_cantera, get_thermo_coolprop = source_functions(args.source_example)
    started = time.perf_counter()
    ideal = ct.Solution(str(mechanism), PHASES["ideal"])
    rk = ct.Solution(str(mechanism), PHASES["rk"])
    pristine = {model: capture_state(phase) for model, phase in (("ideal", ideal), ("rk", rk))}
    preparation_seconds = time.perf_counter()-started
    T, p = 300.0, 1e5*np.linspace(1, 100, 1000)
    restore_state(ideal, pristine["ideal"])
    restore_state(rk, pristine["rk"])
    started = time.perf_counter()
    first = calculation(get_thermo_cantera, get_thermo_coolprop, ideal, rk, T, p)
    first_seconds = time.perf_counter()-started
    reference = stacked(first)
    samples, batch_size = [], 0
    if args.qualification != "validate-only":
        restore_state(ideal, pristine["ideal"])
        restore_state(rk, pristine["rk"])
        started = time.perf_counter()
        warmup = calculation(get_thermo_cantera, get_thermo_coolprop, ideal, rk, T, p)
        warmup_seconds = time.perf_counter()-started
        check_bitwise(reference, warmup)
        batch_size = max(1, min(1000, math.ceil(.020/max(warmup_seconds, 1e-9))))
        gc.collect()
        for _ in range(args.repetitions):
            elapsed = 0.
            for _ in range(batch_size):
                restore_state(ideal, pristine["ideal"])
                restore_state(rk, pristine["rk"])
                started = time.perf_counter()
                repeated = calculation(get_thermo_cantera, get_thermo_coolprop, ideal, rk, T, p)
                elapsed += time.perf_counter()-started
                check_bitwise(reference, repeated)
            samples.append(elapsed/batch_size)
            verify_numerical_threads()
    threads_after = verify_numerical_threads()
    references = independent_references(mechanism, T, p)
    if threads_before != threads_after:
        raise RuntimeError("numerical thread settings changed during measurement")
    if (cantera_library_hashes(ct.__file__) != libraries
            or sha(compiled.__file__) != cantera_extension_sha
            or sha(coolprop_compiled.__file__) != coolprop_extension_sha
            or {p.name: sha(p) for p in loaded_library_paths()
                if "coolprop" in p.name.lower() and p.is_file()} != coolprop_libraries):
        raise RuntimeError("loaded calculation libraries changed during measurement")
    if tuple(sha(path) for path in input_paths) != input_hashes:
        raise RuntimeError("benchmark source or input files changed during measurement")
    verify_core_sources(meta, root)
    expected = {"T": np.array([T]), "pressure": p, "ideal": np.array(first[0]),
                "rk": np.array(first[1]), "helmholtz": np.array(first[2]),
                "ideal_density": references["ideal_density"], "rk_density": references["rk_density"],
                "helmholtz_density": references["helmholtz_density"],
                "helmholtz_phase": references["helmholtz_phase"],
                "saturation_pressure": references["saturation_pressure"]}
    required_arrays = ("T", "pressure", "ideal", "rk", "helmholtz", "ideal_density", "rk_density",
                       "helmholtz_density", "rk_pressure", "helmholtz_pressure", "helmholtz_phase",
                       "saturation_pressure", "seconds", "batch_size", "warm_outputs_checked",
                       "import_seconds", "preparation_seconds", "first_invocation_seconds",
                       "numerical_threads")
    missing_arrays = [key for key in required_arrays if key not in native.files]
    if missing_arrays:
        raise RuntimeError(f"native artifact misses arrays: {missing_arrays}")
    checks = {model: property_check(np.asarray(native[model], dtype=float), expected[model])
              for model in MODELS}
    for key in ("T", "pressure", "ideal_density", "rk_density", "helmholtz_density",
                "saturation_pressure"):
        checks[key] = state_check(np.asarray(native[key], dtype=float), expected[key])
    checks["rk_pressure"] = state_check(np.asarray(native["rk_pressure"], dtype=float), p)
    checks["helmholtz_pressure"] = state_check(np.asarray(native["helmholtz_pressure"], dtype=float), p)
    native_phase = np.asarray(native["helmholtz_phase"])
    phase_shape = native_phase.shape == expected["helmholtz_phase"].shape
    checks["helmholtz_phase"] = {
        "pass": bool(phase_shape and np.array_equal(native_phase, expected["helmholtz_phase"])),
        "relative_tolerance": 0, "absolute_tolerance": 0, "phase_encoding": {"gas": 1, "liquid": 2},
        "coolprop_phase_labels": references["coolprop_phase_labels"],
        "maximum_absolute_error": (float(np.max(np.abs(native_phase - expected["helmholtz_phase"])))
                                   if phase_shape else None)}
    flags = {key: (bool(native[key][0]) if key in native.files else None)
             for key in ("native_checks_pass", "strict_replay_checked")}
    diagnostics = native_diagnostics(native, p)
    native_samples = np.asarray(native["seconds"], dtype=float).tolist()
    native_batch = int(native["batch_size"][0])
    native_checked = int(native["warm_outputs_checked"][0])
    if native_checked != len(native_samples)*native_batch:
        raise RuntimeError("inconsistent native replay count")
    if any(not math.isfinite(v) or v <= 0 for v in samples + native_samples):
        raise RuntimeError("invalid elapsed time")
    thread_names = meta["numerical_thread_names"].splitlines()
    if len(thread_names) != len(native["numerical_threads"]):
        raise RuntimeError("native numerical thread record is malformed")
    native_threads = dict(zip(thread_names, np.asarray(native["numerical_threads"]).tolist()))
    threads_match = (native_threads.get("julia_threads") == native_threads.get("blas_threads") == 1
                     and all(v == 1 for v in native_threads.values())
                     and (args.target != "apple-m4" or native_threads.get("accelerate_threading_mode") == 1))
    ratio = statistics.median(samples)/statistics.median(native_samples) if samples and native_samples else None
    diagnostics_pass = all(d["reference_point_exact_zero"] and d["enthalpy_internal_energy_identity"]["pass"]
                           for m, d in diagnostics.items() if m in MODELS)
    diagnostics_pass = diagnostics_pass and diagnostics["phase_saturation_consistency"]["pass"]
    flags_pass = all(value is True for value in flags.values())
    correct = all(check["pass"] for check in checks.values()) and diagnostics_pass and flags_pass
    controlled = bool(args.qualification == meta["qualification"] == "controlled"
                      and matches_target(hardware, args.target) and matches_target(meta, args.target)
                      and all(hardware[k] == meta[k] for k in ("cpu", "kernel_release"))
                      and threads_match and len(samples) >= 9 and len(native_samples) >= 9
                      and 1 <= batch_size <= 1000 and 1 <= native_batch <= 1000
                      and meta["clapeyron_version"] == EXPECTED_CLAPEYRON
                      and harness_verified and correct)
    report = {"example": "thermo/equations_of_state", "cantera_version": ct.__version__,
              "cantera_source_sha": COMMIT, "coolprop_version": CoolProp.__version__,
              "coolprop_extension_sha256": coolprop_extension_sha,
              "coolprop_loaded_libraries_sha256": coolprop_libraries,
              "published_source_sha256": SOURCE_SHA256,
              "source_url": "https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/thermo/equations_of_state.py",
              "cantera_build_record_sha256": sha(args.cantera_build_record),
              "loaded_libraries_match_build_record": build_matches,
              "cantera_shared_libraries_sha256": libraries,
              "cantera_extension_sha256": cantera_extension_sha,
              "libraries_unchanged_after_measurement": True,
              "harness_sha256": sha(__file__), "paired_julia_harness": harness.name,
              "native_harness_sha256": meta["harness_sha256"], "native_harness_verified": harness_verified,
              "environment_helper_sha256": sha(helper),
              "native_metadata": meta, "native_artifact_sha256": sha(args.julia_result),
              "native_flags": flags, "core_source_files_verified": core_source_count,
              "julia_version": meta["julia_version"], "clapeyron_version": meta["clapeyron_version"],
              "expected_clapeyron_version": EXPECTED_CLAPEYRON,
              "hardware": hardware, "benchmark_target": args.target,
              "numerical_threads": threads_before, "numerical_threads_after": threads_after,
              "native_numerical_threads": native_threads,
              "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "scope": "complete pinned source sweep: original get_thermo_Cantera on CO2-Ideal and CO2-RK plus original get_thermo_CoolProp over all 1000 pressures at 300 K, returning the unmodified (h, u, s, cp, cv) tuples per equation of state; imports, model loading, state restore, plotting and all comparisons excluded",
              "first_call_scope": {"cantera": "cold complete three-EOS source calculation after model preparation; per-call SolutionArray construction and CoolProp initialization inside the original functions included; imports, model loading and pristine state restore excluded",
                                   "julia": meta["first_call_scope"]},
              "equations_of_state": list(MODELS), "property_names": PROPERTIES,
              "property_units": PROPERTY_UNITS, "temperature_kelvin": T, "pressure_points": len(p),
              "input_files": {"mechanism": mechanism.name, "sidecar": sidecar.name, "parameters": parameters.name},
              "mechanism_sha256": sha(mechanism), "sidecar_sha256": sha(sidecar), "parameter_json_sha256": sha(parameters),
              "parameters_match_coolprop_export": parameters_match,
              "checks": checks, "native_diagnostics": diagnostics, "correctness_pass": correct,
              "cantera_import_seconds": IMPORT_SECONDS,
              "cantera_model_preparation_seconds": preparation_seconds,
              "cantera_first_invocation_seconds": first_seconds,
              "julia_import_seconds": float(native["import_seconds"][0]),
              "julia_preparation_seconds": float(native["preparation_seconds"][0]),
              "julia_first_invocation_seconds": float(native["first_invocation_seconds"][0]),
              "cantera_warm_seconds": samples, "julia_warm_seconds": native_samples,
              "sample_definition": "mean timed seconds per complete three-EOS source calculation within each warm batch; state restore and checks excluded",
              "batch_policy": "ceil(0.020 / warmup seconds), clamped to 1..1000 complete calculations",
              "cantera_batch_size": batch_size, "julia_batch_size": native_batch,
              "warm_outputs_checked": len(samples)*batch_size,
              "native_warm_outputs_checked": native_checked,
              "speed_ratio": ratio, "minimum_speed_ratio": .95,
              "requested_qualification": args.qualification,
              "qualification": "controlled" if controlled else "not_qualified",
              "performance_pass": bool(controlled and ratio is not None and ratio >= .95)}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
    print("equations of state", "correctness", correct, "speed ratio", ratio,
          "controlled", controlled, flush=True)
    if not correct:
        raise SystemExit("native CO2 equations of state failed reference checks")


if __name__ == "__main__":
    main()
