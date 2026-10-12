"""Prepare the complete pinned transport/multiprocessing_viscosity example and check native Julia output.

Run with --source-example, --mechanism and --output-dir to stage inputs and
compute the full source reference. After parallel_transport_cases.jl, add
--julia-result to validate all native arrays, thread settings and provenance
hashes. This script checks calculations; it does not qualify performance.
"""
import argparse
import hashlib
import importlib
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path

import cantera as ct
import numpy as np

SOURCE_SHA256 = "6dd426c24a6af5a5b9d2e87d2d9b33e78de83645952b946e81d34e727debc1de"
MECHANISM_SHA256 = "06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345"
SOURCE_URL = ("https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3"
              "/samples/python/transport/multiprocessing_viscosity.py")
ROOT = Path(__file__).resolve().parents[1]
N_PROCS, N_TEMPS = 4, 5000
RTOL, ATOL = 1e-10, 1e-15
PROPERTIES = (("conductivity", "get_thermal_conductivity"), ("viscosity", "get_viscosity"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def loaded_cantera_hashes(environment):
    """Hashes of the mapped Cantera libraries and loaded extension; filenames only."""
    hashes = {p.name: sha(p) for p in environment.loaded_library_paths()
              if "cantera" in p.name.lower() and p.is_file()}
    for module in list(sys.modules.values()):
        file = getattr(module, "__file__", None)
        if file and "cantera" in Path(file).name.lower() \
                and Path(file).suffix in (".so", ".pyd", ".dylib", ".dll") and Path(file).is_file():
            hashes.setdefault(Path(file).name, sha(file))
    return dict(sorted(hashes.items()))


def core_tree():
    base = ROOT / "src"
    paths = sorted(p.relative_to(base).as_posix() for p in base.rglob("*") if p.is_file())
    return paths, [sha(base / p) for p in paths]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-example", type=Path, required=True)
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--julia-result", type=Path)
    args = parser.parse_args()
    if sha(args.source_example) != SOURCE_SHA256:
        parser.error("source example differs from the pinned Cantera calculation")
    if sha(args.mechanism) != MECHANISM_SHA256:
        parser.error("mechanism differs from the pinned Cantera 4 gri30.yaml")
    if not ct.__version__.startswith("4.0"):
        parser.error("Cantera 4.0 development reference required")
    environment = load_module("benchmark_environment", ROOT / "validation" / "benchmark_environment.py")
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    mechanism = out / "gri30.yaml"
    sidecar = Path(str(mechanism) + ".npz")
    transport = Path(str(mechanism) + ".multicomponent.npz")
    if args.julia_result:
        missing = [p.name for p in (mechanism, sidecar, transport) if not p.is_file()]
        if missing:
            parser.error(f"validation-only rerun requires staged inputs: {missing}")
    else:
        if not (mechanism.is_file() and os.path.samefile(args.mechanism, mechanism)):
            shutil.copyfile(args.mechanism, mechanism)
        load_module("export_sidecar", ROOT / "mechanism" / "export_sidecar.py").export(mechanism, sidecar)
        load_module("export_multicomponent", ROOT / "mechanism" / "export_multicomponent.py"
                    ).export_multicomponent(str(mechanism), transport)
    if sha(mechanism) != MECHANISM_SHA256:
        parser.error("staged mechanism differs from the pinned Cantera 4 gri30.yaml")

    sys.path.insert(0, str(args.source_example.resolve().parent))
    original = importlib.import_module(args.source_example.stem)
    if Path(original.__file__).resolve() != args.source_example.resolve():
        parser.error("original module was not imported from the verified source file")
    references = {"T": np.linspace(300.0, 900.0, N_TEMPS)}
    for name, callback in PROPERTIES:
        predicate = getattr(original, callback)
        references[name + "_parallel"] = np.array(original.parallel(str(mechanism), predicate, N_PROCS, N_TEMPS))
        references[name + "_serial"] = np.array(original.serial(str(mechanism), predicate, N_TEMPS))
    np.savez(out / "fullreference.npz", **references)

    checks = {}
    for name, _ in PROPERTIES:
        checks[f"reference_{name}_replay"] = {
            "pass": bool(np.array_equal(references[name + "_parallel"], references[name + "_serial"])),
            "comparison": "exact parallel/serial replay"}
    report = dict(
        example="transport/multiprocessing_viscosity",
        cantera_version=ct.__version__,
        source_sha256=SOURCE_SHA256,
        mechanism_sha256=MECHANISM_SHA256,
        sidecar_sha256=sha(sidecar),
        multicomponent_sha256=sha(transport),
        fullreference_sha256=sha(out / "fullreference.npz"),
        cantera_loaded_sha256=loaded_cantera_hashes(environment),
        host=environment.host_metadata(),
        source_processes=N_PROCS,
        temperatures=N_TEMPS,
        pressure_Pa=101325.0,
        composition="CH4:1.0, O2:1.0, N2:3.76",
        transport_model="multicomponent",
        rtol=RTOL,
        atol=ATOL,
        performance_qualified=False,
        source_url=SOURCE_URL,
        checks=checks,
        native_status="unvalidated")
    passed = all(row["pass"] for row in checks.values())

    if args.julia_result:
        native = np.load(args.julia_result)
        metadata = {k[:-5]: bytes(native[k]).decode() for k in native.files if k.endswith("_utf8")}
        required = sorted(set(references) | {"native_checks_pass", "julia_threads", "blas_threads"})
        missing = [key for key in required if key not in native.files]
        native_checks = {"schema": {"pass": not missing, "missing": missing}}
        for key, expected in references.items():
            if key not in native.files:
                continue
            actual = native[key]
            same_shape = actual.shape == expected.shape
            native_checks[key] = {
                "pass": bool(same_shape and np.isfinite(actual).all()
                             and np.allclose(actual, expected, rtol=RTOL, atol=ATOL)),
                "rtol": RTOL, "atol": ATOL,
                "maximum_absolute_error": float(np.max(np.abs(actual - expected))) if same_shape else None}
        for name, _ in PROPERTIES:
            if name + "_parallel" in native.files and name + "_serial" in native.files:
                native_checks[f"native_{name}_replay"] = {
                    "pass": bool(np.array_equal(native[name + "_parallel"], native[name + "_serial"])),
                    "comparison": "exact parallel/serial replay"}
        if "native_checks_pass" in native.files:
            flag = native["native_checks_pass"]
            native_checks["native_checks_pass"] = {"pass": bool(flag.shape == (1,) and bool(flag[0]))}
        if "julia_threads" in native.files:
            native_checks["julia_threads"] = {"pass": bool(int(native["julia_threads"][0]) == N_PROCS),
                                              "value": int(native["julia_threads"][0])}
        if "blas_threads" in native.files:
            native_checks["blas_threads"] = {"pass": bool(int(native["blas_threads"][0]) == 1),
                                             "value": int(native["blas_threads"][0])}
        provenance = {}
        for key, path in (("example_sha256", ROOT / "example" / "transport" / "multiprocessing_viscosity.jl"),
                          ("harness_sha256", ROOT / "validation" / "parallel_transport_cases.jl"),
                          ("mechanism_sha256", mechanism),
                          ("sidecar_sha256", sidecar),
                          ("transport_sha256", transport)):
            actual = sha(path) if path.is_file() else None
            provenance[key] = {"pass": bool(actual and metadata.get(key) == actual),
                               "file": path.name, "sha256": actual}
        paths, hashes = core_tree()
        provenance["core_source_tree"] = {
            "pass": metadata.get("core_source_paths") == "\n".join(paths)
                    and metadata.get("core_source_sha256") == "\n".join(hashes),
            "files": len(paths)}
        native_pass = (all(row["pass"] for row in native_checks.values())
                       and all(row["pass"] for row in provenance.values()))
        report.update(
            native_status="validated" if native_pass else "failed",
            native_artifact_sha256=sha(args.julia_result),
            native_metadata={key: metadata.get(key)
                             for key in ("julia_version", "cpu", "system", "kernel_release")},
            native_checks=native_checks,
            native_provenance=provenance)
        passed = passed and native_pass

    report["passed"] = bool(passed)
    (out / "checks.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf8")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit("transport reference replay or native comparison failed")


if __name__ == "__main__":
    main()
