"""Validate and time all 50 published gas/graphite HP equilibria with Cantera 4.

Run condensed_equilibrium_timing.jl first, on an otherwise idle target host.
Prepared input models are reused; every call creates and solves each mixture
from its stated inlet composition. Imports, checks and file I/O are not timed.
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import math
import statistics
import time
import tomllib
from benchmark_environment import (cantera_library_hashes, host_metadata,
                                   matches_target, verify_numerical_threads)

started = time.perf_counter()
import numpy as np
import cantera as ct
import cantera._cantera as compiled
IMPORT_SECONDS = time.perf_counter()-started
COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
SOURCE_SHA256 = "d110b027ffc26a3b8417175f8aca047bbde718a63309a8c20090573d87e9eaff"
GAS_SHA256 = "06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345"
GRAPHITE_SHA256 = "6074b6dd691785403dee1eeebf8e583b747681c34a148746bab2256282939cbd"
SHAPES = {"phi": (50,), "T": (50,), "species_moles": (54, 50)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculation(gas, solid, *, refined=False):
    phi = np.linspace(.3, 3.5, 50)
    temperature = np.empty(50)
    moles = np.empty((54, 50))
    for j, value in enumerate(phi):
        gas.set_equivalence_ratio(value, "CH4", "O2:1,N2:3.76")
        mixture = ct.Mixture([(gas, 1.), (solid, 0.)])
        mixture.T = 300.
        mixture.P = 101325.
        options = {"rtol": 1e-12, "max_steps": 3000, "max_iter": 300} if refined else {"max_steps": 1000}
        mixture.equilibrate("HP", solver="gibbs", **options)
        temperature[j] = mixture.T
        moles[:, j] = mixture.species_moles
    return {"phi": phi, "T": temperature, "species_moles": moles}


def identical(left, right):
    return left.keys() == right.keys() and all(np.array_equal(left[k], right[k]) for k in left)


def compare(candidate, reference):
    tolerances = {"phi": (0., 5e-16), "T": (0., 1e-3), "species_moles": (2e-5, 1e-9)}
    checks = {}
    for key, shape in SHAPES.items():
        left, right = candidate[key], reference[key]
        rtol, atol = tolerances[key]
        dimensions = left.shape == right.shape == shape
        difference = np.abs(left-right) if dimensions else None
        checks[key] = {"pass": bool(dimensions and np.all(np.isfinite(left)) and np.all(np.isfinite(right))
                                    and np.allclose(left, right, rtol=rtol, atol=atol)),
                       "maximum_absolute_error": float(np.max(difference)) if dimensions else None,
                       "relative_tolerance": rtol, "absolute_tolerance": atol,
                       "expected_shape": list(shape)}
    return checks


def physical_checks(candidate, gas, solid):
    """Recompute native balances and graphite stability with independent CT properties."""
    elements = gas.element_names
    atoms = np.array([[gas.n_atoms(k, element) for k in range(gas.n_species)] for element in elements])
    solid_atoms = np.array([solid.n_atoms(0, e) if e in solid.element_names else 0. for e in elements])
    carbon = gas.species_index("C")
    balances, enthalpies, affinities, products = [], [], [], []
    for j, value in enumerate(candidate["phi"]):
        gas.TP = 300., 101325.
        gas.set_equivalence_ratio(float(value), "CH4", "O2:1,N2:3.76")
        initial_elements = atoms @ gas.X
        initial_h = gas.enthalpy_mole
        hscale = max(abs(initial_h), 1e6*gas.mean_molecular_weight)
        moles = candidate["species_moles"][:, j]
        if np.any(moles < 0) or not np.all(np.isfinite(moles)) or np.sum(moles[:-1]) <= 0:
            return {"pass": False, "reason": "nonphysical species amounts"}
        gas.TPX = float(candidate["T"][j]), 101325., moles[:-1]
        solid.TP = float(candidate["T"][j]), 101325.
        final_elements = atoms @ moles[:-1] + solid_atoms*moles[-1]
        balances.append(float(np.max(np.abs(final_elements-initial_elements))/np.max(initial_elements)))
        final_h = np.sum(moles[:-1])*gas.enthalpy_mole + moles[-1]*solid.enthalpy_mole
        enthalpies.append(float(abs(final_h-initial_h)/hscale))
        affinity = float((gas.chemical_potentials[carbon]-solid.chemical_potentials[0])/(ct.gas_constant*gas.T))
        affinities.append(affinity)
        products.append(abs(moles[-1]*affinity))
    maxima = {"element_balance": max(balances), "enthalpy_balance": max(enthalpies),
              "graphite_affinity": max(affinities), "complementarity": max(products)}
    # Same predeclared balance limits as the full source/refined prototype check;
    # independent property evaluation permits 1e-8 in dimensionless chemical potential.
    return {"pass": bool(maxima["element_balance"] < 2e-9 and maxima["enthalpy_balance"] < 2e-9
                         and maxima["graphite_affinity"] < 1e-8 and maxima["complementarity"] < 1e-8),
            "maxima": maxima, "balance_limit": 2e-9, "chemical_potential_limit": 1e-8,
            "nonnegative_amounts": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("gas", "graphite", "parameters", "source-example", "julia-result", "output", "cantera-build-record"):
        parser.add_argument("--"+key, type=Path, required=True)
    parser.add_argument("--target", choices=("wsl", "apple-m4"), required=True)
    parser.add_argument("--qualification", choices=("controlled", "informational"), default="informational")
    parser.add_argument("--repetitions", type=int, default=9)
    args = parser.parse_args()
    if args.repetitions < 9:
        parser.error("at least nine complete warm sweeps required")
    root = Path(__file__).resolve().parent.parent
    if (sha(args.gas), sha(args.graphite), sha(args.source_example)) != (GAS_SHA256, GRAPHITE_SHA256, SOURCE_SHA256):
        parser.error("inputs or source example differ from the pinned complete calculation")
    native = np.load(args.julia_result)
    meta = tomllib.loads(Path(str(args.julia_result)+".toml").read_text())
    parameters = json.loads(args.parameters.read_text())
    if parameters.get("source_sha256") != GRAPHITE_SHA256:
        parser.error("condensed parameters were exported from a different mechanism")
    if meta["native_artifact_sha256"] != sha(args.julia_result):
        parser.error("native output was changed after measurement")
    current_inputs = {"gas": sha(args.gas), "sidecar": sha(Path(str(args.gas)+".npz")), "phase": sha(args.parameters)}
    def expected_files():
        files = {p.relative_to(root).as_posix() for p in (root/"src").rglob("*.jl")}
        files.update(("Project.toml", "example/thermodynamics/adiabatic.jl",
                      "validation/condensed_equilibrium_timing.jl", "validation/numerical_threads.jl"))
        return files
    if set(meta["source_hashes"]) != expected_files() or meta["input_hashes"] != current_inputs:
        parser.error("native source inventory or prepared inputs differ")
    reference_inputs = {path: sha(path) for path in (
        Path(__file__), Path(__file__).with_name("benchmark_environment.py"),
        root/"mechanism/export_condensed_phase.py", args.graphite, args.source_example,
        args.cantera_build_record, args.julia_result, Path(str(args.julia_result)+".toml"))}
    def verify_files():
        if set(meta["source_hashes"]) != expected_files():
            raise RuntimeError("native source inventory changed during measurement")
        for file, digest in meta["source_hashes"].items():
            if sha(root/file) != digest:
                raise RuntimeError(f"native source changed after measurement: {file}")
        if {"gas": sha(args.gas), "sidecar": sha(Path(str(args.gas)+".npz")), "phase": sha(args.parameters)} != current_inputs:
            raise RuntimeError("prepared inputs changed after measurement")
        if any(sha(path) != digest for path, digest in reference_inputs.items()):
            raise RuntimeError("reference input, driver or native artifact changed during measurement")
    verify_files()
    build = json.loads(args.cantera_build_record.read_text())
    recorded = {Path(k).name: v for k, v in build["library_hashes"].items()}
    libraries_before = cantera_library_hashes(ct.__file__)
    extension_before = sha(compiled.__file__)
    if not (libraries_before and all(recorded.get(k) == v for k, v in libraries_before.items())
            and recorded.get(Path(compiled.__file__).name) == extension_before
            and build["source"]["commit"] == COMMIT and ct.__version__.startswith("4.0")):
        parser.error("loaded Cantera differs from the pinned development build")
    threads_before = verify_numerical_threads(set_accelerate=True)
    hardware = host_metadata()
    started = time.perf_counter()
    gas, solid = ct.Solution(str(args.gas)), ct.Solution(str(args.graphite))
    preparation_seconds = time.perf_counter()-started
    started = time.perf_counter()
    output = calculation(gas, solid)
    first_seconds = time.perf_counter()-started
    # Each repetition starts from the same thermodynamic state; composition and
    # zero graphite inventory are reset inside each of its 50 solves.
    samples = []
    gc.collect()
    for _ in range(args.repetitions):
        gas.TP = solid.TP = 300., 101325.
        started = time.perf_counter()
        repeated = calculation(gas, solid)
        samples.append(time.perf_counter()-started)
        if not identical(output, repeated):
            raise RuntimeError("complete Cantera sweep changed on repetition")
        verify_numerical_threads()
    refined = calculation(gas, solid, refined=True)
    candidate = {k: native[k] for k in SHAPES}
    checks = {"source": compare(candidate, output), "refined": compare(candidate, refined),
              "physical": physical_checks(candidate, gas, solid)}
    verify_files()
    if libraries_before != cantera_library_hashes(ct.__file__) or extension_before != sha(compiled.__file__):
        raise RuntimeError("loaded Cantera library bytes changed during measurement")
    threads_after = verify_numerical_threads()
    if threads_before != threads_after:
        raise RuntimeError("Cantera thread settings changed")
    native_samples = native["seconds"].tolist()
    if (len(native_samples) < 9 or meta["warm_outputs_checked"] != len(native_samples)
            or not native_samples or any(not math.isfinite(v) or v <= 0 for v in native_samples+samples)):
        raise RuntimeError("invalid complete-workload elapsed samples or replay counts")
    correct = (all(c["pass"] for group in (checks["source"], checks["refined"]) for c in group.values())
               and checks["physical"]["pass"])
    native_threads = meta["numerical_threads"]
    controlled = bool(args.qualification == "controlled" and matches_target(hardware, args.target)
                      and matches_target(meta, args.target) and meta["source_unchanged"] and meta["inputs_unchanged"]
                      and all(hardware[k] == meta[k] for k in ("cpu", "kernel_release"))
                      and native_threads.get("julia_threads") == native_threads.get("blas_threads") == 1
                      and all(v == 1 for v in native_threads.values())
                      and (args.target != "apple-m4" or native_threads.get("accelerate_threading_mode") == 1))
    ratio = statistics.median(samples)/statistics.median(native_samples)
    reference_path = args.output.with_name(args.output.stem+".reference.npz")
    np.savez(reference_path, **{"source_"+k: v for k, v in output.items()}, **{"refined_"+k: v for k, v in refined.items()})
    report = {"example": "thermo/adiabatic", "cases": 50,
              "scope": "All 50 CH4/air phi values 0.3..3.5 at 300 K and 101325 Pa, initial gas 1 kmol and graphite 0 kmol; HP equilibrium temperature and all 54 species amounts. Prepared models; imports, checks, printing and file I/O excluded.",
              "cantera_version": ct.__version__, "cantera_source_sha": COMMIT,
              "published_source_sha256": SOURCE_SHA256, "gas_sha256": GAS_SHA256, "graphite_sha256": GRAPHITE_SHA256,
              "cantera_shared_libraries_sha256": libraries_before, "cantera_extension_sha256": extension_before,
              "cantera_build_record_sha256": sha(args.cantera_build_record), "loaded_libraries_match_build_record": True,
              "harness_sha256": sha(__file__), "environment_helper_sha256": sha(Path(__file__).with_name("benchmark_environment.py")),
              "exporter_sha256": sha(root/"mechanism/export_condensed_phase.py"), "native_metadata": meta,
              "native_artifact_sha256": sha(args.julia_result), "reference_artifact": reference_path.name,
              "reference_artifact_sha256": sha(reference_path), "hardware": hardware, "benchmark_target": args.target,
              "numerical_threads": threads_before, "numerical_threads_after": threads_after,
              "checks": checks, "correctness_pass": bool(correct),
              "cantera_import_seconds": IMPORT_SECONDS, "cantera_preparation_seconds": preparation_seconds,
              "cantera_first_seconds": first_seconds, "julia_first_seconds": float(native["first_seconds"][0]),
              "cantera_warm_seconds": samples, "julia_warm_seconds": native_samples,
              "warm_outputs_checked": args.repetitions, "native_warm_outputs_checked": meta["warm_outputs_checked"],
              "speed_ratio": ratio, "minimum_speed_ratio": .95,
              "qualification": "controlled" if controlled else "not_qualified",
              "performance_pass": bool(controlled and correct and ratio >= .95),
              "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print("graphite sweep", "correctness", correct, "speed ratio", ratio, "controlled", controlled, flush=True)
    if not correct:
        raise SystemExit("native graphite sweep failed independent reference checks")


if __name__ == "__main__":
    main()
