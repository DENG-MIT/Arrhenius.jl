"""Check saved complete transport outputs, worker evidence and conservative timings."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from controlled_reference import validate_workers


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(folder, core_root, driver_root, prior_correctness):
    native = np.load(folder / "native.npz")
    reference = np.load(folder / "reference/reference.npz")
    receipt = json.loads((folder / "reference/summary.json").read_text())
    exact = json.loads((folder / "exact-source/summary.json").read_text())
    prior = json.loads(prior_correctness.read_text())
    fields = ("T", "conductivity_parallel", "conductivity_serial", "viscosity_parallel", "viscosity_serial")
    metadata = {key[:-5]: bytes(native[key]).decode() for key in native.files if key.endswith("_utf8")}
    checks, errors = {}, {}
    for key in fields:
        a, b = native[key], reference[key]
        checks[key] = bool(a.shape == (5000,) and b.shape == a.shape and np.isfinite(a).all()
                           and np.allclose(a, b, rtol=1e-10, atol=1e-15))
        errors[key] = float(np.max(np.abs(a-b)))
    checks["native_internal_checks"] = bool(native["native_checks_pass"].shape == (1,) and native["native_checks_pass"][0])
    checks["reference_internal_checks"] = receipt["passed"] and all(receipt["checks"].values())
    checks["exact_script_checks"] = exact["passed"] and len(exact["runs"]) == 11
    checks["source_and_mechanism"] = (receipt["measured_file_sha256"]["multiprocessing_viscosity.py"] == prior["source_sha256"]
        == exact["hashes"]["source"] and metadata["mechanism_sha256"] == prior["mechanism_sha256"] == exact["hashes"]["mechanism"])
    checks["native_inputs"] = metadata["sidecar_sha256"] == prior["sidecar_sha256"] and metadata["transport_sha256"] == prior["multicomponent_sha256"]
    checks["pinned_cantera_build"] = receipt["cantera_loaded_sha256"] == prior["cantera_loaded_sha256"] and receipt["cantera_version"] == "4.0.0a2"
    checks["native_threads"] = int(native["julia_threads"][0]) == 4 and int(native["blas_threads"][0]) == 1
    checks["host_match"] = all(metadata[key] == receipt["host"][key] for key in ("cpu", "system", "kernel_release"))
    checks["target_host"] = ((receipt["host"]["system"] == "Linux" and "microsoft" in receipt["host"]["kernel_release"].lower())
        or (receipt["host"]["system"] == "Darwin" and "Apple M4" in receipt["host"]["cpu"] and int(native["accelerate_threading_mode"][0]) == 1))
    core_files = sorted(p.relative_to(core_root).as_posix() for p in core_root.rglob("*") if p.is_file())
    core_hashes = {p: sha(core_root / p) for p in core_files}
    checks["tested_core_bytes"] = (metadata["core_source_paths"] == "\n".join(core_files)
                                      and metadata["core_source_sha256"] == "\n".join(core_hashes.values()))
    checks["driver_bytes"] = (metadata["driver_sha256"] == sha(driver_root / "controlled_native.jl")
        and receipt["measured_file_sha256"]["controlled_reference.py"] == sha(driver_root / "controlled_reference.py")
        and exact["hashes"]["driver"] == sha(driver_root / "exact_source_timing.py"))
    checks["native_example_bytes"] = metadata["example_sha256"] == sha(core_root.parent / "example/transport/multiprocessing_viscosity.jl")
    worker_lifecycles = 0
    for run in receipt["runs"]:
        saved = np.load(folder / "reference" / run["label"] / "arrays.npz")
        if not all(np.array_equal(saved[k], reference[k]) for k in fields):
            raise ValueError("saved reference arrays fail exact replay")
        if run["controlled"]:
            assert len(run["pools"]) == 2
            for pool in run["pools"]:
                validate_workers(pool["before"], receipt["cantera_loaded_sha256"])
                validate_workers(pool["after"], receipt["cantera_loaded_sha256"])
                assert pool["before"] == pool["after"] and pool["map_calls"] == 1
                for stage in ("before", "after"):
                    window = pool[stage + "_inspection"]
                    expected_start = max(window["start"], max(r["initializer_completed"] for r in pool[stage])) if stage == "before" else window["start"]
                    assert window["excluded_start"] == expected_start
                    assert window["end"] > expected_start
                    assert window["excluded_seconds"] == window["end"]-expected_start
                worker_lifecycles += 4
            pool_exclusion = sum(pool[s + "_inspection"]["excluded_seconds"] for pool in run["pools"] for s in ("before", "after"))
            assert np.isclose(sum(run["excluded_seconds"]), pool_exclusion, rtol=0, atol=1e-12)
        assert np.array_equal(np.array(run["raw_seconds"])-run["excluded_seconds"], run["adjusted_seconds"])
    checks["worker_lifecycles_and_exclusions"] = worker_lifecycles == 88
    checks["all_reference_arrays_replayed"] = len(receipt["runs"]) == 22
    original, controlled, exact_times = [], [], []
    for i in range(9):
        label = f"warm-{i:02d}"
        row = {r["label"]: r for r in receipt["runs"]}
        original.append(sum(row[label+"-original"]["raw_seconds"]))
        controlled.append(sum(row[label+"-controlled"]["adjusted_seconds"]))
        source = next(r for r in exact["runs"] if r["label"] == label)
        assert source["exit_code"] == 0 and len(source["reported_seconds"]) == 4
        assert source["conservative_seconds"] == [max(0,x-0.0005) for x in source["reported_seconds"]]
        exact_times.append(sum(source["conservative_seconds"]))
    native_seconds = native["warm_seconds"]
    reference_seconds = np.minimum(np.minimum(original, controlled), exact_times)
    checks["timing_vectors"] = bool(native_seconds.shape == (9,) and np.isfinite(native_seconds).all()
        and (native_seconds > 0).all() and np.isfinite(reference_seconds).all() and (reference_seconds > 0).all())
    ratio = float(np.median(reference_seconds) / np.median(native_seconds))
    passed = all(checks.values())
    report = dict(example="transport/multiprocessing_viscosity", passed=passed,
        performance_qualified=bool(passed and ratio >= 0.95), minimum_speed_ratio=0.95,
        speed_ratio=ratio, speed_ratio_definition="median(Cantera seconds) / median(Julia seconds)",
        scope="Four complete 5000-temperature property calls including fresh phases/workspaces and four-worker pool creation/teardown. One first call, one full warmup, nine complete warm repetitions. Julia compilation/imports are outside warm timing. Plotting and output file writes are excluded.",
        reference_rule="Per repetition, minimum of verified-worker adjusted total, uninstrumented original-function total and unchanged source-script timer total. Source-script timers are reduced by 0.0005 s per printed value for rounding. Worker audit subtraction retains initialization through the latest original initializer return timestamp.",
        source_url=prior["source_url"], source_sha256=prior["source_sha256"], cantera_version=receipt["cantera_version"],
        cantera_loaded_sha256=receipt["cantera_loaded_sha256"], native_example_sha256=metadata["example_sha256"],
        native_core_sha256=core_hashes, input_sha256={key:metadata[key+"_sha256"] for key in ("mechanism","sidecar","transport")},
        measured_driver_sha256={"native":metadata["driver_sha256"], "controlled_reference":receipt["measured_file_sha256"]["controlled_reference.py"], "exact_source":exact["hashes"]["driver"]},
        host=receipt["host"], julia_version=metadata["julia_version"], julia_threads=4, source_processes=4, numerical_threads=1,
        worker_lifecycles_verified=worker_lifecycles, multiprocessing_start_method=receipt["multiprocessing_start_method"],
        uninstrumented_child_settings_unverified=True,
        initializer_runtime_settings_unverified=True,
        worker_verification_scope="Controlled workers after original initializer return and before/after the original property map. Unwrapped function calls and exact-script runs are unverified timing caps.",
        exact_script_output_arrays_checked=False,
        exact_script_verification="Pinned source bytes, successful execution and all four original printed timers per call; output-array replay comes from the controlled/unwrapped reference calls.",
        correctness={"rtol":1e-10,"atol":1e-15,"maximum_absolute_errors":errors}, checks=checks,
        native_first_call_seconds=float(native["first_call_seconds"][0]), native_warm_seconds=native_seconds.tolist(),
        original_function_warm_seconds=original, controlled_adjusted_warm_seconds=controlled,
        controlled_raw_property_seconds=[r["raw_seconds"] for r in receipt["runs"] if r["controlled"] and r["label"].startswith("warm-")],
        controlled_excluded_property_seconds=[r["excluded_seconds"] for r in receipt["runs"] if r["controlled"] and r["label"].startswith("warm-")],
        exact_source_rounded_lower_bound_warm_seconds=exact_times, reference_conservative_warm_seconds=reference_seconds.tolist(),
        native_median_seconds=float(np.median(native_seconds)), reference_median_seconds=float(np.median(reference_seconds)),
        evidence_sha256={name:sha(folder / path) for name,path in {"native":"native.npz","reference":"reference/reference.npz",
            "worker_and_timing_receipts":"reference/summary.json","exact_source_receipts":"exact-source/summary.json"}.items()})
    (folder / "paired-summary.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({k:report[k] for k in ("passed","performance_qualified","speed_ratio","native_median_seconds","reference_median_seconds","checks")},indent=2))
    if not passed:
        raise SystemExit("saved transport checks failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("folder", "core-root", "driver-root", "prior-correctness"):
        parser.add_argument("--"+name,type=Path,required=True)
    args = parser.parse_args()
    summarize(args.folder,args.core_root,args.driver_root,args.prior_correctness)
