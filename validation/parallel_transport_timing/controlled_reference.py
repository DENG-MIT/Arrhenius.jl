"""Complete-source timing with verified worker runtimes and conservative audit exclusion.

The original Pool constructor, initializer, map, predicates and teardown run.
We intercept the Pool factory, initializer return and map boundary. Four blocking
inspection jobs prove four distinct live workers before and after each original
map. Inspection time after all original initializers finish is subtracted. The
smaller of adjusted and uninstrumented totals is the reference, retaining the
original phase/pool startup without giving audit overhead a timing advantage.
"""
import argparse
import hashlib
import importlib
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time

import numpy as np

SOURCE_SHA = "6dd426c24a6af5a5b9d2e87d2d9b33e78de83645952b946e81d34e727debc1de"
MECHANISM_SHA = "06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345"
NWORKERS, NPOINTS = 4, 5000
REAL_POOL = multiprocessing.Pool
WORKER_TIMEOUT = 30.0
INITIALIZER_COMPLETED = None


def record_initializer(initializer, initargs):
    """Call the unchanged original initializer, then record its return time."""
    global INITIALIZER_COMPLETED
    initializer(*initargs)
    INITIALIZER_COMPLETED = time.perf_counter()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf8")


def inspect_worker(root, source, mechanism, directory, set_accelerate):
    """Do not touch transport properties or TPX; retain the original initialized phase."""
    root, directory = Path(root), Path(directory)
    environment = load_module("worker_environment", root / "validation/benchmark_environment.py")
    cases = load_module("worker_cases", root / "validation/parallel_transport_cases.py")
    original = importlib.import_module(Path(source).stem)
    phase = original.gases[mechanism]
    report = dict(pid=os.getpid(), initializer_completed=INITIALIZER_COMPLETED,
                  source_sha256=sha(source), mechanism_sha256=sha(mechanism),
                  source_path_matches=Path(original.__file__).resolve() == Path(source).resolve(),
                  transport_model=phase.transport_model,
                  numerical_threads=environment.verify_numerical_threads(set_accelerate=set_accelerate),
                  cantera_loaded_sha256=cases.loaded_cantera_hashes(environment))
    # Rename a complete unique file so the parent never sees partial JSON.
    temporary = directory / f"{os.getpid()}.tmp"
    save_json(temporary, report)
    temporary.rename(directory / f"{os.getpid()}.json")
    deadline = time.monotonic() + WORKER_TIMEOUT
    while not (directory / "release").exists():
        if time.monotonic() >= deadline:
            raise TimeoutError("worker inspection release was not received")
        time.sleep(0.001)
    return report


class ControlledPool:
    def __init__(self, audit, *args, **kwargs):
        self.audit = audit
        # Constructor elapsed is included; do not hoist phase/pool construction.
        if args or set(kwargs) != {"processes", "initializer", "initargs"} or kwargs["processes"] != NWORKERS:
            raise RuntimeError("unexpected original Pool construction")
        self.pool = REAL_POOL(processes=kwargs["processes"], initializer=record_initializer,
                              initargs=(kwargs["initializer"], kwargs["initargs"]))
        self.number = len(audit["pools"])
        self.record = {"before": None, "after": None, "map_calls": 0}
        audit["pools"].append(self.record)
        try:
            self._inspect("before", True)
        except BaseException:
            self.pool.terminate()
            self.pool.join()
            raise

    def _inspect(self, stage, set_accelerate):
        start = time.perf_counter()
        excluded_start = start
        directory = self.audit["directory"] / f"pool-{self.number:02d}-{stage}"
        directory.mkdir()
        try:
            jobs = [self.pool.apply_async(inspect_worker,
                    (str(self.audit["root"]), str(self.audit["source"]), self.audit["mechanism"],
                     str(directory), set_accelerate)) for _ in range(NWORKERS)]
            deadline = time.monotonic() + WORKER_TIMEOUT
            while len(list(directory.glob("*.json"))) != NWORKERS:
                # Surface child errors immediately rather than waiting for missing receipts.
                for job in jobs:
                    if job.ready():
                        job.get(timeout=0)
                if time.monotonic() >= deadline:
                    raise TimeoutError("four distinct initialized workers were not inspected")
                time.sleep(0.001)
            (directory / "release").touch()
            records = sorted((job.get(timeout=WORKER_TIMEOUT) for job in jobs), key=lambda r: r["pid"])
            validate_workers(records, self.audit["libraries"])
            if stage == "before":
                excluded_start = max(start, max(r["initializer_completed"] for r in records))
                if excluded_start > time.perf_counter():
                    raise RuntimeError("inconsistent cross-process initialization clock")
            self.record[stage] = records
            if stage == "after":
                if records != self.record["before"]:
                    raise RuntimeError("worker identity, source, libraries or numerical settings changed")
        finally:
            (directory / "release").touch()
            end = time.perf_counter()
            self.record[stage + "_inspection"] = dict(start=start, end=end,
                  excluded_start=excluded_start, excluded_seconds=end-excluded_start,
                  retained_initializer_wait_seconds=excluded_start-start)
            self.audit["excluded_seconds"] += end - excluded_start

    def __enter__(self):
        self.pool.__enter__()
        return self

    def map(self, *args, **kwargs):
        self.record["map_calls"] += 1
        if self.record["map_calls"] != 1:
            raise RuntimeError("pinned source unexpectedly issued more than one map")
        outputs = self.pool.map(*args, **kwargs)
        self._inspect("after", False)
        return outputs

    def __exit__(self, *args):
        return self.pool.__exit__(*args)


def validate_workers(records, libraries):
    if len(records) != NWORKERS or len({r["pid"] for r in records}) != NWORKERS:
        raise RuntimeError("four distinct worker receipts required")
    for row in records:
        if not isinstance(row["initializer_completed"], float) or row["initializer_completed"] <= 0:
            raise RuntimeError("missing original initializer completion timestamp")
        if not (row["source_path_matches"] and row["source_sha256"] == SOURCE_SHA
                and row["mechanism_sha256"] == MECHANISM_SHA
                and row["transport_model"] == "multicomponent"
                and row["cantera_loaded_sha256"] == libraries):
            raise RuntimeError("worker provenance or original phase mismatch")
        threads = row["numerical_threads"]
        if not threads["threadpools"] or any(p["threads"] != 1 for p in threads["threadpools"]):
            raise RuntimeError("worker numerical runtime was not single threaded")
        if "accelerate_threading_mode" in threads and threads["accelerate_threading_mode"] != 1:
            raise RuntimeError("worker Accelerate mode mismatch")


def sequence(original, mechanism):
    return (
        ("conductivity_parallel", lambda: original.parallel(mechanism, original.get_thermal_conductivity, NWORKERS, NPOINTS)),
        ("conductivity_serial", lambda: original.serial(mechanism, original.get_thermal_conductivity, NPOINTS)),
        ("viscosity_parallel", lambda: original.parallel(mechanism, original.get_viscosity, NWORKERS, NPOINTS)),
        ("viscosity_serial", lambda: original.serial(mechanism, original.get_viscosity, NPOINTS)),
    )


def complete_run(original, args, label, controlled, libraries):
    directory = args.output_dir / label
    directory.mkdir()
    audit = dict(directory=directory, root=args.source_root, source=args.source_example,
                 mechanism=str(args.mechanism), libraries=libraries, excluded_seconds=0.0, pools=[])
    raw_times, excluded_times, outputs = [], [], {"T": np.linspace(300, 900, NPOINTS)}
    if controlled:
        multiprocessing.Pool = lambda *a, **k: ControlledPool(audit, *a, **k)
    try:
        for name, call in sequence(original, str(args.mechanism)):
            before_excluded = audit["excluded_seconds"]
            start = time.perf_counter()
            raw = call()
            elapsed = time.perf_counter() - start
            raw_times.append(elapsed)
            excluded_times.append(audit["excluded_seconds"] - before_excluded)
            outputs[name] = np.asarray(raw)
        np.savez(directory / "arrays.npz", **outputs)
        report = dict(label=label, controlled=controlled, raw_seconds=raw_times,
                      excluded_seconds=excluded_times,
                      adjusted_seconds=(np.array(raw_times)-excluded_times).tolist(), pools=audit["pools"])
        save_json(directory / "timing.json", report)
        if controlled and (len(audit["pools"]) != 2 or any(
                p["before"] is None or p["after"] is None or p["map_calls"] != 1 for p in audit["pools"])):
            raise RuntimeError("incomplete worker lifecycle audit")
        if not all(x > 0 for x in report["adjusted_seconds"]):
            raise RuntimeError("invalid audit subtraction")
        return report, outputs
    except BaseException as exc:
        save_json(directory / "failure.json", dict(error=repr(exc), pools=audit["pools"],
                  raw_seconds=raw_times, excluded_seconds=excluded_times))
        raise
    finally:
        multiprocessing.Pool = REAL_POOL


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source-example", "source-root", "mechanism", "output-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=9)
    args = parser.parse_args()
    for field in ("source_example", "source_root", "mechanism", "output_dir"):
        setattr(args, field, getattr(args, field).resolve())
    if args.repetitions < 1 or sha(args.source_example) != SOURCE_SHA or sha(args.mechanism) != MECHANISM_SHA:
        parser.error("invalid repetitions or pinned source/input hashes")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    import cantera as ct
    if ct.__version__ != "4.0.0a2":
        parser.error("expected pinned Cantera 4.0.0a2 reference")
    environment = load_module("benchmark_environment", args.source_root / "validation/benchmark_environment.py")
    cases = load_module("parallel_transport_cases", args.source_root / "validation/parallel_transport_cases.py")
    parent_threads = environment.verify_numerical_threads(set_accelerate=True)
    libraries = cases.loaded_cantera_hashes(environment)
    paths = [args.source_example, args.mechanism, Path(__file__),
             args.source_root / "validation/benchmark_environment.py",
             args.source_root / "validation/parallel_transport_cases.py"]
    hashes = {p.name: sha(p) for p in paths}
    sys.path.insert(0, str(args.source_example.parent))
    original = importlib.import_module(args.source_example.stem)
    if Path(original.__file__).resolve() != args.source_example:
        parser.error("unexpected source module path")
    runs, reference = [], None
    # Alternate group order across repetitions; neither group reuses pools/phases.
    for label in ["cold", "warmup"] + [f"warm-{i:02d}" for i in range(args.repetitions)]:
        order = [False, True] if len(runs) % 4 == 0 else [True, False]
        for controlled in order:
            name = label + ("-controlled" if controlled else "-original")
            report, outputs = complete_run(original, args, name, controlled, libraries)
            runs.append(report)
            if reference is None:
                reference = outputs
            if not all(np.array_equal(outputs[k], reference[k]) for k in reference):
                raise RuntimeError(f"{name} does not exactly replay original arrays; saved in run directory")
            print(f"PASS {name}: adjusted={sum(report['adjusted_seconds']):.6f}s", flush=True)
    checks = dict(replays_exact=True, worker_lifecycles_verified=True,
                  source_hashes_unchanged=hashes == {p.name: sha(p) for p in paths},
                  libraries_unchanged=libraries == cases.loaded_cantera_hashes(environment),
                  parent_threads_unchanged=parent_threads == environment.verify_numerical_threads(set_accelerate=False),
                  finite_positive=all(v.shape == (NPOINTS,) and np.isfinite(v).all() and (v > 0).all()
                                      for v in reference.values()),
                  serial_parallel_equal=all(np.array_equal(reference[k + "_parallel"], reference[k + "_serial"])
                                            for k in ("conductivity", "viscosity")))
    np.savez(args.output_dir / "reference.npz", **reference)
    summary = dict(example="transport/multiprocessing_viscosity", checks=checks, passed=all(checks.values()),
                   performance_qualified=False, repetitions=args.repetitions, source_processes=NWORKERS,
                   temperatures=NPOINTS, cantera_version=ct.__version__, measured_file_sha256=hashes,
                   cantera_loaded_sha256=libraries, parent_numerical_threads=parent_threads,
                   multiprocessing_start_method=multiprocessing.get_start_method(), host=environment.host_metadata(),
                   uninstrumented_child_settings_unverified=True,
                   timing_scope="All four original calls; pool and phase construction retained through the last original initializer return timestamp. Subsequent worker audit windows subtracted. Use minimum of original and controlled totals per repetition to prevent audit overhead giving the reference a timing disadvantage.",
                   runs=runs)
    save_json(args.output_dir / "summary.json", summary)
    if not summary["passed"]:
        raise RuntimeError("controlled reference verification failed")


if __name__ == "__main__":
    main()
