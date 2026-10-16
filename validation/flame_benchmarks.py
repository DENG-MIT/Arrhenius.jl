"""Correctness and benchmark harness: native Julia flames vs. exact Cantera 4.0.

Orchestrates the two published premixed-flame examples (free flame and
burner-stabilized, plus the optional prescribed-temperature burner), timing
identical scopes on both sides: mechanism/sidecar loading is outside every
timer, flame construction + initialization + adaptive solve are inside. The
first repetition is recorded separately (Julia JIT / Cantera warm-up), followed
by at least five warm repetitions; all raw seconds, medians, grid sizes, and
outputs are retained in the JSON report.

Accuracy criteria (never relaxed to make a case pass):
  * flame speed relative error          <= 1%   (free flame)
  * maximum temperature relative error  <= 1%   (all cases)
  * per-species profile peak error      <= 5% of the species' reference peak,
    with an absolute floor of 1e-7 on the normalizer
  * inlet/outlet elemental mass conservation of the Julia solution <= 1e-6
Free-flame profiles are aligned by temperature phase to remove the arbitrary
translational shift; burner profiles are compared on the physical coordinate.
All errors are reported regardless of gate outcomes.

Performance gate: Cantera median / Julia median >= 0.95, enforced only for a
formal pass (--formal, which verifies --target, default WSL); otherwise
the ratio is reported as informational. Comparisons require Cantera 4.0.x.
If the Cantera build does not expose its source commit, an explicit
--cantera-commit SHA override is required and recorded as provided provenance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path
from benchmark_environment import host_metadata, matches_target, cantera_library_hashes

for _thread_variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_thread_variable] = "1"
import numpy as np

SPEED_TOL = 0.01
TMAX_TOL = 0.01
PROFILE_TOL = 0.05
PROFILE_FLOOR = 1e-7
ELEMENT_TOL = 1e-6
MIN_SPEED_RATIO = 0.95
PROFILE_POINTS = 2000

FIXED_POSITIONS = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 1.0]
FIXED_TEMPERATURES = [373.0, 650.0, 1000.0, 1350.0, 1650.0, 1750.0, 1750.0]

CASES = {
    "freeflame": {
        "T": 300.0, "P_atm": 1.0, "X": "H2:1.1,O2:1,AR:5", "width": 0.03,
        "ratio": 3.0, "slope": 0.06, "curve": 0.12, "kind": "free",
    },
    "burner": {
        "T": 373.0, "P_atm": 0.05, "X": "H2:1.5,O2:1,AR:7", "width": 0.5,
        "mdot": 0.06, "ratio": 3.0, "slope": 0.05, "curve": 0.1, "kind": "burner",
    },
    "burner-fixed": {
        "T": 373.0, "P_atm": 0.05, "X": "H2:1.5,O2:1,AR:7", "width": 0.5,
        "mdot": 0.06, "ratio": 3.0, "slope": 0.05, "curve": 0.1, "kind": "burner-fixed",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_cantera_case(ct, gas, case: dict, reps: int) -> dict:
    """Time construction + initialization + solve, first run then warm reps."""
    seconds, speeds, tmaxes, points = [], [], [], []
    flame = None
    for _ in range(reps + 1):
        start = time.perf_counter()
        gas.TPX = case["T"], case["P_atm"] * ct.one_atm, case["X"]
        if case["kind"] == "free":
            flame = ct.FreeFlame(gas, width=case["width"])
            flame.transport_model = "mixture-averaged"
            flame.set_refine_criteria(
                ratio=case["ratio"], slope=case["slope"], curve=case["curve"])
            flame.solve(loglevel=0, auto=True)
        elif case["kind"] == "burner":
            flame = ct.BurnerFlame(gas, width=case["width"])
            flame.burner.mdot = case["mdot"]
            flame.set_refine_criteria(
                ratio=case["ratio"], slope=case["slope"], curve=case["curve"])
            flame.solve(loglevel=0, auto=True)
        else:
            flame = ct.BurnerFlame(gas, width=case["width"])
            flame.burner.T = case["T"]
            flame.burner.X = case["X"]
            flame.burner.mdot = case["mdot"]
            flame.energy_enabled = False
            flame.flame.set_fixed_temp_profile(FIXED_POSITIONS, FIXED_TEMPERATURES)
            flame.set_refine_criteria(
                ratio=case["ratio"], slope=case["slope"], curve=case["curve"])
            flame.solve(loglevel=0, auto=False)
        seconds.append(time.perf_counter() - start)
        speeds.append(float(flame.velocity[0]))
        tmaxes.append(float(max(flame.T)))
        points.append(len(flame.grid))
    return {
        "seconds": seconds, "speed": speeds, "tmax": tmaxes, "points": points,
        "grid": np.asarray(flame.grid), "T": np.asarray(flame.T),
        "Y": np.asarray(flame.Y), "velocity": np.asarray(flame.velocity),
        "inlet_Y": np.asarray(flame.inlet.Y if case["kind"] == "free" else flame.burner.Y),
        "species_names": list(gas.species_names),
    }


def species_order(names_a, names_b):
    idx = []
    for name in names_a:
        if name not in names_b:
            raise ValueError(f"species {name} missing from reference")
        idx.append(names_b.index(name))
    return idx


def profile_errors(ref, jul, kind: str):
    """Per-species normalized peak errors with an absolute floor."""
    idx = species_order(jul["species_names"], ref["species_names"])
    refY = ref["Y"][idx, :]
    ref_x, jul_x = ref["grid"], jul["grid"]
    if kind == "free":
        # Align a single temperature crossing by translation. Preserve flame
        # thickness and nonmonotonic temperature/species structure.
        level = .5*(max(ref["T"][0],jul["T"][0])+min(max(ref["T"]),max(jul["T"])))
        def crossing(z, T):
            j = int(np.flatnonzero(T >= level)[0])
            if j == 0:
                raise ValueError("temperature profile has no interior alignment crossing")
            return z[j-1]+(level-T[j-1])/(T[j]-T[j-1])*(z[j]-z[j-1])
        ref_x = ref_x-crossing(ref_x,ref["T"])
        jul_x = jul_x-crossing(jul_x,jul["T"])
    lo,hi = max(ref_x[0],jul_x[0]),min(ref_x[-1],jul_x[-1])
    grid = np.unique(np.r_[np.linspace(lo,hi,PROFILE_POINTS),ref_x,jul_x])
    grid = grid[(grid>=lo)&(grid<=hi)]
    ref_vals, jul_vals = refY, jul["Y"]
    temperature_error = float(np.max(np.abs(np.interp(grid,ref_x,ref["T"])-np.interp(grid,jul_x,jul["T"]))))
    temperature = {"max_abs_error_K":temperature_error,
                   "normalized_peak_error":temperature_error/max(ref["T"])}
    species = {}
    for k, name in enumerate(jul["species_names"]):
        r = np.interp(grid, ref_x, ref_vals[k])
        j = np.interp(grid, jul_x, jul_vals[k])
        peak = float(np.max(np.abs(r)))
        err = float(np.max(np.abs(r - j))) / max(peak, PROFILE_FLOOR)
        species[name] = {"peak_reference": peak, "normalized_peak_error": err,
                         "pass": err <= PROFILE_TOL}
    return {"species": species, "temperature": temperature,
            "all_species_pass": all(s["pass"] for s in species.values())}


def element_conservation(gas, jul: dict):
    """Inlet vs outlet elemental mass fractions of the Julia solution."""
    n_elements = gas.n_elements
    idx = species_order(jul["species_names"], list(gas.species_names))
    comp = np.array([[gas.n_atoms(k,m)*gas.atomic_weights[m]/gas.molecular_weights[k]
                      for m in range(n_elements)]
                     for k in range(gas.n_species)])[idx, :]
    names = [gas.element_name(m) for m in range(n_elements)]
    inlet = comp.T @ jul["inlet_Y"]
    outlet = comp.T @ jul["Y"][:, -1]
    diff = np.abs(inlet - outlet)
    return {name: {"inlet": float(i), "outlet": float(o),
                   "abs_error": float(d), "pass": bool(d <= ELEMENT_TOL)}
            for name, i, o, d in zip(names, inlet, outlet, diff)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mechanism", type=Path, help="YAML mechanism with a .npz sidecar")
    parser.add_argument("--output", type=Path, default=Path("flame-benchmark-results"),
                        help="directory for the JSON report and raw NPZ data")
    parser.add_argument("--julia", default="julia", help="Julia executable")
    parser.add_argument("--project", type=Path,
                        default=Path(__file__).resolve().parent.parent,
                        help="Julia project (defaults to this repo)")
    parser.add_argument("--reps", type=int, default=5,
                        help="warm repetitions per case and per side (>=5)")
    parser.add_argument("--cases", default=",".join(CASES),
                        help="comma-separated subset of " + ",".join(CASES))
    parser.add_argument("--cantera-commit", default=None,
                        help="verified Cantera source SHA when the build does "
                             "not report one (recorded as provided provenance)")
    parser.add_argument("--formal", action="store_true",
                        help="formal performance pass: enforce target host and "
                             "the >=0.95 Cantera/Julia median speed ratio")
    parser.add_argument("--target", choices=("wsl", "apple-m4"), default="wsl")
    args = parser.parse_args()
    if args.reps < 5:
        parser.error("--reps must be at least 5")

    import cantera as ct

    errors = []
    report = {"date": time.strftime("%Y-%m-%d %H:%M:%S %z"),
              "criteria": {"speed_rel_tol": SPEED_TOL, "tmax_rel_tol": TMAX_TOL,
                           "profile_peak_rel_tol": PROFILE_TOL,
                           "profile_abs_floor": PROFILE_FLOOR,
                           "element_abs_tol": ELEMENT_TOL,
                           "min_speed_ratio": MIN_SPEED_RATIO},
              "hardware": host_metadata(), "benchmark_target": args.target,
              "threads": {"julia_blas": 1, "julia": 1, "omp": 1}}

    if not ct.__version__.startswith("4.0"):
        parser.error(f"comparisons require Cantera 4.0.x, found {ct.__version__}")
    commit = getattr(ct, "__git_commit__", None)
    provenance = "build-reported"
    if not commit or str(commit).lower() == "unknown":
        if not args.cantera_commit:
            parser.error("this Cantera build does not report its source commit; "
                         "pass --cantera-commit SHA (recorded as provided provenance)")
        commit, provenance = args.cantera_commit, "provided"
    report["versions"] = {"cantera": ct.__version__,
                          "cantera_commit": commit,
                          "cantera_commit_provenance": provenance,
                          "numpy": np.__version__, "python": platform.python_version()}

    if args.formal and not matches_target(report["hardware"],args.target):
        parser.error(f"a formal performance pass requires {args.target}; "
                     f"detected: {report['hardware']}")
    if not args.formal:
        report["note"] = "correctness gates enforced; speed ratio informational only"
    import cantera._cantera as compiled
    report["versions"]["cantera_extension_sha256"] = sha256(Path(compiled.__file__))
    report["versions"]["cantera_shared_libraries_sha256"] = cantera_library_hashes(ct.__file__)

    try:
        source_commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=args.project,
            capture_output=True, text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "status", "--porcelain"], cwd=args.project,
            capture_output=True, text=True, check=True).stdout.strip())
    except Exception:
        source_commit, dirty = None, None
    report["source"] = {"commit": source_commit, "dirty": dirty,
                        "mechanism": str(args.mechanism),
                        "mechanism_sha256": sha256(args.mechanism),
                        "sidecar_sha256": sha256(Path(str(args.mechanism)+".npz")),
                        "julia_sources_sha256": {str(p.relative_to(args.project)):sha256(p)
                            for p in sorted((args.project/"src").rglob("*.jl"))}}

    cases = args.cases.split(",")
    for name in cases:
        if name not in CASES:
            parser.error(f"unknown case {name}; choose from {','.join(CASES)}")

    args.output.mkdir(parents=True, exist_ok=True)

    # --- Cantera reference: mechanism load outside every timed region ---
    ct_results = {}
    for name in cases:
        case = CASES[name]
        gas = ct.Solution(str(args.mechanism))
        gas.TPX = case["T"], case["P_atm"] * ct.one_atm, case["X"]
        ct_results[name] = run_cantera_case(ct, gas, case, args.reps)
        np.savez(args.output / f"cantera-{name}.npz",
                 grid=ct_results[name]["grid"], T=ct_results[name]["T"],
                 Y=ct_results[name]["Y"], velocity=ct_results[name]["velocity"])
        print(f"cantera {name}: first={ct_results[name]['seconds'][0]:.3f}s "
              f"warm_median={statistics.median(ct_results[name]['seconds'][1:]):.3f}s "
              f"points={ct_results[name]['points'][-1]}")

    # --- Julia native runs: same timing scope, single-threaded ---
    julia_npz = args.output / "julia-results.npz"
    env = dict(os.environ, JULIA_NUM_THREADS="1", OMP_NUM_THREADS="1")
    script = Path(__file__).with_name("flame_benchmarks.jl")
    cmd = [args.julia, f"--project={args.project}", str(script),
           str(args.mechanism), str(julia_npz),
           f"--reps={args.reps}", f"--cases={args.cases}"]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        errors.append(f"julia runner exited with code {proc.returncode}")
        jl_raw = None
    else:
        jl_raw = np.load(julia_npz, allow_pickle=False)
        report["julia_host"] = {key: bytes(jl_raw[key+"_utf8"]).decode()
                                for key in ("system","kernel_release","cpu")}
        if args.formal and (not matches_target(report["julia_host"],args.target)
                           or report["julia_host"]["kernel_release"] != report["hardware"]["kernel_release"]):
            errors.append("Julia and Cantera must both run on the selected target host")
        report["versions"]["julia"] = subprocess.run(
            [args.julia, "--version"], capture_output=True, text=True
        ).stdout.strip()

    # --- Comparisons ---
    report["cases"] = {}
    for name in cases:
        case = CASES[name]
        ct_res = ct_results[name]
        entry = {"cantera": {"first_seconds": ct_res["seconds"][0],
                             "warm_seconds": ct_res["seconds"][1:],
                             "warm_median_seconds":
                                 statistics.median(ct_res["seconds"][1:]),
                             "speed": ct_res["speed"][-1],
                             "tmax": ct_res["tmax"][-1],
                             "points": ct_res["points"][-1]},
                 "errors": {}, "gates": {}}
        ct_warm_median = entry["cantera"]["warm_median_seconds"]
        if jl_raw is not None:
            key = name.replace("-", "_")
            jl_seconds = [float(s) for s in jl_raw[key + "_seconds"]]
            jl = {"seconds": jl_seconds,
                  "speed": float(jl_raw[key + "_speed"][-1]),
                  "tmax": float(jl_raw[key + "_tmax"][-1]),
                  "points": int(jl_raw[key + "_points"][-1]),
                  "grid": jl_raw[key + "_grid"], "T": jl_raw[key + "_T"],
                  "Y": jl_raw[key + "_Y"], "velocity": jl_raw[key + "_velocity"],
                  "inlet_Y": jl_raw[key + "_inlet_Y"],
                  "species_names":
                      jl_raw["species_names_utf8"].tobytes().decode().split("\n")}
            jl_warm = jl_seconds[1:]
            entry["julia"] = {"first_seconds": jl_seconds[0],
                              "warm_seconds": jl_warm,
                              "warm_median_seconds": statistics.median(jl_warm),
                              "speed": jl["speed"], "tmax": jl["tmax"],
                              "points": jl["points"]}
            entry["speed_ratio_cantera_over_julia"] = (
                ct_warm_median / statistics.median(jl_warm))

            speed_err = abs(jl["speed"] - ct_res["speed"][-1]) / abs(ct_res["speed"][-1])
            tmax_err = abs(jl["tmax"] - ct_res["tmax"][-1]) / abs(ct_res["tmax"][-1])
            entry["errors"]["speed_rel"] = speed_err
            entry["errors"]["tmax_rel"] = tmax_err
            if case["kind"] == "free":
                entry["gates"]["speed_rel<=1%"] = speed_err <= SPEED_TOL
            entry["gates"]["tmax_rel<=1%"] = tmax_err <= TMAX_TOL

            profiles = profile_errors(ct_res, jl, case["kind"])
            entry["errors"]["profiles"] = profiles
            entry["gates"]["species_peak<=5%"] = profiles["all_species_pass"]
            entry["gates"]["temperature_profile<=1%"] = profiles["temperature"]["normalized_peak_error"] <= TMAX_TOL

            elements = element_conservation(gas, jl)
            entry["errors"]["element_conservation"] = elements
            entry["errors"]["cantera_element_conservation"] = element_conservation(gas,ct_res)
            entry["gates"]["outlet_element_fraction_drift<=1e-6"] = all(
                e["pass"] for e in elements.values())
        report["cases"][name] = entry

    # --- Gates ---
    for name, entry in report["cases"].items():
        for gate, ok in entry["gates"].items():
            if not ok:
                errors.append(f"{name}: gate {gate} failed")
        ratio = entry.get("speed_ratio_cantera_over_julia")
        if ratio is not None:
            entry["gates"]["speed_ratio>=0.95"] = ratio >= MIN_SPEED_RATIO
            if args.formal and ratio < MIN_SPEED_RATIO:
                errors.append(f"{name}: speed ratio {ratio:.3f} < {MIN_SPEED_RATIO} (formal)")

    report["formal"] = args.formal
    report["passed"] = not errors
    if errors:
        report["failures"] = errors

    report_path = args.output / "flame_benchmarks_report.json"
    report_path.write_text(json.dumps(report, indent=2,
        default=lambda obj: obj.item() if isinstance(obj,np.generic) else obj.tolist()) + "\n")

    # --- Concise summary ---
    print("\n=== flame benchmark summary ===")
    for name, entry in report["cases"].items():
        line = f"{name}: "
        if "julia" in entry:
            line += (f"ratio={entry['speed_ratio_cantera_over_julia']:.3f} "
                     f"speed_err={entry['errors']['speed_rel']:.2e} "
                     f"tmax_err={entry['errors']['tmax_rel']:.2e} "
                     f"gates={'PASS' if all(entry['gates'].values()) else 'FAIL'}")
        else:
            line += "julia results unavailable"
        print(line)
    print(f"overall: {'PASS' if report['passed'] else 'FAIL'} "
          f"({'formal '+args.target if args.formal else 'correctness, ratio informational'})")
    print(f"report: {report_path}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
