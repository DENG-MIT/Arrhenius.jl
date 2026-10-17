#!/usr/bin/env python3
"""Compare same-case Arrhenius solver archives for solution neutrality."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

STATE_KEYS = (
    "grid", "T", "Y", "velocity", "inlet_Y", "SL",
    "newton_ok", "native_residual_inf_scaled", "species_names_utf8",
)


def compare(base_dir: Path, candidate_dir: Path, tolerance: float, ignored_counter_keys=()):
    base_files = {p.name: p for p in base_dir.glob("*.npz")}
    candidate_files = {p.name: p for p in candidate_dir.glob("*.npz")}
    missing = sorted(set(base_files) ^ set(candidate_files))
    cases = []
    for name in sorted(set(base_files) & set(candidate_files)):
        with np.load(base_files[name], allow_pickle=False) as base, np.load(candidate_files[name], allow_pickle=False) as cand:
            state_keys = list(STATE_KEYS)
            counter_keys = sorted(key for key in set(base.files) | set(cand.files) if key.startswith("stats_"))
            absent = [key for key in state_keys + counter_keys if key not in base or key not in cand]
            errors = {}
            identical = True
            case_ok = not absent
            for key in state_keys:
                if key in absent:
                    continue
                a, b = np.asarray(base[key]), np.asarray(cand[key])
                if a.shape != b.shape:
                    errors[key] = {"shape_base": list(a.shape), "shape_candidate": list(b.shape)}
                    case_ok = False
                    identical = False
                    continue
                same = np.array_equal(a, b)
                identical &= same
                if key == "species_names_utf8":
                    if not same:
                        errors[key] = {"identical": False}
                        case_ok = False
                    continue
                delta = np.abs(a.astype(np.float64) - b.astype(np.float64))
                scale = max(float(np.max(np.abs(a))) if a.size else 0.0, np.finfo(float).tiny)
                rel = float(np.max(delta) / scale) if delta.size else 0.0
                if rel > tolerance:
                    case_ok = False
                if rel:
                    errors[key] = {"max_abs": float(np.max(delta)), "max_scaled_relative": rel}
            counter_mismatches = {}
            ignored_counter_mismatches = {}
            for key in counter_keys:
                if key in absent:
                    continue
                if not np.array_equal(base[key], cand[key]):
                    target = ignored_counter_mismatches if key in ignored_counter_keys else counter_mismatches
                    target[key] = {
                        "base": np.asarray(base[key]).tolist(),
                        "candidate": np.asarray(cand[key]).tolist(),
                    }
                    if key not in ignored_counter_keys:
                        case_ok = False
            y = np.asarray(cand["Y"], dtype=np.float64)
            species_count = len(bytes(np.asarray(cand["species_names_utf8"], dtype=np.uint8)).decode("utf-8").splitlines())
            species_axis = 0 if y.shape[0] == species_count else -1
            sums = y.sum(axis=species_axis)
            sum_y_error = float(np.max(np.abs(sums - 1.0))) if sums.size else 0.0
            cases.append({
                "id": name.removesuffix(".npz"), "pass": case_ok, "state_bit_identical": identical,
                "max_sumY_error": sum_y_error, "state_errors": errors,
                "counter_mismatches": counter_mismatches,
                "ignored_counter_mismatches": ignored_counter_mismatches, "missing_arrays": absent,
            })
    result = {
        "base_dir": str(base_dir), "candidate_dir": str(candidate_dir), "tolerance": tolerance,
        "case_count": len(cases), "missing_cases": missing, "ignored_counter_keys": sorted(ignored_counter_keys),
        "passed": not missing and bool(cases) and all(c["pass"] for c in cases),
        "bit_identical_cases": sum(c["state_bit_identical"] and not c["counter_mismatches"] for c in cases),
        "max_sumY_error": max((c["max_sumY_error"] for c in cases), default=None), "cases": cases,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--ignore-counter", action="append", default=[],
                        help="report but do not fail on this counter key; may be repeated")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compare(args.base, args.candidate, args.tolerance, args.ignore_counter)
    rendered = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())


