#!/usr/bin/env python3
"""Compare default Newton-policy outputs against a baseline run.

Timing arrays are excluded; all state, scalar, and solver-stat arrays must match exactly.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


def compare_case(case_id, baseline_path, candidate_path):
    if not baseline_path.is_file() or not candidate_path.is_file():
        return {"case": case_id, "equal": False, "missing": [
            str(p) for p in (baseline_path, candidate_path) if not p.is_file()]}
    with np.load(baseline_path) as base, np.load(candidate_path) as candidate:
        keys = (set(base.files) | set(candidate.files)) - {"batch_means"}
        differences = []
        for key in sorted(keys):
            if key not in base or key not in candidate:
                differences.append({"key": key, "reason": "missing array"})
                continue
            left, right = base[key], candidate[key]
            if left.shape != right.shape:
                differences.append({"key": key, "reason": "shape mismatch",
                                    "baseline_shape": left.shape, "candidate_shape": right.shape})
            elif not np.array_equal(left, right):
                delta = np.abs(left.astype(float) - right.astype(float))
                differences.append({"key": key, "reason": "values differ",
                                    "max_abs_difference": float(delta.max(initial=0.0))})
    return {"case": case_id, "equal": not differences, "differences": differences}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases_csv", type=Path)
    parser.add_argument("baseline_dir", type=Path)
    parser.add_argument("candidate_dir", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()
    with args.cases_csv.open(newline="") as stream:
        cases = list(csv.DictReader(stream))
    if not cases or not {"id"} <= set(cases[0]):
        parser.error("cases CSV must contain an id column and at least one row")
    results = [compare_case(row["id"], args.baseline_dir / f"{row['id']}.npz",
                            args.candidate_dir / f"{row['id']}.npz") for row in cases]
    report = {"case_count": len(results), "exact_match_count": sum(r["equal"] for r in results),
              "all_exact": all(r["equal"] for r in results), "cases": results}
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    print(f"exact_matches={report['exact_match_count']}/{report['case_count']} output={args.output_json}")


if __name__ == "__main__":
    main()
