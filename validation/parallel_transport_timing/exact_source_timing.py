"""Time the unchanged source script, including its native spawn/import behavior.

Its four reported timers use three decimal places. Subtract 0.0005 s per timer
to use a conservative rounding bound. No callback, initializer or import wrapper
is present inside these child processes. Property correctness and actual worker
settings are checked separately by controlled_reference.py.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(exist_ok=False)
    mechanism = args.input / "gri30.yaml"
    hashes = {"source": sha(args.source), "mechanism": sha(mechanism), "driver": sha(__file__)}
    assert hashes["source"] == "6dd426c24a6af5a5b9d2e87d2d9b33e78de83645952b946e81d34e727debc1de"
    assert hashes["mechanism"] == "06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345"
    runs = []
    for label in ["cold", "warmup"] + [f"warm-{i:02d}" for i in range(9)]:
        result = subprocess.run([sys.executable, str(args.source.resolve())], cwd=args.input,
                                capture_output=True, text=True, timeout=60)
        (args.output / (label + ".stdout")).write_text(result.stdout)
        (args.output / (label + ".stderr")).write_text(result.stderr)
        result.check_returncode()
        matches = re.findall(r"^(Parallel|Serial): ([0-9.]+) seconds$", result.stdout, re.M)
        assert [name for name, _ in matches] == ["Parallel", "Serial", "Parallel", "Serial"]
        seconds = [float(value) for _, value in matches]
        row = dict(label=label, exit_code=result.returncode, reported_seconds=seconds,
                   conservative_seconds=[max(0, value-0.0005) for value in seconds])
        runs.append(row)
        print(label, sum(row["conservative_seconds"]), flush=True)
    assert hashes == {"source": sha(args.source), "mechanism": sha(mechanism), "driver": sha(__file__)}
    (args.output / "summary.json").write_text(json.dumps(dict(passed=True, runs=runs,
        hashes=hashes, command=[sys.executable, str(args.source)],
        child_environment={k: os.environ.get(k) for k in ("PYTHONPATH", "OPENBLAS_NUM_THREADS",
            "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}), indent=2) + "\n")


if __name__ == "__main__":
    main()
