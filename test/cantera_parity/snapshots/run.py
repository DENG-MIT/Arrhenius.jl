#!/usr/bin/env python3
"""Capture, replay steady stages to first divergence, and run xresid per stage."""
import csv, json, os, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = os.environ.get("RUN_CAPPED") or shutil.which("run_capped")
if not RUN:
    raise SystemExit("run_capped is required on PATH or via RUN_CAPPED")
JL = os.environ.get("JL") or shutil.which("julia")
if not JL:
    raise SystemExit("set JL to the Julia executable")
PROJECT = os.environ.get("ARRHENIUS_PROJECT", str(HERE.parents[2]))
CAP = os.environ.get("SNAPLIB_CAP_SECONDS", "1800")

def capped(argv):
    subprocess.run([RUN, CAP, *map(str, argv)], check=True)

def main(spec, outdir):
    spec = Path(spec).resolve(); outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    capped([sys.executable, HERE / "capture.py", spec, outdir])
    replay = HERE / "replay.jl"
    capped([JL, f"--project={PROJECT}", replay, spec, outdir, outdir / "replay.csv"])
    with open(outdir / "replay.csv", newline="") as stream:
        rows = list(csv.DictReader(line for line in stream if not line.startswith("first_divergence,")))
    manifest = json.loads((outdir / "manifest.json").read_text())
    xresid = []
    for row in rows:
        output_pos = int(row["output"]) - 1  # Julia's one-based stage position
        input_pos = int(row["input"]) - 1
        ct_stage = manifest["stages"][output_pos]
        input_stage = manifest["stages"][input_pos]
        native_state = outdir / f"native_stage_{row['input']}.npz"
        if input_pos == 0:
            # A same-file A/B check makes state transfer at the first replay
            # input explicit before interpreting any solver divergence.
            input_dir = outdir / f"xresid_input_{row['input']:0>4}"
            input_json = input_dir / "report.json"
            input_file = outdir / input_stage["file"]
            capped([sys.executable, HERE / "xresid_adapter.py", spec,
                    input_file, input_file, "--out", input_dir,
                    "--json", input_json, "--julia", JL])
            input_result = json.loads(input_json.read_text())
            xresid.append({"stage": row["input"], "role": "replay_input",
                           **input_result.get("aggregates", {})})
        xdir = outdir / f"xresid_{row['output']:0>4}"
        xjson = xdir / "report.json"
        capped([sys.executable, HERE / "xresid_adapter.py", spec,
                outdir / ct_stage["file"], native_state, "--out", xdir,
                "--json", xjson, "--julia", JL])
        result = json.loads(xjson.read_text())
        xresid.append({"stage": row["output"], "role": "replay_output",
                       **result.get("aggregates", {})})
    report = {"case": str(spec), "capture_stages": len(manifest["stages"]),
              "capture_seconds": manifest["total_seconds"], "replayed_stages": rows,
              "xresid": xresid}
    (outdir / "summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    if len(sys.argv) != 3: raise SystemExit("usage: run.py SPEC.json OUTDIR")
    main(sys.argv[1], sys.argv[2])
