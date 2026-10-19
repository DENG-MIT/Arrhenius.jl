#!/usr/bin/env python3
"""Capture BurnerFlame startup arrays through the production h_v2 build path."""
import os
import sys
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "1"


def run_case(adapter_path, spec_path, out):
    adapter_path = Path(adapter_path)
    source = adapter_path.read_text(encoding="utf-8")
    prefix, marker, _ = source.partition('if "reference" in opts:')
    if not marker:
        raise RuntimeError("h_v2 adapter entry marker not found")

    saved_argv = sys.argv
    sys.argv = [str(adapter_path), spec_path, out]
    namespace = {"__file__": str(adapter_path), "__name__": "h_v2_initial_capture"}
    try:
        exec(compile(prefix, str(adapter_path), "exec"), namespace)
    finally:
        sys.argv = saved_argv

    spec = namespace["spec"]
    if spec["kind"] not in ("fixed", "burner"):
        raise ValueError("unsupported case kind " + spec["kind"])
    flame = namespace["build"]()
    flame.set_initial_guess()
    np = namespace["np"]
    gas = namespace["gas"]
    np.savez_compressed(out, grid=np.asarray(flame.grid), T=np.asarray(flame.T),
        Y=np.asarray(flame.Y), velocity=np.asarray(flame.velocity),
        species_names_utf8=np.frombuffer("\n".join(gas.species_names).encode("utf-8"),
            dtype=np.uint8))
    print("CT_INITIAL case={} kind={} points={} adapter={}".format(
        spec["id"], spec["kind"], len(flame.grid), adapter_path), flush=True)


args = sys.argv[1:]
if len(args) < 3 or (len(args) - 1) % 2:
    raise SystemExit("usage: fixed_t_initial_ct.py h_v2_adapter.py spec.json output.npz [spec.json output.npz ...]")
for i in range(1, len(args), 2):
    run_case(args[0], args[i], args[i + 1])
