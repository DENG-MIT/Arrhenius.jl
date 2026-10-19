# Cantera stage snapshots and native replay

Capture Cantera flame solver stages and replay their transitions in Arrhenius.jl. The committed h2o2 specs use `mechanism/h2o2.yaml` from this repository, so the scripts work from a clean checkout without a harness directory. Cantera is needed for capture and xresid comparisons; Python NumPy and the Julia project dependencies are also required.

## Examples

From the repository root, with `run_capped`, `julia`, and a Python environment containing Cantera and NumPy available on `PATH`:

```sh
run_capped 1800 python3 test/cantera_parity/snapshots/run.py \
  test/cantera_parity/snapshots/specs/h2o2_free.json /tmp/h2o2-free 1
run_capped 1800 python3 test/cantera_parity/snapshots/capture.py \
  test/cantera_parity/snapshots/specs/h2o2_counterflow.json /tmp/h2o2-counterflow
run_capped 1800 julia --project=. test/cantera_parity/snapshots/replay.jl \
  test/cantera_parity/snapshots/specs/h2o2_counterflow.json /tmp/h2o2-counterflow /tmp/replay.csv 1
```

`run.py` accepts an optional maximum transition count for bounded diagnostics. It captures Cantera stages, replays them on matching native grids, and calls the in-repository `xresid.py` to check state transfer and residuals. `capture.py` can be run directly as `capture.py SPEC.json OUTDIR`; `replay.jl` accepts `SPEC.json CAPTURE_DIR REPORT.csv [MAX_TRANSITIONS]`.

The free flame is the executable smoke example. The counterflow spec documents a runnable h2o2 diffusion-flame setup. The replay reports per-transition state agreement and first divergence; xresid reports transfer round-trip and cross-residual agreement.

## Parser regression

The compact Nakamura fixture exercises the Cantera Jacobian-attempt parser without solver execution:

```sh
run_capped 60 python3 test/cantera_parity/snapshots/test_jacobian_trace.py
run_capped 60 julia --project=. test/cantera_parity/snapshots/test_jacobian_trace.jl
```
