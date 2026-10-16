# Cantera stage snapshots and native replay

This directory captures and replays representative Cantera flame-solver stages for free, burner, and axisymmetric diffusion flames. Successful time-step and steady callbacks write compressed NPZ states for every domain, grid, event type, timestep, and available solver statistics.

A second, fixed-grid diagnostic run starts from the same initial guess and records the first transient batch from Cantera's log, including rejected timestep attempts, applied damping, Jacobian refresh counts, and accept/reject outcomes. The report compares the native solver against Cantera's first attempted timestep and stops at the first differing decision or first failed stage.

The diagnostic run saves its own initial state and checks every serialized array against capture stage 1; `ct_decisions.json` reports exact equality and maximum absolute difference before decision counts are compared. The level-2 Newton parser combines refresh markers with the table's cumulative `N_jac` column, so each refresh is attributed to the attempt that triggers it while the cumulative count remains available separately.

## Run

Set `NH3_HARNESS_DIR` to the Cantera harness checkout, `ARRHENIUS_PROJECT` to this repository checkout, and `XRESID_PATH` to the xresid checkout. The harness environment supplies `PY`, `JL`, and `run_capped`:

```bash
: "${NH3_HARNESS_DIR:?set the Cantera harness directory}"
: "${ARRHENIUS_PROJECT:?set the Arrhenius checkout}"
: "${XRESID_PATH:?set the xresid checkout}"
source "$NH3_HARNESS_DIR/env.sh"
cd "$NH3_HARNESS_DIR"
run_capped 1800 "$PY" \
  "$ARRHENIUS_PROJECT/test/cantera_parity/snapshots/run.py" \
  specs_free_ma/alzueta_2023__free__phi0.8_h20_p1__mixture-averaged.json \
  /tmp/free-snapshot-replay
```

`run.py` starts each solver command through `run_capped` (default child cap 1800 seconds; set `SNAPLIB_CAP_SECONDS` lower to reduce it). Outputs include `manifest.json`, `stage_*.npz`, `ct_decisions.log`, `ct_decisions.json`, `replay.csv`, native replay snapshots/logs, per-stage xresid directories, and `summary.json`.

The replay uses the captured input state and grid. At the initial transient stage it uses Cantera's first attempted `dt` (which can be rejected); later successful time-step stages use the callback `dt`. A steady stage invokes one native steady Newton solve. If grids differ, the native refiner runs with the case's ratio/slope/curve/prune criteria and its resulting grid is compared. State tolerance is componentwise relative `1e-8`; grid coordinates use absolute `1e-12 m`.

For the first transient attempt, the report compares accept/reject outcome, applied damping, and Jacobian refresh count. It reports the first differing decision. The state comparison remains against Cantera's successful output snapshot, so a rejected replayed attempt is expected to differ from that output. The xresid cross-residual harness checks state-transfer round trips and evaluates native and Cantera residuals in physical units at both the captured output and native replay state.

Cantera callbacks do not expose every failed Newton trial or the refiner's internal point marks. The rejected-attempt decision trace comes from a separate log-level-4 fixed-grid diagnostic run initialized from the same initial guess, not from the adaptive capture. The tool currently reconstructs the first transient batch only; later rejected attempts and exact refinement point marks are not replayed. Grid agreement is measured from resulting coordinates.


## Parser regression test

The compact Nakamura diffusion fixture records three rejected attempts. The parser reports one Jacobian evaluation for each attempt and preserves Cantera's cumulative counts `[1, 2, 3]` separately. Run the focused test with:

```bash
: "${NH3_HARNESS_DIR:?set the Cantera harness directory}"
: "${ARRHENIUS_PROJECT:?set the Arrhenius checkout}"
source "$NH3_HARNESS_DIR/env.sh"
cd "$ARRHENIUS_PROJECT"
run_capped 60 "$PY" test/cantera_parity/snapshots/test_jacobian_trace.py
run_capped 60 "$JL" --project=. test/cantera_parity/snapshots/test_jacobian_trace.jl
```

This parser test is standalone and does not add solver or Cantera dependencies to `Pkg.test()`.
