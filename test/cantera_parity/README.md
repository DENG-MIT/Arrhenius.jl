# cantera_parity — lockstep native/Cantera comparison oracle

Dumps and compares states, residuals, Jacobians, and solver decisions between
Arrhenius.jl (this tree) and Cantera for a premixed flame case, at matching
states. Written for the stagni_2023 free phi1.2/h2o0.3/p1 lockstep task; works
for any `free`/`burner` spec of the nh3-parity harness.

## Layout contract

Cantera flat vector: per domain, point-major — `x[start + j*ncomp + c]`; flame
components are `velocity, spreadRate, T, Lambda, eField, Uo, Y_1..Y_n`
(unused free-flow components pinned to 0). Native state is a `(n+2) × N`
matrix: row 1 = T/1000, rows 2..n+1 = Y_k, row n+2 = mdot; flattened
column-major `(j-1)*(n+2) + k`. Native residuals carry the
`_flame_timescale = 1e-4` scaling (energy row an extra 1/1000 from kK), so
residual agreement with Cantera is "ratio constant per component", and Newton
steps (scale-invariant) are compared directly.

## Tools

- `ct_dump.py SPEC OUTDIR state` — Cantera initial-guess state (as the
  auto-solve path sets it, incl. the fixed-T grid insertion), steady residual
  per component (`Sim1D.eval(0)` + `Domain1D.residuals`), dense steady Jacobian
  (`eval_jacobian` via a temporary EigenSparseDirectJacobian container) →
  `ct_state.npz` + `ct_meta.json` + `ct_guess.csv`.
- `ct_dump.py SPEC OUTDIR trace` — full solve; stdout log at `CANTERA_LOGLEVEL`
  (default 2; >=4 switches to initial-grid-only, auto=false) → `ct_solve.log`;
  per-successful-timestep states via `set_time_step_callback` → `ct_trace.npz`.
- `native_dump.jl SPEC OUTDIR state CTDIR` — loads `CTDIR/ct_state.npz` into
  native flames (same grid; mdot = rho*u like `set_initial_guess!`), dumps
  residuals, the production Newton step (production `_flame_jacobian` + LAPACK
  banded solve), and the bound fraction for `:conservative` and
  `:finite_difference` → `native_state_cons.npz` / `native_state_fd.npz`.
- `native_dump.jl SPEC OUTDIR trace [scheme]` — full `solve!` with FLAMEDBG
  decision lines (requires `Arrhenius.DEBUG_FLAME[] = true`; hooks on branch
  nh3/32-lockstep print fbound/trials/timestep accept-fail/Jacobian evals) →
  stdout log + `native_trace[_fd].npz` with per-step states and cumulative
  counts.
- `compare_lock.py DIR resid|step|invariants` — Level 1 (per-component residual
  ratio constant-spread), Level 2 (Newton step per component; CT step rebuilt
  from dumped F and J), Level 3 (per-step invariant errors: element balance
  H/N/O, sum Y-1, mass-flux continuity; native-vs-CT state diff while grids
  match). `invariants` needs `LOCKSTEP_SPEC=<spec path>` and cwd with the
  mechanism yaml.
- `parse_decisions.py ct|native LOG` — prints the decision sequence (steady
  attempts, per-iteration alpha/fbound/s0/s1/N_jac, time steps, refines).

Run everything from the nh3-parity directory (mechanism yamls + sidecars), e.g.
`~/work/nh3-parity`, with env.sh sourced (PYTHONPATH to the Cantera 4.0 install
and PY/JL set). Native runs: `julia --project=<this worktree>`.

Reference results (stagni_2023 free phi1.2 h2o0.3 p1, 9-pt initial grid):
FD residual = 1e-4 x Cantera to roundoff; conservative differs structurally at
flame-front points; FD/Cantera first Newton fbound 3.23e-5/3.21e-5; first
decision divergence is the dt growth after time step 1.
