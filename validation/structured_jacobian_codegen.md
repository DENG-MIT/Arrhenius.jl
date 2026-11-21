# Code-generated structured reactor Jacobian

Follow-up to `validation/kinetics_scalar_scaling.md`'s investigation (reordering/kernel work,
scrapped: ~2-11% gains, not worth a second kinetics evaluator). This instead ports the
`LLMImprovments` fork's `@generated`-function Jacobian code generator (`src/Jacobian.jl` there)
onto this codebase's own existing analytical-Jacobian system.

## Context: an existing system, not a blank slate

`example/reactors/structured/` already has a substantial (1284-line), hand-written analytical
`d(wdot)/d[Y,T]` Jacobian for `IdealGasReactor` (`structured_jacobian.jl`), with its own
Cantera-compatible negative-concentration-safe mass-action/gradient math
(`signed_products.jl`/`signed_product_gradient.jl`), a per-call ambiguous-zero-crossing guard
falling back to finite differences (`guarded_jacobian.jl`), and a sparse-plus-rank-2-border
(KLU + Schur complement) linear solve wired into `OrdinaryDiffEqBDF.QNDF()`
(`klu_structured_solver.jl`/`qndf_structured_adapter.jl`, via `example/reactors/ode_solver.jl`'s
`native_bdf`). Scope: `:constant_pressure`/`:adiabatic`, Float64/NASA7 thermo, integer
stoichiometry-equals-order, no Blowers-Masel, no PLOG with concentration colliders.

There was no code generation anywhere in that system before this work -- everything was a
runtime loop over a precomputed `StructuredReactionPlan`/`SignedProductPlan` schedule. That is
architecturally the same tier as the fork's `AnalyticalIsobaricJacobianWorkspace`
(hand-written, schedule-based); the fork's `GeneratedAnalyticalIsobaricJacobianWorkspace` tier
(code-generated, same math, compiled per-mechanism into one flat function) had no counterpart
here. This work adds that missing tier.

## Design: generate only the hot inner step

`StructuredTrialJacobian` gained one new field, `fill_core!::Function` (defaulting to the
existing `_fill_concentration_core!`), and its callable now dispatches through it instead of
calling that function directly. Everything else -- the struct, `_check_mechanism!`, mass_core
scaling, `q`, `rank_vector`, the temperature column (still an AD-based single-direction
`ForwardDiff` dual pass, not code-generated), the energy row, and all of
`guarded_jacobian.jl`/`klu_structured_solver.jl`/`qndf_structured_adapter.jl` -- is completely
unchanged. `generated_structured_jacobian(reactor)` (`generated_jacobian.jl`, new file) builds
an ordinary `structured_trial_jacobian(reactor)` and replaces only `fill_core!` with a
mechanism-specialized, cached, `@generated`-function kernel.

Per reaction, at generation time (not runtime):
- `kind in (1,2,3)` (`SignedProductPlan`'s classification -- 1/2/3 total repeated factors,
  closed-form) forward/reverse mass-action products and their gradients are fully unrolled as
  literal `Expr`s (ported from the fork's `product_expr`/`product_derivative_expr`), with the
  Cantera-compatible "≥2 negative factors → zero" guard baked in as a literal boolean
  expression. Integer order exponents become literal `Int`s, so `Base.literal_pow` turns
  `C[i]^2` into plain multiplication at compile time -- the same effect the fork's
  `k_Pi_C_rectantOrder!` gets by hand-coding `order==1 ? Ci : order==2 ? Ci*Ci : ...`, obtained
  here for free.
- `kind == 4` (fractional orders, or >3 total repeated factors) reactions fall back, per
  reaction, to calling the existing `signed_product_gradient!` and `_ambiguous_anyn_zero` check
  directly -- still one unrolled block per reaction, just reusing the exact existing
  (already-correct) math for that reaction rather than re-deriving it symbolically.
- Three-body/falloff collider derivatives (`_collider_derivative`, reusing
  `src/FlameDerivatives.jl`'s `_flame_base_rate`/`_flame_falloff_derivative`) are likewise left
  as ordinary function calls, not re-derived symbolically -- they only apply to the much
  smaller three-body/falloff subset, so inlining them isn't where the benefit comes from, and
  re-deriving Troe/Lindemann derivatives symbolically would be a correctness risk for no
  measured gain.
- Final scatter into `concentration_core.nzval` uses literal slot indices
  (`StructuredReactionPlan.scatter_slots`, known at generation time), not a runtime
  dependency-index loop.

## Correctness

Compared against the existing `structured_trial_jacobian` (bit-identical required, not just
close) on `h2o2`, `gri30`, the PLOG-bearing (collider-free) `stagni_2023`, and
`nc7_ver3.1_mech` (4846 reactions) -- including repeated random states with a forced slightly
negative trial concentration in a random species each time, to exercise the negative-
concentration guards:

**Every comparison was exactly bit-identical (`max rel diff = 0.0`), including under forced
negative trial concentrations.** Full end-to-end ignition solves (`native_bdf`,
`linear_solver=:generated` vs `:dense`) also matched bit-for-bit in final state, `njacs`, and
total step count on both `h2o2` and `gri30`.

The full package test suite (`test/runtests.jl`) and `test/structured_reactors.jl` (not wired
into `runtests.jl`; has one pre-existing, unrelated failing assertion --
`signed RHS and structured columns`'s `clipped != signed_negative` check, confirmed to fail
identically on the unmodified file before this work) both pass unaffected.

## Performance

**Per-Jacobian-call** (median of 200 samples, 30 for `nc7_ver3.1_mech`, `JULIA_NUM_THREADS=1`):

| Mechanism | reactions | analytic (ns) | generated (ns) | speedup |
|---|---:|---:|---:|---:|
| h2o2 | 28 | 5375 | 4333 | 1.24x |
| stagni_2023 | 203 | 34375 | 24959 | 1.38x |
| gri30 | 325 | 60750-63291 | 46666-47000 | 1.30-1.35x |
| nc7_ver3.1_mech | 4846 | 1,189,646 | 907,751 | 1.31x |

Real, consistent, reproducible across three repeated runs on `gri30` and across four mechanism
sizes spanning 28 to 4846 reactions -- unlike the kinetics-reordering work, this isn't
workload-dependent (it helps on every mechanism, not just low-reversibility ones).

**End-to-end ignition solve** (`native_bdf`, `linear_solver=:dense` vs `:generated`,
`reltol=1e-10`, `abstol=1e-14`): the per-call win above **does not translate into a comparable
end-to-end wall-clock improvement**:

| Mechanism | `destats` (nf, nw, nsolve, njacs, ..., naccept) | dense (s) | generated (s) | end-to-end speedup |
|---|---|---:|---:|---:|
| gri30 | (7337, 123, 5687, **8**, ..., 1647) | 0.0832 | 0.0826 | ~1.01x |
| nc7_ver3.1_mech | (2767, 123, 2079, **29**, ..., 685) | 0.831 | 0.824 | ~1.01x |

Only 8-29 Jacobian evaluations occur across 685-1647 accepted steps -- QNDF/BDF reuses a
factored Jacobian across many steps via its own γ-scaling, so rebuilding it is a small fraction
of total cost. The dominant costs (RHS evaluations: 2767-7337; linear solves: 2079-5687) are
untouched by this work, since `fill_core!`'s replacement only changes how the Jacobian itself
is computed, not the RHS or the linear solve.

## Net assessment

The code generator delivers a real, consistent, mechanism-size-independent ~1.2-1.4x speedup
on Jacobian evaluation itself, fully verified correct including negative-concentration edge
cases and a 4846-reaction mechanism. It is not, on its own, a comparable end-to-end ignition
speedup for this solver configuration (`QNDF` with dense linear algebra), because Jacobian
evaluation is not the dominant cost once a BDF-family solver is reusing factorizations across
steps -- unlike the fork's own report (`gri30 generated ≈107ms vs analytical ≈248ms` end to
end), which compared *supplying an analytic Jacobian at all* against *no Jacobian* (forcing
many more implicit finite-difference/dual evaluations), not two analytic-Jacobian
implementations against each other as this comparison does. Where this would matter more
directly: workloads that evaluate the Jacobian far more frequently relative to the RHS --
Newton-per-iteration contexts (e.g. a steady flame solve) rather than a BDF time-stepper that
reuses it, or direct repeated Jacobian queries (parameter studies, sensitivity/optimization
loops) -- untested here, and a natural next place to look if the ~1.3x per-call number needs
to show up end to end.
