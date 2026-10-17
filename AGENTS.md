# Notes for coding agents

Package tests: `julia --project -e 'using Pkg; Pkg.test()'`.

## Changing a numerical policy

A numerical policy is any solver decision rule: Newton damping, step bounds, norms, convergence tests, time stepping, Jacobian reuse, grid refinement, and the handling of negative or trace species. Most were adopted from Cantera, and a later batched GPU or reduced-precision path will need to change some of them. Matching a rule is not knowing why it exists.

When a change adds, removes or alters such a rule:

1. Read its row and any know-why record in [`docs/src/design_rationale.md`](docs/src/design_rationale.md).
2. Write the predicted effect (solver counts, solutions, which flame kinds) in the PR before running anything.
3. Compare baseline and candidate on the same cases. Use snapshot replay ([`test/cantera_parity/snapshots`](test/cantera_parity/snapshots/README.md)) for the local effect and the full adaptive solve for the total effect. Report `solver_stats` counts, grid sizes and accuracy per case.
4. In the PR, state the prediction, the outcome and any mismatch between them. A wrong prediction is useful; explain it.
5. Update the know-why record. *Unknown* is an allowed status; an invented rationale is not. "No effect on these cases" means unexercised, not unnecessary.
6. If the rule turns out load-bearing, add a regression test that fails when the rule is removed.

Pure refactors are a separate class: they must be bit-identical in states and solver counters, and they do not mix with policy changes in one commit.
