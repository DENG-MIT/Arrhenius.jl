# Notes for coding agents

Package tests: `julia --project -e 'using Pkg; Pkg.test()'`.

## Changing a numerical policy

A numerical policy is any solver decision rule: Newton damping, step bounds, norms, convergence tests, time stepping, Jacobian reuse, grid refinement, and the handling of negative or trace species. Most were adopted from Cantera; see [`docs/src/design_rationale.md`](docs/src/design_rationale.md).

When a change adds, removes or alters such a rule:

1. Write the predicted effect (solver counts, solutions, which flame kinds) in the PR before running anything.
2. Compare baseline and candidate on the same cases. Use snapshot replay ([`test/cantera_parity/snapshots`](test/cantera_parity/snapshots/README.md)) for the local effect and the full adaptive solve for the total effect. Report `solver_stats` counts, grid sizes and accuracy per case.
3. In the PR, state the prediction, the outcome and any mismatch between them.
4. If removing the rule breaks a case, add a regression test that fails when the rule is removed.

Pure refactors are a separate class: they must be bit-identical in states and solver counters, and they do not mix with policy changes in one commit.
