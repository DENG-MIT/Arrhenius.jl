# Validation methods

Arrhenius.jl validation uses two distinct comparisons.

## Same-discrete-problem comparison

To compare solver implementations, both evaluate the same discrete problem: the same grid coordinates, physical state, boundary conditions, mechanism, transport model, tolerances, and solver controls. Residuals are compared in physical units at that shared state. This checks implementation parity for a declared configuration; it does not establish accuracy on a different grid.

## Accuracy against a converged reference

To assess solution accuracy, a flame solution is compared with a separately converged reference. Grid-refinement studies estimate how the result changes as spatial resolution increases. Solver convergence and grid accuracy are checked separately; a small change between two grids alone is not treated as proof of convergence.

Validation cases must represent burning flames when the question concerns a burning solution. Extinguished solutions are reported as extinguished and are not counted as successful burning cases.

These methods describe the validation approach. They do not assert a particular accuracy level for every mechanism, flame family, or operating condition.
