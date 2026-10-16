# Solver options

## Spatial discretization

FreeFlame and BurnerFlame accept :finite_difference and :conservative. The constructor default is :finite_difference. Set discretization=:conservative to use finite-volume species and total-enthalpy balances. Ionized flames also accept both options; their initial conservative solve first uses a frozen finite-difference stage.

## Grid refinement and solve controls

solve! adaptively refines a grid by default. Common controls are:

- ratio, slope, and curve: refinement thresholds.
- prune: remove points whose indicators fall below the threshold; 0 disables pruning.
- max_points and max_time_steps: bounds for the solve.
- initial_time_step, time_step_growth_factor, and time_step_growth_strategy: pseudo-time stepping controls.
- time_step_callback and time_step_regrid: observe or recover from stalled time stepping.
- jacobian_mode=:auto or :finite_difference: use the analytic or numerical Jacobian path where available.
- auto=true: for FreeFlame, use a staged solve and domain expansion; counterflow solvers use their own recovery sequence.

Refinement settings and supported controls vary by flame family. Consult the corresponding constructor and solve! docstrings in the [API reference](api.md).

## Transport

set_transport!(flame, :mixture_averaged) selects mixture-averaged diffusion. :multicomponent requires MultiTransportData exported from the same mechanism. Set soret=true to include thermal diffusion; this also requires matching multicomponent data. For mixture-averaged transport, flux_gradient_basis=:mole or :mass chooses the gradient basis.

Counterflow diffusion flames can enable optically thin radiation with set_radiation!. Radiation and transport options apply only to flame families that expose them.
