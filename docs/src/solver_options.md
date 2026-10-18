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

## Cantera-compatible solve switch

`cantera_compatible=true` selects the following settings without changing other
solve keywords. The resolved settings from the latest call are available through
`solver_options(flame)`.

The resolved Newton preset for every listed kind is:

| `NewtonPolicyOptions` field | Value |
| --- | --- |
| `cantera_compatible` | `true` |
| `damping_trials`, `damping_factor`, `damping_acceptance` | `7`, `sqrt(2)`, `:step_or_reduction` |
| `step_bounds`, `norm_denominator`, `flow_weight`, `steady_convergence` | `:cantera`, `:cantera`, `:velocity`, `:weighted_step` |
| `premixed_steady_maxiters`, `counterflow_steady_maxiters`, `premixed_counterflow_steady_maxiters` | `0` (uncapped) |
| `premixed_transient_maxiters`, `premixed_policy_transient_maxiters`, `counterflow_transient_maxiters` | `0` (uncapped) |
| `jacobian_age`, `jacobian_reevaluations` | `20`, `4` |

| Flame kind | Jacobian mode | Staged ladder selected | Existing `auto` default | Solve path |
| --- | --- | --- | --- | --- |
| `FreeFlame` | `:auto` | Yes | `false` | Premixed staged ladder (uniform 12/24/48-point starts, locked-temperature retry, refinement; free-flame domain expansion when needed) |
| `BurnerFlame` (energy-on) | `:auto` | Yes | `false` | Premixed staged ladder without domain expansion |
| `BurnerFlame` (fixed-T profile) | `:auto` | No | `false` | Cantera `auto=False` refinement loop with compatible Newton/Jacobian settings |
| `CounterflowDiffusionFlame` | `:auto` | No | `true` | Existing counterflow recovery sequence |
| `CounterflowPremixedFlame`, `CounterflowTwinPremixedFlame`, `ImpingingJet` | `:auto` | No | `true` | Existing premixed-counterflow recovery sequence |

For every listed kind, `cantera_compatible=false` retains the current Jacobian,
Newton, and solve-path defaults.

## Transport

set_transport!(flame, :mixture_averaged) selects mixture-averaged diffusion. :multicomponent requires MultiTransportData exported from the same mechanism. Set soret=true to include thermal diffusion; this also requires matching multicomponent data. For mixture-averaged transport, flux_gradient_basis=:mole or :mass chooses the gradient basis.

Counterflow diffusion flames can enable optically thin radiation with set_radiation!. Radiation and transport options apply only to flame families that expose them.
