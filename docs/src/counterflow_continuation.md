# Counterflow flame continuation

Counterflow flames support two fixed-grid continuation drivers. For diffusion flames, `two_point_control!` follows Cantera's flame-control method by lowering the temperatures at two interior flame locations while solving for both inlet mass fluxes. It returns each control temperature, decrement, peak temperature, peak strain rate, inlet flux pair, state, and solve time. The default control spacing and temperature-decrement adaptation follow Cantera's `diffusion_flame_continuation.py` example. Pass `temperature_decrement_schedule` to replay a recorded schedule, or set `initial_time_step` to tune Arrhenius' transient startup for difficult correction steps. Grid refinement is disabled during the trace. This API is specific to `CounterflowDiffusionFlame`; pseudo-arclength remains available for diffusion and premixed counterflow flames. See `example/flames/counterflow_two_point_control.jl` for a runnable example.

Counterflow flames support opt-in pseudo-arclength continuation on a fixed grid:

- `pseudo_arclength!(flame; parameter=:inlet_flux_scale, steps=40, ds=0.05, parameter_scale=nothing, state_scale=1.0)` scales both inlet mass fluxes while preserving their ratio. `:fuel_mass_flux`, `:oxidizer_mass_flux`, `:pressure`, and (for an enabled diffusion flame) `:control_temperature` are also supported. The control-temperature parameter shifts both fixed-grid two-point temperatures together, so arclength can continue on the Cantera control family. For an enabled diffusion control family, `:control_temperature` shifts both control temperatures together; `direction=1` or `-1` selects the initial tangent orientation at a shared state.
- The corrector retains the banded flame Jacobian and solves the scalar parameter/arclength border with a Schur complement. A residual infinity norm is returned at every point; a `turning_point` flag is set when the parameter tangent changes sign.
- The initial flame must be converged. Continuation keeps the grid fixed; repeat on successively refined grids to estimate grid sensitivity. The solver adapts the arc step from corrector iteration count and halves it after a failed correction. `parameter_scale` defaults to the magnitude of the initial parameter (at least 1), and `state_scale` defaults to 1; set them to adjust the relative arc-length metric for parameter and flame-state motion.
- The API covers diffusion, opposed premixed, and twin premixed counterflow flames. It does not support catalytic impinging jets.

Use save_continuation_checkpoint(path, flame, point; parameter) to save the grid, state, parameter, physical tangent, arclength, next step, and metric scales. Resume with pseudo_arclength!(flame; restart=path, steps=...); the loader verifies mechanism, flame kind, species order, and boundary conditions. The H2/O2 example is in `example/flames/counterflow_pseudo_arclength.jl`. The upstream Cantera example uses two-point temperature control to march a stable/unstable branch, adjusting its two control temperatures and recording maximum temperature and strain rate. Cantera's selected control points are tied to grid crossings and its example allows grid refinement during continuation; those choices differ from this fixed-grid pseudo-arclength API.

## Literature

Nishioka, M., Law, C. K., and Takeno, T. (1996). “A flame-controlling continuation method for generating S-curve responses with detailed chemistry.” *Combustion and Flame*, 104(3), 328–342. https://doi.org/10.1016/0010-2180(95)00132-8. The paper reports one- and two-point temperature control and hydrogen-radical control for counterflow flames.

Giovangigli, V. and Smooke, M. D. (1987). “Calculation of extinction limits for premixed laminar flames in a stagnation point flow.” *Journal of Computational Physics*, 68(2), 327–345. https://doi.org/10.1016/0021-9991(87)90061-1.

Allgower, E. L. and Georg, K. (2003). *Introduction to Numerical Continuation Methods*. SIAM. https://doi.org/10.1137/1.9780898719154.

Quadarella, E., Guo, J., Cuoci, A., and Im, H. G. (2022). “Flame-controlling continuation method for extinction of counterflow sooting flames with detailed chemistry.” *AIAA SciTech 2022 Forum*. https://doi.org/10.2514/6.2022-2252. The authors cross-check against OpenSMOKE++.

BibTeX records for these sources are in `tasks/contin/continuation_refs.bib`. Keller's pseudo-arclength attribution was not independently verified and is not used as a source claim.
## BibTeX

```bibtex
@article{Nishioka_1996,
  title={A flame-controlling continuation method for generating S-curve responses with detailed chemistry},
  author={Nishioka, M. and Law, C. K. and Takeno, T.},
  journal={Combustion and Flame}, volume={104}, number={3}, pages={328--342}, year={1996},
  doi={10.1016/0010-2180(95)00132-8}
}
@article{Giovangigli_1987,
  title={Calculation of extinction limits for premixed laminar flames in a stagnation point flow},
  author={Giovangigli, V. and Smooke, M. D.},
  journal={Journal of Computational Physics}, volume={68}, number={2}, pages={327--345}, year={1987},
  doi={10.1016/0021-9991(87)90061-1}
}
@book{Allgower_2003,
  title={Introduction to Numerical Continuation Methods}, author={Allgower, Eugene L. and Georg, Kurt},
  publisher={Society for Industrial and Applied Mathematics}, year={2003},
  doi={10.1137/1.9780898719154}
}
@inproceedings{Quadarella_2022,
  title={Flame-controlling continuation method for extinction of counterflow sooting flames with detailed chemistry},
  author={Quadarella, Erica and Guo, Junjun and Cuoci, Alberto and Im, Hong G.},
  booktitle={AIAA SciTech 2022 Forum}, year={2022}, doi={10.2514/6.2022-2252}
}
```
