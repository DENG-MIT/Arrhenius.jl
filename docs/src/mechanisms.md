# Mechanisms and preprocessing

`CreateSolution("mechanism.yaml")` parses gas species, thermodynamics, and supported reaction data from YAML in Julia. The input path is composed from separate native kinetics, YAML thermo/species, and optional transport-fit readers. A `.yaml.npz` sidecar is optional and supplies pre-fitted transport properties when present. Cantera and Python are not required to load gas chemistry.

Native gas kinetics support elementary and three-body Arrhenius reactions, Lindemann and Troe falloff, pressure-dependent Arrhenius (PLOG) rates, linear-Burke/Lindemann mixture rules, Blowers–Masel rates, explicit reaction orders, and NASA7, NASA9, Shomate, and constant-cp thermodynamic models. SRI, Tsang, chemically-activated, and Chebyshev rate forms are not supported by the native reader.

YAML phase definitions may select local or imported species and reaction sections. Imported files are resolved relative to the mechanism file; optional data_paths can provide additional search directories:

~~~julia
gas = CreateSolution("mechanism.yaml"; data_paths=["shared-data"])
~~~

Transport sidecars are optional for chemistry but should be regenerated after changing the mechanism or an imported YAML file when fitted transport data are needed. For multicomponent diffusion or Soret diffusion, export matching collision-integral data with mechanism/export_multicomponent.py; see [flame transport](flames.md#transport). Use `kinetics_reader=:sidecar` only when comparing against legacy kinetics arrays.

Mechanism checks and an opt-in NASA7 continuity refit are available in Julia:

~~~julia
report = audit_mechanism("mechanism.yaml")
report.issues
refit_nasa7_continuity("mechanism.yaml", "mechanism-refit.yaml")
~~~

The audit checks balance, unused species, duplicate candidates, NASA polynomial continuity, negative pre-exponential declarations, collision-limit screens, and broad transport parameter ranges. Collision and transport thresholds are heuristic review flags. Duplicate groups are candidates and may include false positives because the screen does not fully reproduce Cantera's rate-type and third-body comparisons.

The refit adjusts each high-temperature NASA7 `cp/R` polynomial by constrained least squares so `cp/R` and its first derivative match at `Tmid`, then anchors `h/RT` and `s/R` there. It writes a separate YAML file and reports the largest changes to `cp/R`, `h/RT`, `s/R`, and `g/RT`. It removes `note` annotation fields in the generated copy because YAML.jl does not preserve multiline note indentation.

### Common-midpoint NASA7 refit

For a single shared NASA7 split temperature, use the opt-in refit:

```julia
report, output = refit_nasa7_common_tmid("mechanism.yaml", "mechanism-common-tmid.yaml")
```

It requires the chosen `Tmid` to lie inside the shared interval of all supported NASA7 species. It jointly fits both `cp/R` segments, enforces value and first-derivative continuity, preserves `h/RT` and `s/R` at each species' lower temperature bound, and anchors their continuity at the shared midpoint. The report gives maximum absolute changes in `cp/R`, `h/RT`, `s/R`, and `g/RT`. The output is separate; unsupported NASA7 layouts are listed as skipped. The fit is an opt-in least-squares transformation, so inspect its reported changes before using the generated mechanism.

Mechanism-audit context: Killingsworth et al. describe automated kinetic-mechanism analysis and debugging in [Combustion and Flame 221 (2020), 170–179](https://doi.org/10.1016/j.combustflame.2020.06.010); related Cantera enhancement discussion: [ChemCheck issue #42](https://github.com/Cantera/enhancements/issues/42).

```bibtex
@article{Killingsworth2020CloudMechanismAnalysis,
  author  = {Killingsworth, Nick J. and McNenly, Matthew J. and Whitesides, Russell A. and Wagnon, Scott W.},
  title   = {Cloud based tool for analysis of chemical kinetic mechanisms},
  journal = {Combustion and Flame},
  volume  = {221},
  pages   = {170--179},
  year    = {2020},
  doi     = {10.1016/j.combustflame.2020.06.010}
}
```