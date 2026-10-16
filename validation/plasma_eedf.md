# Standalone electron-energy distribution

The [Julia example](../example/thermodynamics/plasma_eedf.jl) reproduces the
[Cantera air EEDF calculation](https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/thermo/plasma-eedf.py)
at 300 K, 101325 Pa and 200 Td using 43 collision tables and 40 energy cells.
Use the [Phelps air model](https://github.com/Cantera/cantera-example-data/blob/622e02535961079bf4e55dd338294e244a31cffb/air-plasma-Phelps.yaml)
from the pinned example-data revision. The dataset identifies its original
Phelps/LXCat sources and their citation instructions.

[WSL results](results/cantera4_wsl_plasma_eedf.json) and
[Apple M4 results](results/cantera4_m4_plasma_eedf.json) contain all timing samples,
output arrays, accuracy thresholds, and source/build hashes. Warm complete-call
speed ratios are 1.387 and 0.981 respectively (Cantera time / Julia time).
Every call reads the model and completes the solve; the first invocation is
reported separately. These results apply to this standalone fixed-grid case.

The reference fixture [plasma_eedf.npz](reference/plasma_eedf.npz) contains the
independently computed Cantera distribution and the BOLOS comparison points
published in the cited Cantera example. It is used only after native calculations.
The source's coarse-grid distribution differs from linearly interpolated BOLOS
values by about 5% in weighted L1; agreement with Cantera measures implementation
consistency. Mean energies are integrals of the returned distribution.

## Reproduce the measurements

Use Julia 1.12.7 and Cantera 4.0.0a2 built from the revision linked above.
Set numerical thread counts to one before starting either runtime:

```bash
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
julia --project=. validation/plasma_eedf_timing.jl "$PWD" /path/to/air-plasma-Phelps.yaml validation/reference/plasma_eedf.npz native-eedf-results
```

For the Python driver, create a JSON configuration with `model` (the model
path), `expected_cantera_version` (`4.0.0a2`), `expected_model_sha256`
(`be36740e0db01d996668b2c3ff6cbf12a0f296155086a16b7335d26df43f3ce2`),
and `expected_extension_sha256` (the SHA-256 of the pinned build's
`cantera._cantera.__file__`). Then run:

```bash
python validation/plasma_eedf_timing.py eedf-config.json cantera-eedf-results
```

Both output directories must be new. The drivers save one first call and nine
warm calls, retaining every sample. Compare every saved distribution before
dividing the warm medians; the drivers do not independently declare performance
qualification. Plotting and output serialization are outside the timers.


## Phase-selected collisions and continuation

`read_eedf_model(path; phase, data_paths)` resolves selected reaction sections
and imported species or reactions. Root collision tables remain included.
`EEDFState` accepts an explicit total `number_density` in m⁻³ for gases whose
density differs from `P/(kB*T)`.

Pass an earlier result with `solve_eedf(model, state; initial=previous)` to
continue from its center distribution on the same grid. The input remains
unchanged. At or below the low-field threshold, the solver resets to the
current gas-temperature Maxwellian, including valid underflowed zero tails.

[Pulse EEDF checks](results/cantera4_wsl_pulse_eedf_primitives.json) cover all
34 collision tables and the original methane-discharge example's initial and
190 Td states. [Full mechanism thermochemistry](plasma_thermochemistry.md) is validated separately;
the complete pulse trajectory remains unvalidated.
