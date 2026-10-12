# Complete parallel transport timing

These drivers compare all four 5,000-temperature calculations in Cantera's
`transport/multiprocessing_viscosity.py`: parallel and serial multicomponent
conductivity, followed by parallel and serial viscosity. Each call constructs
its phases and workspaces. Python pool creation and teardown are included.
Use the pinned Cantera 4.0 development source and the inputs prepared by
`validation/parallel_transport_cases.py`.

Run on an otherwise idle host, with four Julia threads, four Python workers
and single-thread numerical libraries. Set `OPENBLAS_NUM_THREADS`,
`OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`,
`NUMEXPR_NUM_THREADS` and `BLIS_NUM_THREADS` to `1` before starting Python.
The drivers also query the actual loaded numerical libraries, including
Accelerate on macOS. Reuse a Julia environment containing Arrhenius.

```sh
julia --threads=4 --project=JULIA_ENV validation/parallel_transport_timing/controlled_native.jl SOURCE_ROOT INPUT_DIR OUTPUT_DIR/native.npz
python validation/parallel_transport_timing/controlled_reference.py --source-root SOURCE_ROOT --source-example CANTERA_SOURCE_EXAMPLE --mechanism INPUT_DIR/gri30.yaml --output-dir OUTPUT_DIR/reference
python validation/parallel_transport_timing/exact_source_timing.py --source CANTERA_SOURCE_EXAMPLE --input INPUT_DIR --output OUTPUT_DIR/exact-source
python validation/parallel_transport_timing/summarize_controlled.py --folder OUTPUT_DIR --core-root ARRHENIUS_PACKAGE_ROOT/src --driver-root SOURCE_ROOT/validation/parallel_transport_timing --prior-correctness CORRECTNESS_REPORT
```

Create `OUTPUT_DIR` first; the reference drivers require fresh subdirectories.
`CORRECTNESS_REPORT` is the same host's full property-validation report.
`SOURCE_ROOT` contains the example and validation drivers; the loaded Arrhenius
core must match `ARRHENIUS_PACKAGE_ROOT`. The summary checks their actual bytes.

After one first call and one full warmup, nine complete calls are measured.
Every repeat is checked against the initial output. The worker audit observes
four distinct initialized workers before and after each original property map.
It subtracts inspection windows while retaining phase initialization through
the last original initializer's return timestamp. Original initializer and
property functions are called unchanged.

Two additional reference measurements guard against audit overhead: calls
without inspection, and separate executions of the unchanged original script.
Its printed times are reduced by half their rounding unit. For each repetition,
the reference is the smallest of those two totals and the adjusted audited
total. The reported speed ratio is the reference median divided by the Julia
median. The extra uninstrumented runs alone do not verify worker settings;
qualification also requires the audited workers and all numerical checks.

Imports, plotting, result-file writes and comparisons are outside timers.
The Julia first call includes specialization and is reported separately.
Worker inspection rejection checks can be run without Cantera:

```sh
python validation/parallel_transport_timing/controlled_guard_checks.py
```
