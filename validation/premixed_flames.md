# Premixed-flame validation

The complete free, burner and prescribed-temperature examples were checked with
Julia 1.12.7 and Cantera 4.0.0a2 at revision
[726522be](https://github.com/Cantera/cantera/tree/726522be4e2a13454d8415b7ef799d621f665cf3).
[Results and all timing samples](results/cantera4_flames_10dad59.json) identify
native revision 10dad59 and the source, mechanism and reference hashes.

| Complete source calculation | WSL speed ratio | Conservative ratio | Speed qualification |
| --- | ---: | ---: | --- |
| Free flame, four transport stages | 1.147 | 1.099 | Passed |
| Burner flame, two transport stages | 1.980 | 1.731 | Inconclusive: timing spread exceeds 5% |
| Prescribed-temperature flame, two transport stages | 1.728 | 1.699 | Inconclusive: timing spread exceeds 5% |

Each implementation performs 96 warm complete calculations, grouped into three
consecutive batches of 32. The speed ratio divides the median Cantera batch mean
by the median Julia batch mean. The conservative ratio divides the smallest
Cantera batch mean by the largest Julia batch mean. Both ratios must reach 0.95,
and each implementation's batch-mean spread must remain within 5%.

Each repetition reuses its loaded mechanism and starts from fresh flame state.
Timers include construction, initialization, solver setup, continuation, adaptive
solves and required numerical profiles. Imports, compilation, mechanism/sidecar
loading, validation and file I/O are excluded. One complete warmup and one explicit
collection follow loading; both are excluded from the warm samples. Automatic
garbage collection remains enabled and included in measured calculations. All
samples are retained, with one numerical thread per implementation. Linux elapsed
timers use `CLOCK_MONOTONIC_RAW`; both runtime clocks were checked against an
independent high-resolution host clock before measurement.

Every transport stage passes the full-profile and elemental-conservation checks
against independently refined Cantera references. Acceptance limits are 1% for
temperature profiles and free-flame speed, 5% for species profiles with a 1e-7
reference-peak floor, and 1e-6 for elemental conservation. Numerical replay checks
every first and warm result at relative tolerance 1e-12 and absolute tolerance
1e-14. Local tests pass 2,901 assertions at the measured native revision.

The [shared Julia calculation](../example/flames/source_flame_sequence.jl) is used
by the examples and the [paired benchmark](conservative_source_bench.py). Run the
benchmark with `--help` for its provenance and reference inputs. These results
cover the listed source conditions and transport sequences.
[Earlier measurements](results/cantera4_flames_5ca6641.json) retain their original
source revision and timing scope.