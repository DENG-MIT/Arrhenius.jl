# Scalar kinetics performance: fork vs. upstream, across mechanism size

Prepared while comparing `TJP-Karpowski/Arrhenius.jl` (`LLMImprovments` branch) against
`DENG-MIT/Arrhenius.jl` main, ahead of a possible `Kinetics.jl` performance PR. Not a merge
request by itself — this is the measurement and root-cause evidence for one.

## Question

Does the fork's independently-rewritten kinetics core (`StandardKinetics.jl` /
`ReorderedKinetics.jl` / `PackedKinetics.jl`) evaluate chemistry faster than upstream's
`Kinetics.jl`, and if so, does the gap depend on mechanism size?

## Method

- Machine: Apple M4 Pro, macOS, Julia 1.12.1, single-threaded.
- Mechanisms are byte-identical on both sides (sha256-verified): `h2o2.yaml`, `gri30.yaml`,
  `nc7_ver3.1_mech.yaml` (all predate the 2021-12-20 fork point and were never touched by
  either side since).
- Same tolerances both sides: `reltol=1e-6`, `abstol=1e-9`. No custom/analytical Jacobian on
  either side — both let `CVODE_BDF` (Sundials) use its own internal finite-difference
  Jacobian, so the comparison isolates the RHS/rate-evaluation cost, not either side's
  Jacobian infrastructure.
- RHS-only timing: 200 timed calls after 1 untimed warmup call, reusable in-place buffers
  (zero allocation on both sides), reporting the median.
- Full ignition timing: constant-pressure `CVODE_BDF` integration to the stated `t_end`, with
  `alg_hints=[:stiff]`, 1 untimed warmup solve, then 15 timed solves (`h2o2`, `gri30`) or 5
  timed solves (`nc7_ver3.1_mech`, due to per-solve cost at this size), reporting the median.

## Results

| Mechanism | species | reactions | upstream RHS (ns) | fork RHS (ns) | **RHS speedup** | upstream ignition | fork ignition | **ignition speedup** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `h2o2.yaml` | 9 | 28 | 667 | 708 | 0.94x | 1.121 ms / 173.3 KB | 1.190 ms / 187.7 KB | 0.94x |
| `gri30.yaml` | 53 | 325 | 6,708 | 5,792 | 1.16x | 16.99 ms / 327.1 KB | 14.82 ms / 334.2 KB | 1.15x |
| `nc7_ver3.1_mech.yaml` | 631 | 4,846 | 87,854 | 56,125 | **1.57x** | 3,745.9 ms / 3,663.4 KB | 2,936.2 ms / 3,682.1 KB | **1.28x** |

Conditions: `h2o2` — 1000 K, 1 atm, H2:O2 = 2:1. `gri30` — 1500 K, 10 atm, CH4:O2:N2 =
1:2:7.52. `nc7_ver3.1_mech` — 1200 K, 10 atm, NC7H16:O2:N2 = 1:11:41.36.

**Finding:** the two are at parity on a small mechanism (fork is actually ~6% slower on
`h2o2`), but the fork pulls increasingly ahead as `n_reactions` grows — reaching a **57%
faster scalar RHS** and a **28% faster full stiff-solve ignition** at 631-species/4,846-
reaction scale. The trend across all three points (0.94x → 1.16x → 1.57x) is monotonic in
mechanism size, not noise.

## Root cause (read-only code comparison, not yet acted on)

- Fork: `src/ReorderedKinetics.jl:261` (`wdot_func!`), called from `isobaricReactor`'s RHS
  closure in `src/ODE_functions.jl:302`. Third-body path at `ReorderedKinetics.jl:~225`.
- Upstream: `src/Kinetics.jl:463` (`_rate_factors!`), called from `wdot!` at
  `src/Kinetics.jl:682`. Third-body/falloff path at `Kinetics.jl:574-575,596`
  (`for i in reaction.index_three_body`, `dot(@view(reaction.efficiencies_coeffs[:,i]), C)`).

The likely driver is **reaction layout, not the sparse dot product itself** — a microbenchmark
isolating `dot(@view(sparse_col), C)` against a manual `nzrange`/`rowval` loop showed no
scaling difference with `n_species` (both flat, 83-125 ns), ruling that specific call out.

The fork's `ReorderedKinetics` buckets reactions by type **at parse time** — plain Arrhenius,
then third-body, then falloff, then Troe, each a contiguous index range
(e.g. `tB_range = reactions.index_three_body_start:reactions.index_fallOff_start-1`). The hot
loop walks each bucket as a plain contiguous range: branch-free, sequential memory access,
trivially vectorizable.

Upstream keeps reactions in original mechanism order and instead carries explicit **index
lists** (`reaction.index_three_body`, `reaction.index_falloff`) iterated with
`for i in index_three_body`. That's a gather/scatter access pattern into `kf[i]` and
`efficiencies_coeffs[:, i]` rather than a contiguous sweep. As a mechanism grows, those index
lists get sparser relative to the full reaction range, which plausibly explains a cost that
grows faster than linearly in `n_reactions` for upstream specifically — matching the observed
trend.

Secondary, unconfirmed: upstream's `_mass_action_product` (`Kinetics.jl:312`) re-derives a
"standard C1/C2/C3 order" classification from `orders`/`stoich` on **every call**, inside the
hot loop; the fork's bucketing encodes the equivalent classification once at parse time. Minor
per-reaction branching overhead — not confirmed as a major contributor, worth a look if the
primary fix doesn't close the whole gap.

## Suggested next step (not done here)

Reorder/bucket `Kinetics.jl`'s reactions by type at parse/load time into contiguous ranges
(mirroring the fork's `ReorderedKinetics` approach), replacing the `index_three_body` /
`index_falloff` gather-loop indirection, while keeping upstream's existing reader/Cantera-
semantics validation (which is more complete than the fork's and should not be replaced).
This is a pure-refactor candidate under upstream's own `AGENTS.md` convention: it must stay
bit-identical in states and solver counters, and should land as its own commit separate from
any policy change, with a baseline-vs-candidate `n_reactions` sweep (at least `h2o2`, `gri30`,
and a 500+-species mechanism) as the PR's evidence, same as this report.

## Caveats

- Single-threaded CPU only; GPU dispatch (`KernelAbstractions`/Metal) and the batched/turbo
  path are untracked by this report.
- Five samples on `nc7_ver3.1_mech` (vs. 15 on the smaller cases) due to per-solve cost
  (~3-4 s); min/median agreed closely in both cases (+/-1%), so this is unlikely to be noise.
- This has not been run on any mechanism between 325 and 4,846 reactions, so the exact shape
  of the scaling curve (linear-in-excess vs. some other growth) between `gri30` and
  `nc7_ver3.1_mech` is not characterized — only that it continues to grow.
