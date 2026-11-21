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

## Follow-up: isolating reordering from the equilibrium-skip fix (2026-11-21)

`eab69a7` (wiring `_ReversibleRatePlan` into `wdot!` by default) landed after this report and,
per its own commit message, already closed most of the fork-vs-upstream gap measured above with
a 3-line, non-reordering change. This follow-up isolates what reordering itself still has left
to contribute on top of that fix, on branch `reactions-reorder` (`upstream-port/kinetics-performance`
+ `ReorderedReactions`, a bucketed-by-kind reaction layout carrying an original-order permutation).

**Reordering alone (bucketing `index_three_body`/`index_falloff` into contiguous ranges,
generic evaluator unchanged): noise-level**, ±0-3% on `h2o2`, `gri30`, `stagni_2023`, and
`nc7_ver3.1_mech` — consistent with `eab69a7`'s own dropped `BucketedEfficiencyCache` finding:
these mechanisms' third-body/falloff fraction is too small for removing that gather/scatter to
matter once the equilibrium work is already skipped for irreversible reactions. The
`ReorderedKinetics`-vs-upstream delta measured above was overwhelmingly the equilibrium-skip
effect, not reordering.

**The fork's kernel (`wdot_func!`) is still consistently 7-12% faster than this repo's generic
evaluator (`_rate_factors!`/`_mass_action!`/`wdot!`) at the raw kernel level**, flat across
mechanism size (28 to 4,846 reactions) — not a locality/scaling effect, a roughly constant
per-call difference:

| Mechanism | reactions | fork `wdot_func!` (ns) | upstream kernel (ns) | fork/upstream |
|---|---:|---:|---:|---:|
| h2o2 | 28 | 541 | 584 | 0.93x |
| h2_li_2004 | 21 | 458 | 500 | 0.92x |
| gri30 | 325 | 5959 | 6750 | 0.88x |
| nc7_ver3.1_mech | 4846 | 58771 | 63416 | 0.93x |

Ruled out: generic-path overhead from supporting `log_rate_data`/`temperature_cache` when
unused. `@code_typed` on `_rate_factors!` called the way this benchmark calls it (those
keywords omitted) shows zero residual `isnothing`/`Union` branches — Julia already specializes
them away at compile time, so there is no runtime-branch cost to eliminate there.

One real, explained piece of the gap was identified and ported: the fork's `ReorderedKinetics`
carries a separate `vk_rev` matrix sized to just the reversible reactions, letting
`wdot_func!` skip the delta_s/delta_h/equilibrium-constant computation for irreversible
reactions via plain contiguous-range access rather than an `is_reversible`/`_ReversibleRatePlan`
lookup per reaction. Porting this (`ReorderedReactionsKernels.jl`: `ReorderedReactions` now
sub-buckets its plain-reaction group by reversibility too, and dispatches a specialized
`_rate_factors!` that walks `index_rev_start:n_reactions` directly) reproduces this specific
piece of the gap, and *only* this piece:

| Mechanism | reversible fraction | kernel speedup from the specialized path |
|---|---:|---:|
| h2o2 | 100% | ~1.0x (nothing to skip) |
| gri30 | ~95% | ~1.0x (nothing to skip) |
| h2_li_2004 | lower | ~1.08-1.09x |
| nc7_ver3.1_mech | ~4.5% | ~1.04-1.10x (noisier at this mechanism's per-call cost) |

This is structurally expected, not a bug: on a mechanism that's almost entirely reversible
(`h2o2`, `gri30`), there is no irreversible prefix to skip, matching `eab69a7`'s own observation
that its equilibrium-skip fix likewise made no difference on those two. It also means `gri30`'s
share of the fork's measured 0.88x is **not** explained by the reversible-prefix-skip — some
other, smaller difference (candidates: the fork's `k_Pi_C_rectantOrder!`-style concentration
product vs. this repo's C1/C2/C3-classified `_mass_action_product`, or the single fused
`wdot_func!` vs. the generic path's three-function split) accounts for the residual there, and
was not further isolated.

**Net assessment**: reordering by itself is not worth shipping for raw speed on top of
`eab69a7`. The reversible-prefix-skip kernel is a real, explained, correctly-scoped ~4-10% win
specifically on mechanisms with a meaningful irreversible fraction (detailed/large mechanisms
are more likely to have one), shipped alongside — not in place of — the generic evaluator via
ordinary multiple dispatch, with every other kinetics feature (PLOG, linear-Burke,
Blowers-Masel, `log_rate_data`, `temperature_cache`, `rate_multipliers`) unchanged. It is not a
uniform ~10% win across all mechanisms, and does not fully explain the fork's own measured
advantage on reversible-heavy mechanisms like `gri30`.

## Follow-up: profiling the residual gap, and a fused kernel (2026-11-21, continued)

Profiled both kernels directly (`Profile.jl`, 200,000 calls each, equal ~1.1M-sample budgets,
`gri30`) to look for the cause of the still-open `gri30` gap. No single dominant cause: the
generic path shows real `kwcall`/`NamedTuple` construction self-time from its keyword-argument
interface (absent from the fork's all-positional `wdot_func!`), but it's small in absolute
terms (~6,948 self samples / 200,000 calls, a few ns/call) — not enough alone to explain a
~700ns/call gap. The rest is spread thin across many small implementation-style differences
(`getindex`/`_setindex!`/`*`/`reinterpret` self-time differ in both directions; no standout
line). Conclusion: no isolated fix would close the rest of the gap on its own.

Re-reading the fork's `k_Pi_C_rectantOrder!` line-by-line against this repo's
`_mass_action_product` surfaced something a pure performance read would have missed: they are
not equivalent algorithms. For a repeated-species order-2 term with a negative concentration,
Cantera's rule (and `_mass_action_product`'s, exercised by `test/negative_mass_action.jl`'s 61
tests) gives `0.0`; the fork's `order==2 ? Ci*Ci : ...` gives `Ci^2` (positive) — the opposite
of the tested, physically-motivated behavior for solver overshoot. **Decision (confirmed with
the task owner): keep `_mass_action_product` unchanged, port the fork's structure (not its
math) around it.**

Added a second specialization, `wdot!(reaction::ReorderedReactions, ...)`
(`ReorderedReactionsKernels.jl`), fusing the generic path's three separate full-mechanism
passes (`_rate_factors!` → `_mass_action!` writing a full `rates_of_progress` array →
`wdot!`'s own scatter reading it back) into one pass — mirroring the fork's single-function
`wdot_func!` shape, with `_mass_action_product` reused exactly as-is. `kf[i]` is overwritten in
place to hold the net rate of progress for the rest of the function, the same proxy-reuse
trick the fork uses for `kf .-= kr; vk_kf!(wdot,vk,kf)` (documented in the method's docstring).

Result (median of 200 samples, 30 for `nc7_ver3.1_mech`, full `IdealGasReactor` RHS, three
repeated runs):

| Mechanism | before fusion | after fusion |
|---|---:|---:|
| h2o2 | ~1.0x | ~1.05-1.06x |
| gri30 | ~1.0x | **~1.02-1.03x** |
| h2_li_2004 | ~1.08-1.09x | ~1.06x |
| nc7_ver3.1_mech | ~1.04-1.10x | ~1.09-1.11x |

Real and repeatable (all four mechanisms improved, consistently, across three separate runs —
not noise), but `gri30` is still far from the fork's full 0.88x (would need ~1.14x to match).
This is the expected cost of the correctness decision above: the remaining gap most plausibly
lives in `_mass_action_product`'s extra negative-concentration branching relative to the fork's
simpler (and, in that one case, wrong) algorithm, which this port deliberately did not adopt.

**Updated net assessment**: the fused kernel is a modest, broad, genuinely explained
improvement (not a partial fix limited to low-reversibility mechanisms, like the first
specialization). It still does not close `gri30`'s full gap, and closing the rest would require
either accepting the fork's less-robust mass-action algorithm (rejected) or a different,
not-yet-identified safe optimization of `_mass_action_product` itself.
