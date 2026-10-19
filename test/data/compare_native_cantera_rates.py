#!/usr/bin/env python3
"""Compare Arrhenius and Cantera point rate constants on selected gas mechanisms."""
import os
import subprocess
import sys
from pathlib import Path
import cantera as ct
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
STATES = ((500.0, 0.1 * ct.one_atm),
          (1200.0, ct.one_atm),
          (2500.0, 10.0 * ct.one_atm))
paths = [Path(p).resolve() for p in sys.argv[1:]]
if not paths:
    raise SystemExit("usage: compare_native_cantera_rates.py mechanism.yaml ...")
jl = os.environ.get("JL", "julia")
project = os.environ.get("ARRHENIUS_PROJECT", os.environ.get("AJ", str(ROOT)))
result = subprocess.run(
    [jl, f"--project={project}", str(ROOT / "test/data/native_rate_dump.jl"),
     *map(str, paths)],
    check=False, text=True, capture_output=True,
)
if result.returncode:
    sys.stderr.write(result.stderr)
    raise SystemExit(result.returncode)
dump = result.stdout
native = {}
for line in dump.splitlines():
    if not line or line.startswith("#"):
        continue
    mech, T, P, index, forward, reverse = line.split("|")
    native[(mech, float(T), float(P), int(index))] = (float(forward), float(reverse))

print("mechanism,reactions,states,compared,max_rel_forward,max_rel_reverse,diff_gt_1e-12,worst_reaction")
for path in paths:
    gas = ct.Solution(str(path))
    X = np.full(gas.n_species, 1.0 / gas.n_species)
    count = 0
    max_forward = max_reverse = 0.0
    worst = (0.0, 0)
    over = 0
    for T, P in STATES:
        gas.TPX = T, P, X
        can_f = gas.forward_rate_constants
        can_r = gas.reverse_rate_constants
        tb = gas.third_body_concentrations
        for i, rxn in enumerate(gas.reactions(), start=1):
            factor = 1.0
            third_body = getattr(rxn, "third_body", None)
            if third_body is not None and third_body.mass_action:
                factor = tb[i - 1]
            expected = (float(can_f[i - 1] * factor), float(can_r[i - 1] * factor))
            actual = native[(path.name, T, P, i)]
            rels = [abs(a - c) / max(abs(a), abs(c), 1.0e-300)
                    for a, c in zip(actual, expected)]
            max_forward = max(max_forward, rels[0])
            max_reverse = max(max_reverse, rels[1])
            this = max(rels)
            if this > worst[0]:
                worst = (this, i)
            over += this > 1.0e-12
            count += 1
    print(f"{path.name},{gas.n_reactions},{len(STATES)},{count},"
          f"{max_forward:.6g},{max_reverse:.6g},{over},{worst[1]}:{worst[0]:.6g}")
