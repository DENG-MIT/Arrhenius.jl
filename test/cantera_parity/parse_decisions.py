#!/usr/bin/env python
"""Parse and align Cantera vs native solver-decision logs (test/cantera_parity).

Usage:
    parse_decisions.py ct <ct_solve.log>       Cantera log (loglevel >= 3)
    parse_decisions.py native <native_solve.log>
    parse_decisions.py both <ct_solve.log> <native_solve.log>

Prints one line per solver decision: grid phase, steady-attempt outcome,
Newton iteration (damping alpha, bound fraction, s0/s1 or native s0), Jacobian
evaluations (cumulative), time-step accept/fail with dt, refine events.  The
`both` mode prints the aligned opening of the two sequences.
"""
import re
import sys

CT_ROW = re.compile(r"^\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/(\d+)\s*$")


def parse_ct(path):
    events = []
    njac = 0
    for line in open(path):
        s = line.strip()
        if s.startswith("Solving on"):
            events.append(("phase", s))
        elif s.startswith("Attempt Newton solution"):
            events.append(("steady_start", None))
        elif s.startswith("Attempt") and "timesteps" in s:
            events.append(("ts_batch", s))
        elif m := CT_ROW.match(line):
            it, alpha, fbound, lss, ls0, ls1, nj, age, agemax = m.groups()
            njac = int(nj)
            events.append(("iter", dict(m=int(it), alpha=float(alpha),
                fbound=float(fbound), logss=float(lss), logs0=float(ls0),
                logs1=float(ls1), njac=njac, age=f"{age}/{agemax}")))
        elif s.startswith("Timestep (") and "dt=" in s:
            events.append(("ts_step", s))
        elif s.startswith("Timestep failed") or s.startswith("Timestep ("):
            events.append(("ts_event", s))
        elif "No damped step" in s or "Re-evaluating Jacobian" in s or \
                "Maximum Jacobian age" in s or "Damping coefficient found" in s or \
                "Resetting negative" in s:
            events.append(("note", s))
        elif s.startswith("Newton steady-state solve"):
            events.append(("steady_end", "succeeded" in s))
        elif s.startswith("Refining grid") or s.startswith("no new points"):
            events.append(("refine", s))
    return events


def parse_native(path):
    events = []
    njac = 0
    for line in open(path):
        s = line.strip()
        if s.startswith("Solving on"):
            events.append(("phase", s))
        elif s.startswith("FLAMEDBG jacobian"):
            njac += 1
            events.append(("jac", s))
        elif s.startswith("FLAMEDBG steady_attempt"):
            events.append(("steady_end", "success=true" in s))
        elif s.startswith("FLAMEDBG steady_trial"):
            events.append(("steady_trial", s))
        elif s.startswith("FLAMEDBG timestep"):
            events.append(("ts_step", s))
        elif s.startswith("FLAMEDBG iter="):
            m = re.match(r"FLAMEDBG iter=(\d+) fbound=(\S+) resnorm=(\S+)", s)
            if m:
                events.append(("iter", dict(it=int(m.group(1)),
                    fbound=float(m.group(2)), resnorm=float(m.group(3)), njac=njac)))
            m2 = re.match(r"FLAMEDBG iter=(\d+) trials alpha_end=(\S+) accepted=(\S+) s0=(\S+)", s)
            if m2:
                events.append(("trial", dict(it=int(m2.group(1)),
                    alpha=float(m2.group(2)), accepted=m2.group(3) == "true",
                    s0=float(m2.group(4)))))
        elif s.startswith("Retrying") or s.startswith("Time stepping failed") or \
                s.startswith("Expanding domain"):
            events.append(("note", s))
    return events


def show(events, limit=None):
    n = 0
    for kind, data in events:
        if kind == "phase":
            print(f"\n== {data} ==")
        elif kind == "steady_start":
            print("  steady attempt:")
        elif kind == "iter":
            print(f"    iter {data.get('it', data.get('m'))}: {data}")
            n += 1
        elif kind == "trial":
            print(f"      -> {data}")
        elif kind in ("ts_step", "ts_batch", "ts_event"):
            print(f"  {data}")
        elif kind == "steady_end":
            print(f"  steady {'SUCCESS' if data else 'FAILED'}")
        elif kind == "steady_trial":
            print(f"  {data}")
        elif kind == "refine":
            print(f"== {data} ==")
        elif kind == "jac":
            pass
        else:
            print(f"  [{kind}] {data}")
        if limit and n > limit:
            break


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "ct":
        show(parse_ct(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else None)
    elif mode == "native":
        show(parse_native(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else None)
    else:
        print("use ct|native")
