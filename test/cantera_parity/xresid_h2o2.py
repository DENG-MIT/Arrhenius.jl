#!/usr/bin/env python3
"""h2o2 regression test for xresid: round trip + converged-state residual ~0.

For four flame kinds (free, burner, fixed, stagnation -- both families plus
the imposed-temperature path) this builds an h2o2 mini, dumps the initial
guess (A) and the Cantera-converged solution (B) as canonical states, runs
xresid.py on them, and asserts:

  * the round trip passes (the tool itself refuses residuals otherwise), and
  * both codes' residuals agree on every shared row type at BOTH states
    (agreement_A/agreement_B <= 1e-9), and
  * the converged state's common relative residual is small in BOTH codes
    (rel_max_B_* <= 5e-2; measured ~5e-3 at Cantera's default tolerances).

Run:  python3 xresid_h2o2.py [--keep]      (julia via $XRESID_JULIA or PATH)
Exit 0 = all kinds pass; any failure prints the offending numbers and exits 1.
"""
import argparse
import json
import os
import shlex
import subprocess
import sys
import tempfile

for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np
import cantera as ct

HERE = os.path.dirname(os.path.abspath(__file__))
MECH = os.path.join(HERE, "..", "..", "mechanism", "h2o2")
XRESID = os.path.join(HERE, "xresid.py")
AGREE_MAX = 1e-9        # cross-code agreement on shared row types
RELB_MAX = 5e-2         # converged-state common relative residual
RTOL = 1e-11            # round-trip gate (xresid refuses above this)

BASE = {"mechanism": MECH, "transport": "mixture-averaged",
        "T": 300, "P_atm": 1, "X": "H2:2,O2:1,AR:3.76"}
SPECS = {
    "free": dict(BASE, id="h2o2_free_mini", kind="free", width=0.03),
    "burner": dict(BASE, id="h2o2_burner_mini", kind="burner", width=0.02,
                   mdot=0.1),
    "stagnation": dict(BASE, id="h2o2_stagnation_mini", kind="stagnation",
                       width=0.02, mdot=0.1, T_surface=300),
    # profile filled in at run time from the burner solution
    "fixed": dict(BASE, id="h2o2_fixed_mini", kind="fixed", width=0.02,
                  mdot=0.1, profile_z=None, profile_T=None),
}


def build(spec):
    gas = ct.Solution(spec["mechanism"] + ".yaml")
    kind = spec["kind"]
    if kind == "free":
        gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
        f = ct.FreeFlame(gas, width=spec["width"])
    elif kind in ("burner", "fixed"):
        gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
        f = ct.BurnerFlame(gas, width=spec["width"])
        f.burner.mdot = spec["mdot"]
    elif kind == "stagnation":
        gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
        f = ct.ImpingingJet(gas, width=spec["width"])
        f.inlet.mdot = spec["mdot"]
        f.surface.T = spec["T_surface"]
    else:
        raise ValueError(kind)
    f.transport_model = spec["transport"]
    # Guess with the energy equation ON (Sim1D.solve's auto path); the fixed
    # kind turns it off AFTER guessing (guessing energy-off leaves u = 0, a
    # degenerate state whose inlet row differs from Cantera's by the
    # state-vs-spec mdot anchor on Y_in).
    if kind != "stagnation":
        f.set_initial_guess()
    else:
        f.set_initial_guess("equil")
    f.energy_enabled = kind != "fixed"
    if kind == "fixed":
        z = np.asarray(spec["profile_z"], float)
        f.flame.set_fixed_temp_profile(z / z[-1],
                                       np.asarray(spec["profile_T"], float))
    return f


def dump(f, spec, path, imposed=False):
    d = f.flame
    T = np.asarray(d.values("T"), float)
    if imposed:
        g = np.asarray(d.grid, float)
        zz = (g - g[0]) / (g[-1] - g[0])
        pz = np.asarray(spec["profile_z"], float)
        T = np.interp(zz, pz / pz[-1], np.asarray(spec["profile_T"], float))
    kw = dict(grid=np.asarray(d.grid, float), T=T,
              velocity=np.asarray(d.values("velocity"), float),
              Y=np.column_stack([np.asarray(d.values(s), float)
                                 for s in f.gas.species_names]),
              species_names_utf8="\n".join(f.gas.species_names).encode())
    if spec["kind"] in ("stagnation",):
        kw["spreadRate"] = np.asarray(d.values("spreadRate"), float)
        kw["Lambda"] = np.asarray(d.values("Lambda"), float)
    np.savez(path, **kw)


def run_case(kind, spec, tmpdir):
    f = build(spec)
    dump(f, spec, os.path.join(tmpdir, f"{kind}_A.npz"),
         imposed=spec["kind"] == "fixed")
    f.solve(loglevel=0, auto=spec["kind"] != "fixed")
    dump(f, spec, os.path.join(tmpdir, f"{kind}_B.npz"))
    out = os.path.join(tmpdir, f"out_{kind}")
    res = os.path.join(tmpdir, f"{kind}_res.json")
    specfile = os.path.join(tmpdir, f"spec_{kind}.json")
    with open(specfile, "w") as fh:
        json.dump(spec, fh)
    julia = shlex.split(os.environ.get("XRESID_JULIA", "julia"))
    proc = subprocess.run(
        [sys.executable, XRESID, specfile,
         os.path.join(tmpdir, f"{kind}_A.npz"), os.path.join(tmpdir, f"{kind}_B.npz"),
         "--out", out, "--json", res, "--rtol", str(RTOL), "--julia", *julia],
        capture_output=True, text=True)
    if proc.returncode:
        print(proc.stdout[-2000:])
        print(proc.stderr[-2000:])
        return None, f"xresid exited {proc.returncode}"
    return json.load(open(res)), None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep", action="store_true", help="keep generated states and reports")
    parser.add_argument("--kinds", nargs="+", choices=tuple(SPECS), default=list(SPECS),
                        help="select cases; defaults to all four flame kinds")
    args = parser.parse_args()
    tmpdir = tempfile.mkdtemp(prefix="xresid_h2o2_")
    failures = []
    # The fixed kind's imposed profile comes from the burner solution.
    if "fixed" in args.kinds:
        fb = build(SPECS["burner"])
        fb.solve(loglevel=0, auto=True)
        gz = np.asarray(fb.flame.grid, float)
        SPECS["fixed"]["profile_z"] = ((gz - gz[0]) / (gz[-1] - gz[0])).tolist()
        SPECS["fixed"]["profile_T"] = np.asarray(
            fb.flame.values("T"), float).tolist()
        del fb
    for kind in args.kinds:
        spec = SPECS[kind]
        data, err = run_case(kind, spec, tmpdir)
        if err:
            print(f"FAIL {kind}: {err}")
            failures.append(kind)
            continue
        agg = data["aggregates"]
        rt = max(v for chain in data["roundtrip"].values()
                 for v in chain.values())
        checks = {
            "roundtrip": (rt, RTOL),
            "agreement_A": (agg["agreement_A"], AGREE_MAX),
            "agreement_B": (agg["agreement_B"], AGREE_MAX),
            "rel_max_B_native": (agg["rel_max_B_native"], RELB_MAX),
            "rel_max_B_cantera": (agg["rel_max_B_cantera"], RELB_MAX),
        }
        bad = {k: v for k, (v, lim) in checks.items() if v > lim}
        line = " ".join(f"{k}={v:.2e}" for k, (v, _) in checks.items())
        if bad:
            print(f"FAIL {kind}: {line} (limits exceeded: {sorted(bad)})")
            failures.append(kind)
        else:
            print(f"pass {kind}: {line}")
    if args.keep:
        print(f"kept {tmpdir}")
    if failures:
        print(f"xresid_h2o2: FAILED kinds: {failures}")
        sys.exit(1)
    print(f"xresid_h2o2: all selected kinds pass: {args.kinds}")


if __name__ == "__main__":
    main()
