#!/usr/bin/env python3
"""xresid — cross-residual check: both codes' residuals at both states.

Usage (from a directory where the mechanism YAML resolves):
    xresid.py <spec.json> <stateA.npz> <stateB.npz> [--scales B|A]
              [--out DIR] [--json FILE] [--rtol 1e-11] [--top 30]
              [--discretization finite_difference|conservative]

States are canonical npz (physical units): grid (N,), T (N,) [K], Y (N,n)
points-major, velocity (N,) [m/s], spreadRate/Lambda (N,) for counterflow
kinds; ct-style dumps (values/<name>, (n,N) Y) also load. stateB is normally
the converged solution and provides the scales.

The tool verifies the state transfer round-trips through both codes to --rtol
per component: native canonical->install->extract, canonical->Cantera
set_values->read back, native-extract->Cantera->read back, and
Cantera-readback->native->extract. If any chain fails the tool REFUSES to
report residuals (exit 2).

Report per component (T, each species, V, Lambda, mass flux): max |r| in
physical units (K/s, 1/s, kg/m2/s, ...) and where (point index j, z, species),
plus the common relative norm

    rel(row, eqtype) = max_j |r_row(j)| / S(row, eqtype),
    S = max over points AND additive terms of |term| (physical units)
        evaluated on the scales state (default stateB),

the same S for both codes, so the norms are comparable across codes.
Cross-code agreement per row is max_j |r_native - r_cantera| / S over
row/point pairs where both codes solve the same equation type; structural
differences (native boundary excess-species closure row, the vacuous Cantera
free-flame inlet mass-flux row, anchor index disagreement, twin-wall species
form) are listed and excluded from the agreement, never silently dropped.
"""
import argparse
import json
import os
import shlex
import subprocess
import sys

for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np
import cantera as ct

HERE = os.path.dirname(os.path.abspath(__file__))          # .../test/cantera_parity
PROJECT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
CF_KINDS = ("diffusion", "premixed", "twin", "stagnation")
# native row order: T, species..., [V, Lambda,] mdot ("mdot" = CT "velocity",
# native V = CT "spreadRate" component in ct 4.0)
CT_COMPONENT = {"mdot": "velocity", "V": "spreadRate", "Lambda": "Lambda", "T": "T"}


def mechanism_path(spec):
    path = os.fspath(spec["mechanism"])
    if not os.path.isabs(path):
        path = os.path.join(PROJECT, path)
    return path + ".yaml"


def load_spec(path):
    spec = json.load(open(path))
    if "kind" not in spec and "spec" in spec:
        spec = spec["spec"]          # reproducer jsons nest under "spec"
    return spec


def canonical(path, species):
    state = np.load(path)
    z = np.asarray(state["grid"], float)
    T = (np.asarray(state["T"], float) if "T" in state
         else np.asarray(state["values/T"], float))
    u = (np.asarray(state["velocity"], float) if "velocity" in state
         else np.asarray(state["values/velocity"], float))
    n, N = len(species), len(z)
    if "Y" in state:
        Y = np.asarray(state["Y"], float)
        if Y.shape == (n, N) and Y.shape != (N, n):
            Y = Y.T
    else:
        Y = np.column_stack([np.asarray(state["values/" + s], float)
                             for s in species])
    if Y.shape != (N, n):
        raise SystemExit(f"Y shape {Y.shape} != ({N}, {n}) for {path}")
    out = {"grid": z, "T": T, "velocity": u, "Y": Y}
    for key in ("spreadRate", "Lambda"):
        if key in state:
            out[key] = np.asarray(state[key], float)
        elif "values/" + key in state:
            out[key] = np.asarray(state["values/" + key], float)
    return out


def save_canonical(path, can, species):
    np.savez(path, grid=can["grid"], T=can["T"], velocity=can["velocity"],
             Y=can["Y"],
             species_names_utf8=np.frombuffer("\n".join(species).encode(), np.uint8),
             **{k: can[k] for k in ("spreadRate", "Lambda") if k in can})


def rtol(a, b):
    """Componentwise round-trip error: max|a-b| / max(max|a|, max|b|)."""
    denom = max(np.max(np.abs(a)), np.max(np.abs(b)), 1e-300)
    return float(np.max(np.abs(a - b)) / denom)


def roundtrip_all(can, back):
    rt = {}
    for key in ("T", "velocity", "spreadRate", "Lambda"):
        if key in back:
            rt[key] = rtol(can[key], back[key])
    rt["Y"] = rtol(can["Y"], back["Y"])
    return rt


# ---------------------------------------------------------------- Cantera side

def build(spec, grid):
    """Cantera flame on the exact grid (constructor kwarg; resizing corrupts)."""
    gas = ct.Solution(mechanism_path(spec))
    kind = spec["kind"]
    if kind == "diffusion":
        gas.TP = spec["T"], spec["P_atm"] * ct.one_atm
        f = ct.CounterflowDiffusionFlame(gas, grid=grid)
        f.fuel_inlet.mdot, f.fuel_inlet.X, f.fuel_inlet.T = (
            spec["mdot_fuel"], spec["fuel"], spec["T"])
        f.oxidizer_inlet.mdot, f.oxidizer_inlet.X, f.oxidizer_inlet.T = (
            spec["mdot_oxidizer"], spec["oxidizer"],
            spec.get("T_oxidizer", spec["T"]))
    elif kind == "premixed":
        gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
        f = ct.CounterflowPremixedFlame(gas, grid=grid)
        f.reactants.mdot = spec["mdot"]
        gas.equilibrate("HP")
        f.products.X, f.products.T = gas.X, gas.T
        f.products.mdot = spec["mdot_products"]
    elif kind == "twin":
        gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
        f = ct.CounterflowTwinPremixedFlame(gas, grid=grid)
        f.reactants.mdot = spec["mdot"]
    elif kind == "stagnation":
        gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
        f = ct.ImpingingJet(gas, grid=grid)
        f.inlet.mdot = spec["mdot"]
        f.surface.T = spec["T_surface"]
        f.set_initial_guess("equil")
    elif kind == "free":
        gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
        f = ct.FreeFlame(gas, grid=grid)
        f.set_initial_guess()
    elif kind in ("burner", "fixed"):
        gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
        f = ct.BurnerFlame(gas, grid=grid)
        f.burner.mdot = spec["mdot"]
        f.set_initial_guess()
        if kind == "fixed":
            z, T = np.asarray(spec["profile_z"]), np.asarray(spec["profile_T"])
            # normalize ABSOLUTE z (val specs) or already-normalized z alike
            f.flame.set_fixed_temp_profile((z - z[0]) / (z[-1] - z[0]), T)
    else:
        raise SystemExit(f"unsupported kind {kind}")
    f.transport_model = spec["transport"]
    # Flow1D::setupGrid->resize leaves energy disabled past the first point;
    # the fixed kind wants it off everywhere (imposed profile).
    f.energy_enabled = kind != "fixed"
    if kind == "fixed":
        # Flow1D::_finalize fills m_fixedtemp (linearInterp of the profile,
        # Flow1D.cpp:287-297) but no guess/eval path calls it for a non-free
        # flame -- only Sim1D::solve does. Sim1D::setFixedTemperature ends with
        # resize()+finalize() (Sim1D.cpp:615-617) and leaves a non-free
        # domain's grid and state untouched: the python-reachable trigger.
        f.fixed_temperature = float(np.asarray(spec["profile_T"])[-1])
    return f


def install(f, can):
    d = f.flame
    profiles = {"T": can["T"], "velocity": can["velocity"]}
    profiles.update({s: can["Y"][:, k] for k, s in enumerate(f.gas.species_names)})
    for key in ("spreadRate", "Lambda"):
        if key in can:
            profiles[key] = can[key]
    for name, values in profiles.items():
        d.set_values(name, np.asarray(values, float))


def readback(f, spread):
    """spread=True for counterflow kinds: read V/Lambda back too. (The flow
    domain reports the full component layout even for components it does not
    use, so the kind decides, not the component names.)"""
    d = f.flame
    out = {"grid": np.asarray(d.grid, float),
           "T": np.asarray(d.values("T"), float),
           "velocity": np.asarray(d.values("velocity"), float)}
    out["Y"] = np.column_stack([np.asarray(d.values(s), float)
                                for s in f.gas.species_names])
    if spread:
        out["spreadRate"] = np.asarray(d.values("spreadRate"), float)
        out["Lambda"] = np.asarray(d.values("Lambda"), float)
    return out


def ct_residuals(f, rows):
    """Steady residual of the flow domain in native row order, (B, N)."""
    f.eval(rdt=0.0)
    return np.column_stack([np.asarray(f.flame.residuals(CT_COMPONENT.get(r, r)), float)
                            for r in rows]).T


def ct_types(f, spec, R, can):
    """Mirror of the native row_eqtypes with Cantera's structural exceptions."""
    kind = spec["kind"]
    n, N = len(f.gas.species_names), len(can["grid"])
    types = np.empty((R.shape[0], N), dtype=object)
    cf = kind in CF_KINDS
    # Cantera's excess species = argmax of the boundary mass fractions,
    # re-selected on every property update (Flow1D::updateProperties); numpy
    # argmax returns the first maximum, same as ranges::max_element.
    excess_left = int(np.argmax(can["Y"][0]))
    excess_right = int(np.argmax(can["Y"][-1]))
    # Cantera's free-flame anchor: T - T_fixed sits in the mass-flux row
    # (Flow1D::evalContinuity, m_isFree branch).
    anchor = None
    if kind == "free" and np.isfinite(f.fixed_temperature):
        d = R[-1] - (can["T"] - f.fixed_temperature)
        hit = np.nonzero(d == 0.0)[0]
        anchor = int(hit[0]) if len(hit) else None
    for j in range(N):
        boundary = j in (0, N - 1)
        types[0, j] = "T_bc" if boundary or kind == "fixed" else "energy"
        for k in range(n):
            if boundary and k == (excess_left if j == 0 else excess_right):
                types[1 + k, j] = "closure_ct"
            else:
                types[1 + k, j] = ("flux_bc" if boundary else "species")
        if cf:
            types[n + 1, j] = "zg_V" if boundary else "momentum"
            types[n + 2, j] = "lambda_left" if j == 0 else "lambda"
            types[-1, j] = "mdot_bc" if j == N - 1 else "continuity"
        elif kind in ("burner", "fixed"):
            types[-1, j] = "mdot_bc" if j == 0 else "continuity"
        elif kind == "free":
            if j == anchor:
                types[-1, j] = "T_anchor"
            elif j == N - 1:
                types[-1, j] = "zero_grad"
            else:
                types[-1, j] = "continuity"
    if kind == "twin":
        # Symmetry plane: V/T zero-gradient; species stay flux-form in Cantera.
        types[0, N - 1] = "zg_T"
        types[n + 1, N - 1] = "zg_V"
    return types, anchor, excess_left, excess_right


# ------------------------------------------------------------------ reporting

def group_rows(R, types, rows):
    """(row, eqtype) -> (max |r|, j (1-based))."""
    groups = {}
    for i, row in enumerate(rows):
        for t in np.unique(types[i]):
            js = np.nonzero(types[i] == t)[0]
            vals = np.abs(R[i, js])
            m = int(np.argmax(vals))
            groups[(row, t)] = (float(vals[m]), int(js[m]) + 1)
    return groups


def fmt(v):
    return f"{v:.2e}" if v is not None and np.isfinite(v) else "-"


def col(v, width, spec="e"):
    return "-" if v is None else f"{v:.2{spec}}" if spec == "e" else f"{v:.{spec}}".rjust(width)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("stateA")
    ap.add_argument("stateB")
    ap.add_argument("--scales", choices=("B", "A"), default="B")
    ap.add_argument("--out", default="xresid_out")
    ap.add_argument("--json")
    ap.add_argument("--rtol", type=float, default=1e-11)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--discretization", default="finite_difference")
    ap.add_argument("--julia", default=os.environ.get("XRESID_JULIA", "julia"))
    args = ap.parse_args()

    spec = load_spec(args.spec)
    spread = spec["kind"] in CF_KINDS
    os.makedirs(args.out, exist_ok=True)
    julia = shlex.split(args.julia) + [f"--project={PROJECT}",
                                       os.path.join(HERE, "xresid_native.jl")]

    def run_native(jobs, tag):
        path = os.path.join(args.out, f"jobs_{tag}.json")
        with open(path, "w") as fh:
            json.dump(jobs, fh)
        proc = subprocess.run(julia + [path], capture_output=True, text=True)
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        if proc.returncode:
            raise SystemExit(f"xresid_native {tag} failed ({proc.returncode})")

    def jload(tag, job_id):
        suffix = {"residual": "res", "extract": "extract"}[tag]
        return json.load(open(os.path.join(args.out, f"{job_id}_{suffix}.json")))

    # -- 1. canonical states; native self round trips (julia run 1)
    species = list(ct.Solution(mechanism_path(spec)).species_names)
    A = canonical(args.stateA, species)
    B = canonical(args.stateB, species)
    if spec["kind"] in CF_KINDS and ("spreadRate" not in A or "spreadRate" not in B):
        raise SystemExit("counterflow states need spreadRate and Lambda")
    save_canonical(os.path.join(args.out, "A_can.npz"), A, species)
    save_canonical(os.path.join(args.out, "B_can.npz"), B, species)

    jobs = {"spec": args.spec, "discretization": args.discretization,
            "extract": [
                {"id": "A", "state": f"{args.out}/A_can.npz",
                 "out": f"{args.out}/A_nat.npz", "out_json": f"{args.out}/A_extract.json"},
                {"id": "B", "state": f"{args.out}/B_can.npz",
                 "out": f"{args.out}/B_nat.npz", "out_json": f"{args.out}/B_extract.json"}]}
    run_native(jobs, "extract")

    # -- 2. Cantera self round trips + native->CT chain (on A)
    ct_roundtrips = {}
    fA = build(spec, A["grid"])
    install(fA, A)
    ct_roundtrips["ct_self_A"] = roundtrip_all(A, readback(fA, spread))
    fB = build(spec, B["grid"])
    install(fB, B)
    backB = readback(fB, spread)
    ct_roundtrips["ct_self_B"] = roundtrip_all(B, backB)
    A_nat = canonical(os.path.join(args.out, "A_nat.npz"), species)
    f2 = build(spec, A_nat["grid"])
    install(f2, A_nat)
    ct_roundtrips["native_to_ct_A"] = roundtrip_all(A_nat, readback(f2, spread))

    # -- 3. CT->native chain (on B) and both residuals (julia run 2)
    save_canonical(os.path.join(args.out, "B_ct.npz"), backB, species)
    jobs = {"spec": args.spec, "discretization": args.discretization,
            "extract": [{"id": "B2", "state": f"{args.out}/B_ct.npz",
                         "out": f"{args.out}/B_nat2.npz",
                         "out_json": f"{args.out}/B2_extract.json"}],
            "residual": [
                {"id": "A", "state": f"{args.out}/A_can.npz",
                 "out_npz": f"{args.out}/A_res.npz", "out_json": f"{args.out}/A_res.json"},
                {"id": "B", "state": f"{args.out}/B_can.npz",
                 "out_npz": f"{args.out}/B_res.npz", "out_json": f"{args.out}/B_res.json"}]}
    run_native(jobs, "residual")
    resA, resB = jload("residual", "A"), jload("residual", "B")
    rt = {"native_self_A": resA["roundtrip"],
          "native_self_B": resB["roundtrip"],
          "ct_to_native_B": jload("extract", "B2")["roundtrip"],
          **ct_roundtrips}

    worst = max((v for chain in rt.values() for v in chain.values()), default=0.0)
    print(f"== xresid {spec.get('id', spec['kind'])} kind={spec['kind']} "
          f"N_A={len(A['grid'])} N_B={len(B['grid'])} scales={args.scales}")
    print("-- round trip (max componentwise relative)")
    for chain, vals in rt.items():
        print(f"   {chain:16s}: " + " ".join(f"{k}={v:.2e}" for k, v in sorted(vals.items())))
    if worst > args.rtol:
        print(f"REFUSING to report residuals: round trip {worst:.2e} > {args.rtol:.0e}")
        if args.json:
            json.dump({"refusal": "round trip", "roundtrip": rt, "rtol": args.rtol},
                      open(args.json, "w"), indent=1)
        sys.exit(2)

    # -- 4. residuals in both codes at both states, common relative norm
    rows = ["T"] + species + (["V", "Lambda", "mdot"] if spec["kind"] in CF_KINDS
                              else ["mdot"])
    Rct = {"A": ct_residuals(fA, rows), "B": ct_residuals(fB, rows)}
    out, types_ct, meta = {}, {}, {}
    for tag, f, can, res in (("A", fA, A, resA), ("B", fB, B, resB)):
        npz = np.load(os.path.join(args.out, f"{tag}_res.npz"))
        tnames = bytes(npz["type_names"]).decode().splitlines()
        tnat = np.array(tnames, dtype=object)[npz["type_codes"] - 1]  # 1-based
        tct, anchor, el, er = ct_types(f, spec, Rct[tag], can)
        out[tag] = (npz["physical"], tnat)
        types_ct[tag] = tct
        meta[tag] = {"anchor_ct": anchor, "excess_left": el, "excess_right": er,
                     "anchor_native": res["anchor"],
                     "dependent_species": res["dependent_species"]}
    scales = (resB if args.scales == "B" else resA)["scales"]
    species_index = {s: k + 1 for k, s in enumerate(species)}
    state = {"A": A, "B": B}
    # Rows whose every term vanishes on the scales state (a species with zero
    # inlet fraction, or a frozen solution) get S ~ 1e-19 that turns transport
    # roundoff into rel ~ 1e6; floor such S at SCALE_FLOOR x the largest scale
    # in the problem so numerically-zero rows stay numerically-zero in rel.
    # The floor can only make rel smaller, so it never hides a real residual.
    SCALE_FLOOR = 1e-10
    global_max = max(scales.values(), default=0.0)

    def scale_for(row, eqtype):
        if eqtype in ("closure", "closure_ct"):
            return scales.get("closure", 1.0)
        if eqtype in ("species", "flux_bc", "zg_Y") and row in species_index:
            S = scales.get(f"{eqtype}/{species_index[row]}", np.inf)
        else:
            S = scales.get(eqtype, np.inf)
        return S if not np.isfinite(S) else max(S, SCALE_FLOOR * global_max)

    groups = {t: (group_rows(out[t][0], out[t][1], rows),
                  group_rows(Rct[t], types_ct[t], rows)) for t in ("A", "B")}
    keys = {k for t in ("A", "B") for g in groups[t] for k in g}

    def worst_rel(key):
        S = scale_for(*key)
        vals = [g[key][0] / S for g in groups.values() if key in g
                and np.isfinite(S) and S > 0]
        return -max(vals, default=0.0)

    table, agree = [], {"A": [], "B": []}
    for key in sorted(keys, key=worst_rel):
        row, eqtype = key
        i = rows.index(row)
        S = scale_for(row, eqtype)
        entry = {"row": row, "eqtype": eqtype, "scale": S, "max": {}, "rel": {},
                 "agree_A": None, "agree_B": None, "j": None, "z": None}
        for t in ("A", "B"):
            Rn, tnat = out[t]
            tct = types_ct[t]
            gn = groups[t][0].get(key)
            gc = groups[t][1].get(key)
            entry["max"][f"{t}_native"] = gn[0] if gn else None
            entry["max"][f"{t}_cantera"] = gc[0] if gc else None
            entry["rel"][f"{t}_native"] = gn[0] / S if gn and S > 0 else None
            entry["rel"][f"{t}_cantera"] = gc[0] / S if gc and S > 0 else None
            if gn or gc:
                j = (gn or gc)[1]
                entry["j"], entry["z"] = j, float(state[t]["grid"][j - 1])
            # agreement over points where both codes solve this equation type
            js = [j for j in range(Rn.shape[1]) if tnat[i, j] == eqtype
                  and tct[i, j] == eqtype]
            if eqtype == "closure":     # native sum(Y)-1 vs cantera 1-sum(Y)
                js += [j for j in range(Rn.shape[1]) if tnat[i, j] == "closure"
                       and tct[i, j] == "closure_ct"]
            if js:
                sign = -1.0 if eqtype == "closure" else 1.0
                d = float(np.max(np.abs(Rn[i, js] - sign * Rct[t][i, js])) / S)
                agree[t].append(d)
                entry[f"agree_{t}"] = d
        table.append(entry)
    structural = set()
    for t in ("A", "B"):
        _, tnat = out[t]
        for i, row in enumerate(rows):
            for j in range(tnat.shape[1]):
                a, b = tnat[i, j], types_ct[t][i, j]
                if a != b and not (a == "closure" and b == "closure_ct"):
                    structural.add(f"{row}@j{j + 1}: native={a} cantera={b}")
    aggregates = {f"rel_max_{t}_{c}":
                  max((e["rel"][f"{t}_{c}"] for e in table
                       if e["rel"][f"{t}_{c}"] is not None), default=None)
                  for t in ("A", "B") for c in ("native", "cantera")}
    aggregates.update({f"agreement_{t}": max(agree[t], default=None)
                       for t in ("A", "B")})

    print("-- residual: max |r| physical; rel = |r|/S (common scales); "
          "agree = |native-CT|/S on shared rows")
    print(f"   {'row':14s} {'eqtype':12s} {'S':>9s} {'|r|n@A':>9s} {'|r|c@A':>9s} "
          f"{'|r|n@B':>9s} {'|r|c@B':>9s} {'relA_n':>8s} {'relA_c':>8s} "
          f"{'relB_n':>8s} {'relB_c':>8s} {'agrA':>8s} {'agrB':>8s} {'j':>5s} {'z':>9s}")
    for e in table[:args.top]:
        mx = " ".join(f"{fmt(e['max'][f'{t}_{c}']):>9s}"
                      for t in ("A", "B") for c in ("native", "cantera"))
        rl = " ".join(f"{fmt(e['rel'][f'{t}_{c}']):>8s}"
                      for t in ("A", "B") for c in ("native", "cantera"))
        ag = " ".join(f"{fmt(e[f'agree_{t}']):>8s}" for t in ("A", "B"))
        j = "-" if e["j"] is None else f"{e['j']:5d}"
        z = "-" if e["z"] is None else f"{e['z']:9.4f}"
        print(f"   {e['row'][:14]:14s} {e['eqtype'][:12]:12s} {fmt(e['scale']):>9s} "
              f"{mx} {rl} {ag} {j} {z}")
    print("-- aggregates")
    for k, v in aggregates.items():
        print(f"   {k:24s}: {fmt(v) if v is not None else 'n/a'}")
    print("-- structural (agreement skipped)")
    for s in sorted(structural)[:20]:
        print(f"   {s}")
    if len(structural) > 20:
        print(f"   ... {len(structural) - 20} more")
    for t in ("A", "B"):
        m = meta[t]
        anc = "-" if m["anchor_ct"] is None else m["anchor_ct"] + 1  # 1-based
        print(f"-- {t}: anchor native={m['anchor_native']} cantera={anc} "
              f"excess L/R={m['excess_left']}/{m['excess_right']} "
              f"dependent={m['dependent_species']}")
    if spec["kind"] in ("burner", "fixed"):
        # Cantera's Inlet1D adds spec-mdot*Y_in to the left species rows
        # (Boundary1D.cpp Inlet1D::eval); native's unstrained row multiplies
        # Y_in by the STATE mass flux (Flames.jl, burner j=1 row). The rows
        # are the same equation exactly when the inlet mass-flux row holds;
        # report the anchor gap so a large left flux_bc disagreement on an
        # arbitrary guess state is attributable, not mysterious.
        for t, f, can in (("A", fA, A), ("B", fB, B)):
            f.gas.TPY = can["T"][0], spec["P_atm"] * ct.one_atm, can["Y"][0]
            gap = abs(f.gas.density * can["velocity"][0] - spec["mdot"])
            print(f"-- {t}: inlet anchor gap |rho*u(0)-mdot| = {gap:.3e} "
                  f"(left species-row offset <= {gap * can['Y'][0].max():.3e})")

    if args.json:
        json.dump({"spec": spec.get("id", spec["kind"]), "kind": spec["kind"],
                   "roundtrip": rt, "scales_source": args.scales, "scales": scales,
                   "rows": table, "aggregates": aggregates,
                   "structural": sorted(structural), "meta": meta},
                  open(args.json, "w"), indent=1)


if __name__ == "__main__":
    main()
