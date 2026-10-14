#!/usr/bin/env python
"""Lockstep comparison of Cantera vs Arrhenius.jl dumps (test/cantera_parity).

Usage:
    compare_lock.py DIR resid       Level 1: residual on the same state
    compare_lock.py DIR step        Level 2: one Newton step from the same state
    compare_lock.py DIR invariants  Level 3: per-step invariants of both traces

DIR must contain ct_state.npz + ct_meta.json (ct_dump.py state) and
native_state_cons.npz / native_state_fd.npz (native_dump.jl state);
for invariants, ct_trace.npz/ct_trace_meta.json and native_trace.npz.

Native rows are [T/1000; Y_1..Y_n; mdot]; Cantera flame components are
[velocity, spreadRate, T, Lambda, eField, Uo, Y_1..Y_n]. Native residuals are
scaled by _flame_timescale=1e-4 (species, per-1/(rho cp) for T), so residual
agreement is judged by the per-point ratio native/Cantera being CONSTANT, and
the Newton step (scale-invariant) by relative difference.
"""
import json
import os
import sys

import numpy as np


def load(dir):
    ct = np.load(os.path.join(dir, "ct_state.npz"))
    meta = json.load(open(os.path.join(dir, "ct_meta.json")))
    nat = {}
    for tag in ("cons", "fd"):
        p = os.path.join(dir, f"native_state_{tag}.npz")
        if os.path.exists(p):
            nat[tag] = np.load(p)
    return ct, meta, nat


def interior(grid):
    return slice(1, len(grid) - 1)


def report_ratio(name, ct_res, nat_res, grid):
    """Report how close native_res/ct_res is to a constant per point."""
    n = len(ct_res)
    sel = np.abs(ct_res) > 1e-3 * max(np.max(np.abs(ct_res)), 1e-300)
    sel[0] = sel[-1] = False  # boundary rows hold different equations
    if not np.any(sel):
        print(f"  {name:8s} max|ct|={np.max(np.abs(ct_res)):9.3e} "
              f"max|nat|={np.max(np.abs(nat_res)):9.3e}  (no interior signal)")
        return
    ratio = nat_res[sel] / ct_res[sel]
    med = np.median(ratio)
    spread = np.max(np.abs(ratio / med - 1))
    imax = int(np.argmax(np.abs(ratio / med - 1)))
    idx = np.flatnonzero(sel)[imax]
    print(f"  {name:8s} max|ct|={np.max(np.abs(ct_res)):9.3e} "
          f"max|nat|={np.max(np.abs(nat_res)):9.3e} ratio med={med:9.4e} "
          f"spread={spread:9.2e} at z={grid[idx]:.4g}")


def cmd_resid(dir):
    ct, meta, nat = load(dir)
    grid = ct["grid"]
    species = meta["species_names"]
    for tag, nd in nat.items():
        print(f"== Level 1 residual, native :{ 'conservative' if tag=='cons' else 'finite_difference'} ==")
        nr = nd["residual"]  # rows [T; Y_k; M], cols points
        # state agreement first
        dT = np.max(np.abs(nd["state"][0, :] * 1000 - ct["values/T"]))
        print(f"  state: max|dT|={dT:.3e} K  anchor={int(nd['anchor'][0])}")
        report_ratio("T", ct["residual/T"], nr[0, :], grid)
        for k, s in enumerate(species):
            report_ratio(s, ct[f"residual/{s}"], nr[k + 1, :], grid)
        report_ratio("M/u", ct["residual/velocity"], nr[-1, :], grid)


def flat_ct_residual(ct, meta):
    """Assemble Cantera's flat residual: domain order, point-major."""
    parts = []
    for d in meta["domain_names"]:
        lay = meta["layout"][d]
        nc, np_ = lay["n_components"], lay["n_points"]
        if d == "flame":
            for j in range(np_):
                for c in lay["components"]:
                    parts.append(ct[f"residual/{c}"][j])
        else:
            for j in range(np_):
                for c in lay["components"]:
                    parts.append(ct[f"residual/{d}:{c}"][j])
    return np.array(parts)


def flame_step_components(step, meta):
    lay = meta["layout"]["flame"]
    nc, np_ = lay["n_components"], lay["n_points"]
    start = sum(meta["layout"][d]["n_components"] * meta["layout"][d]["n_points"]
                for d in meta["domain_names"][:meta["domain_names"].index("flame")])
    block = step[start:start + nc * np_].reshape(np_, nc).T
    return {c: block[i] for i, c in enumerate(lay["components"])}


def cmd_step(dir):
    ct, meta, nat = load(dir)
    species = meta["species_names"]
    F = flat_ct_residual(ct, meta)
    J = ct["jacobian"]
    step = np.linalg.solve(J, -F)
    comps = flame_step_components(step, meta)
    grid = ct["grid"]
    print(f"== Level 2 Newton step from shared state (n={len(F)}) ==")
    print(f"  CT: max|step T|={np.max(np.abs(comps['T'])):.4e} K  "
          f"max|step u|={np.max(np.abs(comps['velocity'])):.4e}")
    for tag, nd in nat.items():
        ns = nd["newton_step"]  # rows [T/1000; Y; M]
        dT = ns[0, :] * 1000 - comps["T"]
        relT = np.max(np.abs(dT)) / np.max(np.abs(comps["T"]))
        worst = ("", 0.0)
        for k, s in enumerate(species):
            ref = np.max(np.abs(comps[s]))
            rel = np.max(np.abs(ns[k + 1, :] - comps[s])) / max(ref, 1e-300)
            if rel > worst[1]:
                worst = (s, rel)
        print(f"  native:{tag:5s} rel step diff T={relT:.3e}  worst species "
              f"{worst[0]}={worst[1]:.3e}  fbound={nd['fbound'][0]:.4e}")
        # where is the worst T mismatch?
        j = int(np.argmax(np.abs(dT)))
        print(f"    max dT at z={grid[j]:.4g}: native={ns[0, j] * 1000:.6e} "
              f"ct={comps['T'][j]:.6e}")


def invariants(grid, T, Y, mdot, species, gas_MW, elements):
    """Element balance H/N/O error, sumY-1, mass-flux continuity."""
    out = {}
    out["sumY-1"] = float(np.max(np.abs(Y.sum(axis=0) - 1)))
    out["mdot_cont"] = float(np.max(np.abs(np.diff(mdot)))) if len(mdot) > 1 else 0.0
    for el in ("H", "N", "O"):
        col = elements[el]
        per = Y * col[:, None]
        inlet = per[:, 0].sum() * 1.0
        err = np.max(np.abs(per.sum(axis=0) - inlet)) / max(abs(inlet), 1e-300)
        out[f"elem_{el}"] = float(err)
    return out


def cmd_invariants(dir):
    """Per-time-step invariants for both traces; needs Cantera for elements."""
    import cantera as ctmod
    spec = json.load(open(os.environ["LOCKSTEP_SPEC"]))
    gas = ctmod.Solution(spec["mechanism"] + ".yaml")
    gas.TPX = spec["T"], spec["P_atm"] * ctmod.one_atm, spec["X"]
    species = list(gas.species_names)
    MW = gas.molecular_weights
    eix = {el: gas.element_index(el) for el in ("H", "N", "O")}
    # elemental mass content per species: atoms * MW_el / MW_species
    elem = {el: np.array([gas.n_atoms(k, eix[el]) * gas.atomic_weight(eix[el])
                          / MW[k] for k in range(gas.n_species)])
            for el in ("H", "N", "O")}
    R = 8314.46261815324

    def inv(tag, grid, T, Y, mdot):
        inlet = {el: float(np.dot(elem[el], Y[:, 0])) for el in elem}
        out = {"sumY-1": float(np.max(np.abs(Y.sum(axis=0) - 1))),
               "mdot_cont": float(np.max(np.abs(np.diff(mdot))))}
        for el in elem:
            prof = elem[el] @ Y
            out[f"elem_{el}"] = float(np.max(np.abs(prof - inlet[el])) /
                                      max(abs(inlet[el]), 1e-300))
        return out

    ct = np.load(os.path.join(dir, "ct_trace.npz"))
    nat = np.load(os.path.join(dir, "native_trace.npz"))
    nct, nnat = int(np.ravel(ct["n_steps"])[0]), int(np.ravel(nat["n_steps"])[0])
    print(f"steps: ct={nct} native={nnat}")
    hdr = f"{'i':>3s} {'side':>5s} {'dt':>10s} {'sumY-1':>9s} {'H':>9s} {'N':>9s} {'O':>9s} {'mdot_cont':>9s}"
    print(hdr)
    rows_ct, rows_nat = [], []
    for i in range(nct):
        grid = ct[f"step_{i}/grid"]; T = ct[f"step_{i}/T"]
        Y = ct[f"step_{i}/Y"]; u = ct[f"step_{i}/u"]
        MWmix = 1.0 / (Y.T @ (1.0 / MW))
        rho = spec["P_atm"] * ctmod.one_atm * MWmix / (R * T)
        rows_ct.append((i, ct["step_dt"][i], inv("ct", grid, T, Y, rho * u)))
    for i in range(nnat):
        st = nat[f"step_{i}/state"]
        Y = st[1:-1, :]
        rows_nat.append((i, nat["step_dt"][i], inv("nat", None, None, Y, st[-1, :])))
    for i in range(min(nct, nnat, 25)):
        _, dtc, rc = rows_ct[i]
        _, dtn, rn = rows_nat[i]
        print(f"{i:3d} {'ct':>5s} {dtc:10.3e} {rc['sumY-1']:9.2e} {rc['elem_H']:9.2e} {rc['elem_N']:9.2e} {rc['elem_O']:9.2e} {rc['mdot_cont']:9.2e}")
        print(f"{i:3d} {'nat':>5s} {dtn:10.3e} {rn['sumY-1']:9.2e} {rn['elem_H']:9.2e} {rn['elem_N']:9.2e} {rn['elem_O']:9.2e} {rn['mdot_cont']:9.2e}")
        # state difference while grids match
        gc, gn = ct[f"step_{i}/grid"], nat[f"step_{i}/grid"]
        if len(gc) == len(gn) and np.allclose(gc, gn):
            Tc, Tn = ct[f"step_{i}/T"], nat[f"step_{i}/state"][0, :] * 1000
            Yc, Yn = ct[f"step_{i}/Y"], nat[f"step_{i}/state"][1:-1, :]
            dT = np.max(np.abs(Tn - Tc)) / np.max(np.abs(Tc))
            dY = np.max(np.abs(Yn - Yc) / (np.abs(Yc) + 1e-30))
            print(f"{'':3s} diff  relT={dT:9.2e} relY(max over species)={dY:9.2e}")


def main():
    dir = sys.argv[1]
    cmd = sys.argv[2] if len(sys.argv) > 2 else "resid"
    {"resid": cmd_resid, "step": cmd_step, "invariants": cmd_invariants}[cmd](dir)


if __name__ == "__main__":
    main()
