#!/usr/bin/env python
"""Cantera-side dump tool for lockstep comparison with Arrhenius.jl.

Usage: call from any directory; bundled mechanism paths resolve from the repository root:
    ct_dump.py spec.json outdir state    # initial-guess state + steady residual + Jacobian
    ct_dump.py spec.json outdir trace    # full solve; per-timestep states + stdout log

state mode writes outdir/ct_state.npz:
    grid                          (N,) flame-domain grid [m]
    component_names               (C,) flame-domain component names
    values/<name>                 (N,) solution component profiles (u, V, T, Lambda, Y_*)
    residual/<name>               (N,) steady residual per component after eval(rdt=0)
    jacobian                      (n, n) dense steady Jacobian at the same state
    inlet/outlet residuals        boundary-domain residuals (flat)
state mode also writes outdir/ct_guess.csv (grid, T, u, Y_* full precision) for
loading the identical state into native via set_initial_guess! or NPZ.

trace mode writes:
    ct_solve.log                  raw solver stdout at the requested loglevel
    ct_trace.npz                  stacked per-successful-timestep arrays:
        step_grid (S,) grid index per step, step_dt (S,), and per step
        step_<i>/T, step_<i>/u, step_<i>/Y (n_species, N), step_<i>/grid
    debug_sim1d.yaml is left in the CWD by Cantera when loglevel > 6
    (NewtonFail / Timestepping / NewtonSuccess states with residuals).

Environment: CANTERA_LOGLEVEL (default 2 for trace; state mode evaluates quietly).
"""
import json
import os
import sys
from pathlib import Path

for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np
import cantera as ct


def build(spec):
    mechanism = Path(spec["mechanism"])
    if not mechanism.is_absolute():
        mechanism = Path(__file__).resolve().parents[2] / mechanism
    gas = ct.Solution(str(mechanism) + ".yaml")
    kind = spec["kind"]
    if kind in ("diffusion", "premixed", "twin", "stagnation"):
        # Counterflow family; case setup mirrors the harness Cantera driver.
        if kind == "diffusion":
            gas.TP = spec["T"], spec["P_atm"] * ct.one_atm
            f = ct.CounterflowDiffusionFlame(gas, width=spec["width"])
            f.fuel_inlet.mdot, f.fuel_inlet.X, f.fuel_inlet.T = (
                spec["mdot_fuel"], spec["fuel"], spec["T"])
            f.oxidizer_inlet.mdot, f.oxidizer_inlet.X, f.oxidizer_inlet.T = (
                spec["mdot_oxidizer"], spec["oxidizer"], spec.get("T_oxidizer", spec["T"]))
        else:
            gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
            if kind == "premixed":
                f = ct.CounterflowPremixedFlame(gas, width=spec["width"])
                f.reactants.mdot = spec["mdot"]
                gas.equilibrate("HP")
                f.products.X, f.products.T = gas.X, gas.T
                f.products.mdot = spec["mdot_products"]
            elif kind == "twin":
                f = ct.CounterflowTwinPremixedFlame(gas, width=spec["width"])
                f.reactants.mdot = spec["mdot"]
            else:  # stagnation
                f = ct.ImpingingJet(gas, width=spec["width"])
                f.inlet.mdot = spec["mdot"]
                f.surface.T = spec["T_surface"]
                # Mirrors nh3_cantera.build(): the adapter guesses "equil" at
                # build time; without it the domain keeps its flat 300 K
                # default profile (a frozen state, not the validation guess).
                f.set_initial_guess("equil")
        f.transport_model = spec["transport"]
        prune = spec.get("prune", 0.0)
        f.set_refine_criteria(ratio=spec["ratio"], slope=spec["slope"],
                              curve=spec["curve"], prune=min(prune, spec["slope"],
                              spec["curve"]))
        return f
    gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
    if kind == "free":
        f = ct.FreeFlame(gas, width=spec["width"])
    elif kind == "burner":
        f = ct.BurnerFlame(gas, width=spec["width"])
        f.burner.mdot = spec["mdot"]
    else:
        raise ValueError(f"unsupported kind {spec['kind']}")
    f.transport_model = spec["transport"]
    f.set_refine_criteria(ratio=spec["ratio"], slope=spec["slope"],
                          curve=spec["curve"], prune=0.0)
    return f


def state_arrays(f):
    """Per-component solution profiles of the flow domain plus grid."""
    dom = f.flame
    names = [dom.component_name(i) for i in range(dom.n_components)]
    out = {"grid": np.asarray(dom.grid)}
    for name in names:
        out[f"values/{name}"] = np.asarray(dom.values(name))
    return names, out


def write_state(f, outdir, tag, jacobian=True):
    names, out = state_arrays(f)
    # String metadata goes to a JSON sidecar (NPZ.jl cannot parse str arrays).
    meta = {"component_names": names, "species_names": list(f.gas.species_names)}
    # Domain layout of the flat solution vector: for each domain, its
    # n_components, n_points and component names (x is ordered point-major:
    # x[start + j*ncomp + c]).
    doms = list(f.domains)
    meta["domain_names"] = [d.name for d in doms]
    meta["layout"] = {d.name: {
        "n_components": d.n_components, "n_points": d.n_points,
        "components": [d.component_name(i) for i in range(d.n_components)]}
        for d in doms}
    # Steady residual on the identical state
    f.eval(rdt=0.0)
    for name in names:
        out[f"residual/{name}"] = np.asarray(f.flame.residuals(name))
    for dom in doms:
        if dom.name == "flame":
            continue
        for i in range(dom.n_components):
            cname = dom.component_name(i)
            out[f"residual/{dom.name}:{cname}"] = np.asarray(dom.residuals(cname))
    if jacobian:
        # Steady Jacobian at the same state. The default banded solver does not
        # expose its matrix, so evaluate through an Eigen sparse container; the
        # residual-differentiation path (OneDim::evalJacobian) is the same.
        f.linear_solver = ct.EigenSparseDirectJacobian()
        f.eval_jacobian()
        J = f.linear_solver.jacobian
        out["jacobian"] = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    os.makedirs(outdir, exist_ok=True)
    np.savez(os.path.join(outdir, f"{tag}.npz"), **out)
    with open(os.path.join(outdir, f"{tag}_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    return names, out, meta


def dump_state(spec, outdir):
    f = build(spec)
    # The auto solve path sets the initial guess itself (Sim1D.solve with
    # auto=True -> set_initial_guess on the initial grid); mirror that here.
    if spec["kind"] != "stagnation":
        # build() already set the stagnation initial guess ("equil"), matching
        # the harness adapter; re-guessing would only rebuild the same state.
        f.set_initial_guess()
    # FlowBase.__init__ enables the energy equation, but Flow1D::setupGrid ->
    # resize() extends m_do_energy from the 1-point constructor grid, so all
    # points past the first come back DISABLED. FreeFlame only recovers because
    # its guess calls setFixedTemperature -> Sim1D::finalize ->
    # Flow1D::_finalize -> solveEnergyEqn(); the counterflow guesses never do.
    # energy_enabled reads point 0 and keeps reporting True, and a state dump
    # would then hold rsd[T,j] = T - T_fixed (garbage) at interior points.
    # solve(auto=True) turns the energy equation on before its first steady
    # attempt; mirror that here so the dumped residual is the one CT solves.
    f.energy_enabled = True
    names, out, _ = write_state(f, outdir, "ct_state",
                                jacobian=spec["kind"] in ("free", "burner"))
    # Full-precision CSV for native set_initial_guess!
    hdr = ["grid", "temperature_K", "velocity_m_s"] + [
        f"Y_{s}" for s in f.gas.species_names]
    cols = [out["grid"], out["values/T"], out["values/velocity"]]
    cols += [out[f"values/{s}"] for s in f.gas.species_names]
    with open(os.path.join(outdir, "ct_guess.csv"), "w") as fh:
        fh.write(",".join(hdr) + "\n")
        for row in zip(*cols):
            fh.write(",".join(f"{v:.17g}" for v in row) + "\n")
    print(f"ct_dump state: {len(out['grid'])} points, {len(names)} components")


def dump_solved(spec, outdir):
    f = build(spec)
    try:
        f.solve(loglevel=0, auto=True)
        status = "ok"
    except Exception as e:
        status = str(e).split("\n")[0][:80]
    names, out, _ = write_state(f, outdir, "ct_solved", jacobian=False)
    print(f"ct_dump solved: {status} {len(out['grid'])} points, "
          f"{len(names)} components")


def dump_trace(spec, outdir, loglevel):
    f = build(spec)
    os.makedirs(outdir, exist_ok=True)
    steps = []
    state = {"n": 0}

    def cb(dt):
        i = state["n"]
        state["n"] += 1
        steps.append((i, float(dt), np.asarray(f.flame.grid).copy(),
                      np.asarray(f.flame.values("T")).copy(),
                      np.asarray(f.flame.values("velocity")).copy(),
                      np.asarray(f.Y).copy()))
        return 0.0

    f.set_time_step_callback(cb)
    # Capture C++ writelog output (fd 1) into a file. std::cout flushes go to
    # whatever fd 1 points at at flush time, so keep the redirect until exit
    # and report on the saved original fd.
    logpath = os.path.join(outdir, "ct_solve.log")
    logfd = os.open(logpath, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    old = os.dup(1)
    os.dup2(logfd, 1)
    err = "none"
    try:
        if loglevel >= 4:
            # Verbose logs are only affordable on the initial grid: solve with
            # auto=False (no refinement stages), which still runs the full
            # steady-attempt + time-stepping loop on the first grid.
            f.solve(loglevel=loglevel, refine_grid=False, auto=False)
        else:
            f.solve(loglevel=loglevel, auto=True)
    except Exception as e:
        err = str(e).split("\n")[0][:120]
    stats = f.solver_stats
    njac = sum(stats.get("jacobian_evals", []))
    out = {"converged": np.array(err == "none"), "error": np.array(err),
           "n_steps": np.array(len(steps)),
           "step_index": np.array([s[0] for s in steps]),
           "step_dt": np.array([s[1] for s in steps])}
    with open(os.path.join(outdir, "ct_trace_meta.json"), "w") as fh:
        json.dump({"species_names": list(f.gas.species_names)}, fh)
    for i, dt, grid, T, u, Y in steps:
        out[f"step_{i}/grid"] = grid
        out[f"step_{i}/T"] = T
        out[f"step_{i}/u"] = u
        out[f"step_{i}/Y"] = Y
    np.savez(os.path.join(outdir, "ct_trace.npz"), **out)
    with open(os.path.join(outdir, "ct_trace_meta.json"), "w") as fh:
        json.dump({"species_names": list(f.gas.species_names)}, fh)
    os.write(old, f"ct_dump trace: converged={err == 'none'} steps={len(steps)} "
                 f"jac={njac} log={os.path.getsize(logpath)} bytes\n".encode())
    os.close(old)


def main():
    spec = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else "state"
    if mode == "state":
        dump_state(spec, outdir)
    elif mode == "solved":
        dump_solved(spec, outdir)
    elif mode == "trace":
        dump_trace(spec, outdir, int(os.environ.get("CANTERA_LOGLEVEL", "2")))
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
