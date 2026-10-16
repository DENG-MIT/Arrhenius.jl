"""Complete canonical ionized flame sequences for independent validation."""

import time

import numpy as np
import cantera as ct

REACTANTS = "CH4:1,O2:2,N2:7.52"
WIDTH = 0.05  # m
ION_SPECIES = ("E", "H3O+", "HCO+")


def _validate_case(case):
    if case not in ("free", "burner"):
        raise ValueError("case must be free or burner")


def _validate_mechanism(gas):
    # Outside the timer: pin the mechanism to the ionized gri30 variant.
    if gas.n_species != 56 or gas.n_reactions != 331:
        raise ValueError("supply the source gri30_ion mechanism")
    if gas.transport_model != "ionized-gas":
        raise ValueError("supply ionized-gas transport")
    names = set(gas.species_names)
    missing = [s for s in ION_SPECIES if s not in names]
    if missing:
        raise ValueError("mechanism lacks ionized species: " + ", ".join(missing))


def run_ion_source_sequence(gas, case, *, after_stage=None,
                            clock_ns=time.monotonic_ns, loglevel=0):
    """Run both ionized stages; return seconds, points, snapshots.

    Stage "frozen": electric field off, canonical ``solve(auto=True)``.
    Stage "field":  electric field on, canonical ``solve`` without ``auto``.

    The gas mechanism object is reused; a fresh flame is built per call
    and is never reset between stages. Any solve failure propagates --
    no retries, no fallbacks, no tolerance changes.
    """
    _validate_case(case)
    _validate_mechanism(gas)
    species_names = tuple(gas.species_names)  # pinned outside the timer
    Tin = 300.0 if case == "free" else 600.0

    seconds = [0.0, 0.0]
    points = [0, 0]
    snapshots = {}
    f = None
    for i, stage in enumerate(("frozen", "field")):
        started = clock_ns()
        if i == 0:
            gas.TPX = Tin, ct.one_atm, REACTANTS
            # Capture inlet mass fractions immediately after the canonical
            # reset, before any property query can disturb the phase state.
            inlet_Y = np.array(gas.Y, copy=True)
            if case == "free":
                f = ct.FreeFlame(gas, width=WIDTH)
            else:
                mdot = 0.15 * gas.density
                f = ct.BurnerFlame(gas, width=WIDTH)
                f.burner.mdot = mdot
            f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
        f.electric_field_enabled = (i == 1)
        if i == 0:
            f.solve(loglevel=loglevel, auto=True)
        else:
            f.solve(loglevel=loglevel)
        snapshot = {
            "grid": np.array(f.grid, copy=True),
            "T": np.array(f.T, copy=True),
            # Raw mass fractions straight from the flow domain.
            "Y": np.array([f.flame.values(name) for name in species_names],
                          copy=True),
            "X": np.array(f.X, copy=True),
            "E": np.array(f.E, copy=True),  # no silent zero fallback
            "velocity": np.array(f.velocity, copy=True),
            "rho": np.array(f.density, copy=True),
            "qdot": np.array(f.heat_release_rate, copy=True),
            "inlet_Y": inlet_Y.copy(),
            "P": np.array([f.P]),
        }
        seconds[i] = (clock_ns() - started) / 1e9
        # Timer stops before validation and the after_stage callback.
        snapshots[stage] = snapshot
        for name, arr in snapshot.items():
            if not np.all(np.isfinite(arr)):
                raise RuntimeError(
                    f"nonfinite ionized flame output: {case} {stage} {name}")
        points[i] = len(f.grid)
        if after_stage is not None:
            after_stage(f, stage)
    return seconds, points, snapshots
