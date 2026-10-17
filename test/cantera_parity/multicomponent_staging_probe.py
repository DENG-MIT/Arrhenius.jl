"""Probe Cantera FreeFlame auto-stage transport and final flame speed."""
import argparse
from pathlib import Path

import cantera as ct
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--soret", choices=("on", "off"), required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    gas = ct.Solution(str(root / "mechanism" / "h2o2.yaml"))
    gas.TPX = 300.0, ct.one_atm, "H2:1.1, O2:1, AR:5"
    flame = ct.FreeFlame(gas, width=0.03)
    flame.flame.grid = 0.03 * np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    flame.transport_model = "multicomponent"
    flame.soret_enabled = args.soret == "on"
    flame.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.10, prune=0.05)

    transport_sequence = []
    print(f"initial_transport={flame.transport_model}", flush=True)

    def observe(_):
        model = flame.transport_model
        if not transport_sequence or transport_sequence[-1] != model:
            transport_sequence.append(model)
            print(f"steady_callback_transport={model}", flush=True)
        return 0.0

    flame.set_steady_callback(observe)
    flame.solve(loglevel=0, auto=True, refine_grid=True)
    print(f"final_transport={flame.transport_model}", flush=True)
    print(
        f"soret={args.soret} Su_cm_s={100 * flame.velocity[0]:.14g} "
        f"points={len(flame.grid)} width={flame.grid[-1] - flame.grid[0]:.14g} "
        f"transport_sequence={' -> '.join(transport_sequence)}"
    )


if __name__ == "__main__":
    main()
