"""Generate Cantera reference initial grids/guesses for the native counterflow
initial-guess regression tests (test/counterflow_initial_guess_native.jl).

Reads the repository's mechanism/h2o2.yaml so both codes see the identical
mechanism, builds each counterflow variant with the same inputs as the test,
and stores grid, T, velocity, spreadRate, Lambda and Y.

Usage: generate_counterflow_init_refs.py <repo_root>
"""
import sys
from pathlib import Path
import numpy as np
import cantera as ct


def save(path, f, gas):
    np.savez(path, grid=f.grid.copy(), T=f.T.copy(), velocity=f.velocity.copy(),
             spreadRate=f.flame.spread_rate.copy(), Lambda=f.flame.Lambda.copy(),
             Y=f.Y.copy(),
             species_names=np.frombuffer("\n".join(gas.species_names).encode(),
                                         dtype=np.uint8))


root = Path(sys.argv[1])
yaml = str(root / "mechanism" / "h2o2.yaml")
out = root / "test" / "data"

gas = ct.Solution(yaml)
gas.TP = 300.0, ct.one_atm
f = ct.CounterflowDiffusionFlame(gas, width=0.02)
f.fuel_inlet.mdot, f.fuel_inlet.X, f.fuel_inlet.T = 0.24, "H2:1,AR:1", 300.0
f.oxidizer_inlet.mdot, f.oxidizer_inlet.X, f.oxidizer_inlet.T = 0.72, "O2:.2,AR:.8", 300.0
f.set_initial_guess()
save(out / "cantera_counterflow_diffusion_init.npz", f, gas)

composition = "H2:1.6,O2:1,AR:7"
P = 0.05 * ct.one_atm
gas.TPX = 373.0, P, composition
f = ct.CounterflowPremixedFlame(gas, width=0.2)
f.reactants.mdot = 0.12
gas.equilibrate("HP")
f.products.X, f.products.T = gas.X, gas.T
f.products.mdot = 0.06
gas.TPX = 373.0, P, composition
f.set_initial_guess()
save(out / "cantera_counterflow_premixed_init.npz", f, gas)

gas.TPX = 373.0, P, composition
f = ct.CounterflowTwinPremixedFlame(gas, width=0.2)
f.reactants.mdot = 0.12
f.set_initial_guess()
save(out / "cantera_counterflow_twin_init.npz", f, gas)

gas.TPX = 373.0, P, composition
f = ct.ImpingingJet(gas, width=0.2)
f.inlet.mdot = 0.06
f.surface.T = 500.0
f.set_initial_guess("equil")
save(out / "cantera_counterflow_stagnation_init.npz", f, gas)

print("wrote 4 reference files to", out)
