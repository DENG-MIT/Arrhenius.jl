#!/usr/bin/env python3
import csv, json, os, sys
import cantera as ct
import numpy as np
source, destination = sys.argv[1:3]
phase = ct.Solution(source, "ar-plasma", transport_model=None)
electron_x = os.environ.get("PLASMA_ELECTRON_X", "1e-30")
phase.TPX = 300.0, 133.32, f"Ar:1.0,Electron:{electron_x}"
phase.reduced_electric_field = 100e-21
phase.update_electron_energy_distribution()
energy = np.asarray(phase.electron_energy_levels, dtype=float)
values = np.asarray(phase.electron_energy_distribution, dtype=float)
if len(energy) < 3 or energy[0] != 0.0 or energy[-1] != 100.0:
    raise RuntimeError(f"unexpected fixed grid: n={len(energy)} endpoints={energy[[0,-1]]}")
mobility = float(phase.electron_mobility)
if not np.all(np.isfinite(values)) or np.any(values <= 0.0) or not np.isfinite(mobility) or mobility <= 0.0:
    raise RuntimeError("Cantera returned a non-positive or non-finite EEDF / mobility")
normalization = float(np.trapezoid(values * np.sqrt(energy), energy))
with open(destination, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["energy_eV", "eedf"]); w.writerows(zip(energy, values))
with open(destination.replace(".csv", "-metadata.json"), "w") as f:
    json.dump({"cantera_mobility_m2_V_s": mobility, "normalization": normalization,
               "version": ct.__version__, "points": len(energy)}, f, indent=2)
print(f"CANTERA_OK version={ct.__version__} points={len(energy)} electron_X={electron_x} mobility={mobility:.17g} normalization={normalization:.17g}")

