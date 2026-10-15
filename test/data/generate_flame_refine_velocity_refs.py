"""Generate one-pass Cantera 4.0 refinement fixtures for planar flames."""
from pathlib import Path
import sys
import cantera as ct
import numpy as np
repo = Path(__file__).resolve().parents[2]
gas_path = repo / "mechanism" / "h2o2.yaml"
records = {}
for kind in ("free", "burner", "fixedT"):
    gas = ct.Solution(str(gas_path))
    gas.TPX = 300.0, ct.one_atm, "H2:1.6,O2:1,AR:7"
    if kind == "free":
        f = ct.FreeFlame(gas, width=0.03)
    else:
        f = ct.BurnerFlame(gas, width=0.03)
        f.burner.mdot = 0.04
    f.transport_model = "mixture-averaged"
    f.set_initial_guess()
    if kind == "fixedT":
        profile_x = np.array([0.0, 0.25, 0.6, 1.0])
        profile_T = np.array([300.0, 420.0, 670.0, 900.0])
        f.flame.set_fixed_temp_profile(profile_x, profile_T)
        f.energy_enabled = False
        f.flame.set_profile("T", profile_x, profile_T)
    f.set_refine_criteria(ratio=3.0, slope=0.03, curve=0.06, prune=0.0)
    pre = {"grid": np.asarray(f.grid), "T": np.asarray(f.T), "Y": np.asarray(f.Y),
           "velocity": np.asarray(f.velocity)}
    if kind == "free":
        pre["fixed_temperature"] = np.asarray(f.fixed_temperature)
        pre["fixed_point"] = np.asarray(f.fixed_temperature_location)
    for i in range(f.flame.n_components):
        name = f.flame.component_name(i)
        pre[f"component/{name}"] = np.asarray(f.flame.values(name))
    f.refine(loglevel=0)
    post = {"grid": np.asarray(f.grid), "T": np.asarray(f.T), "Y": np.asarray(f.Y),
            "velocity": np.asarray(f.velocity)}
    for key, value in pre.items():
        records[f"{kind}/pre/{key}"] = value
    for key, value in post.items():
        records[f"{kind}/post/{key}"] = value
    for i in range(f.flame.n_components):
        name = f.flame.component_name(i)
        records[f"{kind}/post/component/{name}"] = np.asarray(f.flame.values(name))
    f.eval(rdt=0.0)
    for i in range(f.flame.n_components):
        name = f.flame.component_name(i)
        if name in ("velocity", "T", *f.gas.species_names):
            records[f"{kind}/post/residual/{name}"] = np.asarray(f.flame.residuals(name))
    print(f"{kind}: {len(pre['grid'])}->{len(post['grid'])} points")
np.savez(sys.argv[1], **records)