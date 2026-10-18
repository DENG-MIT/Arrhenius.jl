#!/usr/bin/env python3
import json
import sys
import numpy as np
import cantera as ct

def capture(spec_path, out):
    spec = json.load(open(spec_path, encoding="utf-8"))
    gas = ct.Solution(spec["mechanism"] + ".yaml")
    gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
    f = ct.BurnerFlame(gas, width=spec["width"])
    f.burner.mdot = spec["mdot"]
    f.transport_model = spec["transport"]
    if spec["kind"] == "fixed":
        f.set_initial_guess()
        f.flame.set_fixed_temp_profile(np.asarray(spec["profile_z"])/spec["width"], spec["profile_T"])
        f.energy_enabled = False
    elif spec["kind"] != "burner":
        raise ValueError("unsupported kind " + spec["kind"])
    captured = []
    def on_steady(_):
        if not captured:
            captured.append((np.asarray(f.grid).copy(), np.asarray(f.T).copy(),
                np.asarray(f.Y).copy(), np.asarray(f.velocity).copy()))
        return 0.0
    f.set_steady_callback(on_steady)
    f.solve(loglevel=0, auto=(spec["kind"] == "burner"), refine_grid=False)
    if not captured:
        raise RuntimeError("no successful steady callback for " + spec["id"])
    grid, temperature, mass_fractions, velocity = captured[0]
    names = "\n".join(gas.species_names).encode("utf-8")
    np.savez_compressed(out, grid=grid, T=temperature, Y=mass_fractions, velocity=velocity,
        species_names_utf8=np.frombuffer(names, dtype=np.uint8))
    print("CT_FIRST_STEADY case={} points={}".format(spec["id"], len(grid)), flush=True)

args=sys.argv[1:]
if len(args)<2 or len(args)%2:
    raise SystemExit("usage: round9_first_ct.py spec.json output.npz [pairs ...]")
for i in range(0,len(args),2):
    capture(args[i],args[i+1])
