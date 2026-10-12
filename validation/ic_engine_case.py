"""Exact pinned Cantera 4 ic_engine.py calculation and native mechanism preparation.

https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/reactors/ic_engine.py
Keeps the source's 8-revolution loop, one-degree maximum output interval,
20 K advance limit, physical parameters, and integration tolerances.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import cantera as ct
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output",type=Path)
    parser.add_argument("--accepted-states",action="store_true",
        help="record every accepted step to check source output-quadrature convergence")
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if not ct.__version__.startswith("4."):
        raise RuntimeError("This reference requires Cantera 4")
    phase = ct.Solution("nDodecane_Reitz.yaml","nDodecane_IG")
    prepared = args.output / "dodecane_IG.yaml"
    phase.write_yaml(str(prepared))
    exporter_path = Path(__file__).resolve().parents[1] / "mechanism" / "export_sidecar.py"
    spec = importlib.util.spec_from_file_location("sidecar_exporter",exporter_path)
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    exporter.export(prepared,Path(str(prepared)+".npz"))
    frequency,displacement,compression,diameter = 50.0,0.5e-3,20.0,0.083
    clearance = displacement/(compression-1)
    area = np.pi*diameter**2/4
    stroke = displacement/area
    def gas(T,P,X):
        g = ct.Solution("nDodecane_Reitz.yaml","nDodecane_IG")
        g.TPX = T,P,X
        return g
    cylinder = ct.IdealGasReactor(gas(300,1.3e5,"o2:1,n2:3.76"),volume=clearance,clone=True)
    inlet = ct.Reservoir(gas(300,1.3e5,"o2:1,n2:3.76"),clone=True)
    injector = ct.Reservoir(gas(300,1600e5,"c12h26:1"),clone=True)
    outlet = ct.Reservoir(gas(300,1.2e5,"o2:1,n2:3.76"),clone=True)
    ambient = ct.Reservoir(gas(300,1e5,"o2:1,n2:3.76"),clone=True)
    crank = lambda t:np.remainder(2*np.pi*frequency*t,4*np.pi)
    gate = lambda t,opening,closing:np.mod(crank(t)-np.deg2rad(opening),4*np.pi)<np.mod(np.deg2rad(closing-opening),4*np.pi)
    inlet_valve = ct.Valve(inlet,cylinder,K=1e-6)
    inlet_valve.time_function = lambda t:gate(t,-18,198)
    injection_rate = 3.2e-5/((365-350)/360/frequency)
    injector_mfc = ct.MassFlowController(injector,cylinder,mdot=lambda t:injection_rate*gate(t,350,365))
    outlet_valve = ct.Valve(cylinder,outlet,K=1e-6)
    outlet_valve.time_function = lambda t:gate(t,522,18)
    piston_speed = lambda t:-stroke*np.pi*frequency*np.sin(crank(t))
    wall = ct.Wall(ambient,cylinder,A=area,velocity=piston_speed)
    net = ct.ReactorNet([cylinder])
    net.rtol,net.atol = 1e-12,1e-16
    net.max_steps = 1000000
    if args.accepted_states:
        net.max_time_step = 1/(360*frequency)
    cylinder.set_advance_limit("temperature",20.0)
    keys = ("time","crank_angle","temperature","pressure","volume","mass","mean_molecular_weight",
            "entropy_mass","internal_energy","mdot_in","mdot_fuel","mdot_out","work_rate","heat_release_rate")
    rows,X,Y = [],[],[]
    def record():
        g = cylinder.phase
        t = net.time
        rows.append([t,crank(t),g.T,g.P,cylinder.volume,cylinder.mass,g.mean_molecular_weight,
            g.entropy_mass,cylinder.mass*g.int_energy_mass,inlet_valve.mass_flow_rate,
            injector_mfc.mass_flow_rate,outlet_valve.mass_flow_rate,
            -(g.P-1e5)*area*piston_speed(t),g.heat_release_rate*cylinder.volume])
        X.append(g.X.copy())
        Y.append(g.Y.copy())
    net.initialize()
    record()
    while net.time < 8/frequency:
        if args.accepted_states:
            if 8/frequency-net.time <= 1/(360*frequency):
                net.advance(8/frequency,apply_limit=False)
            else:
                net.step()
        else:
            net.advance(net.time+1/(360*frequency))
        record()
    arrays = {key:values for key,values in zip(keys,np.asarray(rows).T)}
    arrays.update(X=np.asarray(X),Y=np.asarray(Y))
    np.savez(args.output/"ic_engine.reference.npz",**arrays)
    trapezoid = getattr(np,"trapezoid",None) or np.trapz
    first = 0 if args.accepted_states else 1
    source_times = arrays["time"][first:]
    heat = trapezoid(arrays["heat_release_rate"][first:],source_times)
    work = trapezoid(arrays["work_rate"][first:],source_times)
    weights = arrays["mean_molecular_weight"][first:]*arrays["mdot_out"][first:]
    emission = trapezoid(weights*arrays["X"][first:,phase.species_index("co")],source_times)/trapezoid(weights,source_times)
    metadata = dict(cantera=ct.__version__,source_commit="726522be4e2a13454d8415b7ef799d621f665cf3",
        points=len(rows),species=phase.n_species,reactions=phase.n_reactions,
        rtol=1e-12,atol=1e-16,temperature_advance_limit_K=20.0,
        heat_J=float(heat),work_J=float(work),efficiency=float(work/heat),CO_ppm=float(1e6*emission),
        end_time=float(net.time),source_output_initial_state=False)
    metadata["output_sampling"] = "accepted_states" if args.accepted_states else "source_advance_limit"
    if args.accepted_states:
        selected = np.unique(np.r_[np.arange(0,len(rows),2),len(rows)-1])
        metadata["every_second_state_heat_J"] = float(trapezoid(arrays["heat_release_rate"][selected],arrays["time"][selected]))
    (args.output/"manifest.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(metadata,flush=True)


if __name__ == "__main__":
    main()
