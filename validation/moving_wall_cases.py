"""Pinned Cantera 4 moving-wall examples, with exact physical inputs and output times.

Sources: https://github.com/Cantera/cantera/tree/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/reactors
reactor2.py, piston.py, custom2.py. Tolerances are tightened to make trajectory
comparison resolve Julia errors rather than reference integration error.
Cantera is used only for preparation and independent validation.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import shutil

import cantera as ct
import numpy as np


def prepare(output):
    directory = output / "mechanisms"
    directory.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).resolve().parents[1] / "mechanism" / "export_sidecar.py"
    spec = importlib.util.spec_from_file_location("sidecar_exporter", source)
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    for name in ("air.yaml", "gri30.yaml", "h2o2.yaml"):
        source = next(Path(d) / name for d in ct.get_data_directories() if (Path(d)/name).is_file())
        target = directory / name
        shutil.copyfile(source, target)
        exporter.export(target, Path(str(target)+".npz"))
    return directory


class InertialWallReactor(ct.ExtensibleIdealGasReactor):
    """The pinned custom2.py extension, with the same acceleration coefficient."""
    def __init__(self, *args, neighbor, **kwargs):
        super().__init__(*args, **kwargs)
        self.v_wall = 0.0
        self.neighbor = neighbor
        self.n_vars += 1
        self.i_wall = self.n_vars - 1

    def after_get_state(self, y):
        y[self.i_wall] = self.v_wall

    def after_update_state(self, y):
        self.v_wall = y[self.i_wall]
        self.walls[0].velocity = self.v_wall

    def after_eval(self, t, LHS, RHS):
        RHS[self.i_wall] = 0.01*(self.phase.P-self.neighbor.phase.P)


def build(name, directory):
    def gas(file, T, P, X):
        phase = ct.Solution(str(directory / (file+".yaml")))
        phase.TPX = T, P, X
        return phase
    if name == "reactor2":
        argon = ct.IdealGasReactor(gas("air",1000,20*ct.one_atm,"AR:1"),clone=True)
        reacting = ct.IdealGasReactor(gas("gri30",500,0.2*ct.one_atm,"CH4:1.1,O2:2,N2:7.52"),clone=True)
        env = ct.Reservoir(ct.Solution(str(directory / "air.yaml")),clone=True)
        piston = ct.Wall(reacting,argon,A=1,K=0.5e-4,U=100)
        external = ct.Wall(reacting,env,A=1,U=500)
        return [argon,reacting], [piston,external], np.linspace(0,0.12,301), lambda:[piston.expansion_rate]
    if name == "piston":
        left = ct.IdealGasReactor(gas("h2o2",900,ct.one_atm,"H2:2,O2:1,AR:20"),clone=True,volume=0.5)
        right = ct.IdealGasReactor(gas("gri30",900,ct.one_atm,"CO:2,H2O:0.01,O2:5"),clone=True,volume=0.1)
        wall = ct.Wall(left,right,velocity=lambda t:0.0 if t<0.1 else 1e-4*(left.phase.P-right.phase.P))
        return [left,right], [wall], np.linspace(0,0.2,201), lambda:[wall.expansion_rate]
    if name == "custom2":
        phase = ct.Solution(str(directory / "h2o2.yaml"))
        phase.TPY = 920,ct.one_atm,"H2:1,O2:1,N2:3.76"
        env = ct.Reservoir(phase,clone=True)
        reactor = InertialWallReactor(phase,neighbor=env,clone=True)
        wall = ct.Wall(reactor,env)
        return [reactor], [wall], np.linspace(0,0.5,101), lambda:[reactor.v_wall]
    raise ValueError(name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output",type=Path)
    parser.add_argument("--cases",nargs="+",default=["reactor2","piston","custom2"])
    args = parser.parse_args()
    if not ct.__version__.startswith("4."):
        raise RuntimeError("Final moving-wall references require Cantera 4")
    directory = prepare(args.output)
    for name in args.cases:
        reactors,walls,times,velocity = build(name,directory)
        net = ct.ReactorNet(reactors)
        net.rtol,net.atol = 1e-11,1e-20
        net.max_steps = 1000000
        records = {key:[] for key in ("temperature","pressure","volume","mass","internal_energy","velocity")}
        ys = [[] for _ in reactors]
        for t in times:
            if t > 0:
                net.advance(float(t))
            for key,values in (
                ("temperature",[r.T for r in reactors]),("pressure",[r.phase.P for r in reactors]),
                ("volume",[r.volume for r in reactors]),("mass",[r.mass for r in reactors]),
                ("internal_energy",[r.mass*r.phase.int_energy_mass for r in reactors]),("velocity",velocity())):
                records[key].append(values)
            for r,states in zip(reactors,ys):
                states.append(r.phase.Y.copy())
        arrays = {key:np.asarray(values) for key,values in records.items()}
        arrays["time"] = times
        arrays.update({f"Y{i+1}":np.asarray(states) for i,states in enumerate(ys)})
        np.savez(args.output/(name+".reference.npz"),**arrays)
        print(name, "final T", arrays["temperature"][-1],"V",arrays["volume"][-1],flush=True)
    (args.output / "manifest.json").write_text(json.dumps({
        "cantera":ct.__version__,"source_commit":"726522be4e2a13454d8415b7ef799d621f665cf3",
        "cases":args.cases,"rtol":1e-11,"atol":1e-20,
        "source":"https://github.com/Cantera/cantera/tree/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/reactors",
        "notes":"Exact physical source configurations and output coordinates, including initial state; tighter reference tolerances."
    },indent=2)+"\n")


if __name__ == "__main__":
    main()
