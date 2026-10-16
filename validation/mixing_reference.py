"""Pinned complete mix1 stationary calculation, native inputs and physical checks."""
import argparse
import ast
import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil

import cantera as ct
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def calculation(source, refined=False):
    tree = ast.parse(source.read_text())
    solve = next(i for i, node in enumerate(tree.body) if isinstance(node, ast.Expr)
                 and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Attribute)
                 and node.value.func.attr == "solve_steady")
    # All original construction statements are retained. Only the terminal
    # report printing and diagram construction are outside the calculation.
    setup = compile(ast.Module(body=tree.body[:solve], type_ignores=[]), str(source), "exec")
    last = compile(ast.Module(body=[tree.body[solve]], type_ignores=[]), str(source), "exec")
    namespace = {}
    exec(setup, namespace)
    initial = namespace["mixer"].phase.Y.copy()
    if refined:
        namespace["sim"].rtol, namespace["sim"].atol = 1e-12, 1e-20
    exec(last, namespace)
    phase, reactor = namespace["mixer"].phase, namespace["mixer"]
    flows = np.array([namespace[k].mass_flow_rate for k in ("mfc1", "mfc2", "outlet")])
    incoming = [namespace[k].phase for k in ("res_a", "res_b")]
    mapped = []
    for upstream in incoming:
        mapped.append(np.array([upstream[name].Y[0] if name in upstream.species_names else 0
                                for name in phase.species_names]))
    species_rate = sum(flows[i]*mapped[i] for i in range(2)) - flows[2]*phase.Y + phase.net_production_rates*phase.molecular_weights*reactor.volume
    enthalpy_rate = sum(flows[i]*incoming[i].enthalpy_mass for i in range(2)) - flows[2]*phase.enthalpy_mass
    values = dict(Y=phase.Y, X=phase.X, initial_Y=initial, flows=flows, species_rate=species_rate,
                  source_mapped_Y=np.array(mapped), wdot=phase.net_production_rates,
                  temperature=np.array([phase.T]), pressure=np.array([phase.P]),
                  density=np.array([phase.density]), mass=np.array([reactor.mass]),
                  mean_molecular_weight=np.array([phase.mean_molecular_weight]),
                  enthalpy=np.array([phase.enthalpy_mass]), internal_energy=np.array([phase.int_energy_mass]),
                  entropy=np.array([phase.entropy_mass]), gibbs=np.array([phase.gibbs_mass]),
                  cp=np.array([phase.cp_mass]), cv=np.array([phase.cv_mass]),
                  enthalpy_rate=np.array([enthalpy_rate]))
    record = dict(reactor_species=phase.n_species, air_species=incoming[0].n_species,
                  species_names=phase.species_names, chemistry_enabled=reactor.chemistry_enabled,
                  rtol=namespace["sim"].rtol, atol=namespace["sim"].atol,
                  maximum_species_mass_residual_per_s=float(np.max(np.abs(species_rate))/reactor.mass),
                  enthalpy_balance_W=float(enthalpy_rate), report=phase.report(threshold=0))
    return values, record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root",type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument("--cantera-source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args = p.parse_args()
    source = args.cantera_source / "samples/python/reactors/mix1.py"
    if sha(source) != "87ccbfa629cb6cd494afae82d206c60b215261e16982675264db017be4314e11":
        raise ValueError("mix1 source differs from pinned Cantera commit")
    expected = {name:sha(args.cantera_source / "data" / name) for name in ("gri30.yaml","air.yaml")}
    if expected != {"gri30.yaml":"06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345",
                    "air.yaml":"c352ce80a7d9e673d3436968ed5833799a6dc874f9ea3fee9e01506632d9ceca"}:
        raise ValueError("mixer mechanism files differ from pinned Cantera sources")
    args.output.mkdir(parents=True,exist_ok=False)
    inputs = args.output / "inputs"
    inputs.mkdir()
    for name in expected:
        shutil.copyfile(args.cantera_source / "data" / name, inputs / name)
    exporter = load(args.source_root / "mechanism/export_sidecar.py", "exporter")
    for name in expected:
        exporter.export(inputs/name,inputs/(name+".npz"))
    environment = load(args.source_root / "validation/benchmark_environment.py","environment")
    threads = environment.verify_numerical_threads(set_accelerate=True)
    cases = load(args.source_root / "validation/parallel_transport_cases.py", "cases")
    libraries = cases.loaded_cantera_hashes(environment)
    before = dict(source=sha(source),driver=sha(__file__),inputs={p.name:sha(p) for p in inputs.iterdir()})
    oldcwd = Path.cwd()
    try:
        os.chdir(inputs)
        rows = {}
        for name,refined in (("original",False),("refined",True)):
            values,row = calculation(source,refined)
            np.savez(args.output / (name+".npz"),**values)
            rows[name] = row
    finally:
        os.chdir(oldcwd)
    checks = dict(cantera4=ct.__version__ == "4.0.0a2",
        source_inputs_unchanged=before == dict(source=sha(source),driver=sha(__file__),inputs={p.name:sha(p) for p in inputs.iterdir()}),
        libraries_unchanged=libraries == cases.loaded_cantera_hashes(environment),
        numerical_threads_unchanged=threads == environment.verify_numerical_threads(),
        separate_mechanisms=all(row["air_species"] == 8 and row["reactor_species"] == 53 for row in rows.values()),
        chemistry_enabled=all(row["chemistry_enabled"] for row in rows.values()))
    report = dict(passed=all(checks.values()),checks=checks,source_url="https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/reactors/mix1.py",
        provenance=before,cantera_loaded_sha256=libraries,numerical_threads=threads,host=environment.host_metadata(),cases=rows)
    (args.output / "reference.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    assert report["passed"]


if __name__ == "__main__":
    main()
