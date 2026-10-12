"""Prepare the complete pinned CO2 example and check native Julia output.

Run with --source-example, --mechanism and --output-dir to prepare inputs.
After equations_of_state_cases.jl, add --julia-result to validate all states.
This script checks calculations; it does not qualify performance.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import shutil
from pathlib import Path

import cantera as ct
import CoolProp
import numpy as np
from CoolProp.CoolProp import PropsSI, PhaseSI

SOURCE_SHA256 = "5a80d530ec9908e10ae45ea3653caa816a8e26047da73c856b9460ae1e968cd0"
ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exporter(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "mechanism" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.export


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-example", type=Path, required=True)
    parser.add_argument("--mechanism", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--julia-result", type=Path)
    args = parser.parse_args()
    if sha(args.source_example) != SOURCE_SHA256:
        parser.error("source example differs from the pinned Cantera calculation")
    if not ct.__version__.startswith("4.0"):
        parser.error("Cantera 4.0 development reference required")
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    mechanism = out / "co2-thermo.yaml"
    if mechanism.resolve() != args.mechanism.resolve():
        shutil.copyfile(args.mechanism, mechanism)
    exporter("export_sidecar")(mechanism, Path(str(mechanism) + ".npz"))
    exporter("export_helmholtz")("CO2", out / "carbon-dioxide.json")
    selected = [n for n in ast.parse(args.source_example.read_text()).body
                if isinstance(n, ast.FunctionDef)
                and n.name in ("get_thermo_Cantera", "get_thermo_CoolProp")]
    assert len(selected) == 2
    namespace = dict(ct=ct, np=np, PropsSI=PropsSI)
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(args.source_example), "exec"), namespace)
    T, p = 300.0, 1e5 * np.linspace(1, 100, 1000)
    arrays = dict(T=np.array([T]), pressure=p)
    for key, phase in (("ideal", "CO2-Ideal"), ("rk", "CO2-RK")):
        gas = ct.Solution(str(mechanism), phase)
        arrays[key] = np.array(namespace["get_thermo_Cantera"](gas, T, p))
        states = ct.SolutionArray(ct.Solution(str(mechanism), phase), len(p))
        states.TPX = T, p, "CO2:1"
        arrays[key + "_density"] = states.density
    arrays["helmholtz"] = np.array(namespace["get_thermo_CoolProp"](T, p))
    arrays["helmholtz_density"] = np.array([PropsSI("D", "P", pi, "T", T, "CO2") for pi in p])
    phases = [PhaseSI("P", pi, "T", T, "CO2") for pi in p]
    # CoolProp labels the compressed liquid above Pc "supercritical_liquid"
    # even though this 300 K sweep remains below Tc.
    assert set(phases) == {"gas", "liquid", "supercritical_liquid"}
    arrays["helmholtz_phase"] = np.array([1 if phase == "gas" else 2 for phase in phases])
    arrays["saturation_pressure"] = np.array([PropsSI("P", "T", T, "Q", 0, "CO2")])
    np.savez(out / "reference.npz", **arrays)
    report = dict(example="thermo/equations_of_state", cantera_version=ct.__version__,
        coolprop_version=CoolProp.__version__, source_sha256=SOURCE_SHA256,
        mechanism_sha256=sha(mechanism), parameter_json_sha256=sha(out / "carbon-dioxide.json"),
        reference_sha256=sha(out / "reference.npz"), performance_qualified=False,
        property_names=["relative_enthalpy", "relative_internal_energy", "relative_entropy", "cp", "cv"],
        property_units=["kJ/kg", "kJ/kg", "kJ/kg/K", "kJ/kg/K", "kJ/kg/K"],
        source_url="https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/thermo/equations_of_state.py")
    if args.julia_result:
        native = np.load(args.julia_result)
        metadata = {k[:-5]: bytes(native[k]).decode() for k in native.files if k.endswith("_utf8")}
        if metadata["parameter_sha256"] != report["parameter_json_sha256"]:
            raise ValueError("native Helmholtz parameters differ from the reference export")
        if metadata["example_sha256"] != sha(ROOT / "example/thermodynamics/equations_of_state.jl"):
            raise ValueError("native example source changed after execution")
        checks = {}
        for key, expected in arrays.items():
            actual = native[key]
            rtol, atol = (1e-7, 1e-7) if key in ("ideal", "rk", "helmholtz") else (1e-8, 1e-8)
            if key == "helmholtz_phase":
                passed = np.array_equal(actual, expected)
            else:
                passed = actual.shape == expected.shape and np.allclose(actual, expected, rtol=rtol, atol=atol)
            checks[key] = dict(passed=bool(passed), rtol=rtol, atol=atol,
                maximum_absolute_error=float(np.max(np.abs(actual - expected))))
            if key in ("ideal", "rk", "helmholtz"):
                checks[key]["maximum_absolute_error_by_property"] = np.max(np.abs(actual - expected), axis=1).tolist()
        report.update(checks=checks, correctness_pass=all(c["passed"] for c in checks.values()),
                      native_artifact_sha256=sha(args.julia_result), native_metadata=metadata)
    (out / "checks.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf8")
    print(json.dumps(report, indent=2))
    if report.get("correctness_pass") is False:
        raise SystemExit("native CO2 comparison failed")


if __name__ == "__main__":
    main()
