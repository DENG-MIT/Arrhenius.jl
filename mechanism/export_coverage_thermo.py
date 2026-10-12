"""Export coverage-dependent surface thermodynamic parameters in SI units.

Usage: python mechanism/export_coverage_thermo.py covdepsurf.yaml covdep_lin
       --output covdep_lin.coverage.json

Cantera resolves species definitions and input units during preprocessing.
The archive contains physical model parameters, with no evaluated state tables.
Julia evaluates the temperature and coverage dependence at runtime.
"""
import argparse
import hashlib
import json
from pathlib import Path

import cantera as ct


def converted(node, key, units, default=0.):
    return node.convert(key, units) if key in node else default


def interpolation(node, property_name, units, model):
    if model == "piecewise-linear":
        keys = [property_name+suffix for suffix in ("-low", "-change", "-high")]
        if not any(key in node for key in keys):
            return {"coverages": [0., 1.], "values": [0., 0.]}
        low = node.convert(keys[0], units)
        change = float(node[keys[1]])
        high = node.convert(keys[2], units)
        if not 0 < change < 1:
            raise ValueError("piecewise-linear coverage breakpoint must be inside (0,1)")
        return {"coverages": [0., change, 1.],
                "values": [0., change*low, change*low+(1-change)*high]}
    values_key = "enthalpies" if property_name == "enthalpy" else "entropies"
    knots_key = property_name+"-coverages"
    knots = node.get(knots_key, [0., 1.])
    values = converted(node, values_key, units, [0., 0.])
    # Cantera's surface implementation supplies zero endpoints before inserting
    # user knots. Preserve those endpoints and its last-value duplicate rule.
    table = {0.: 0., 1.: 0.}
    table.update(zip(knots, values, strict=True))
    ordered = sorted(table)
    if ordered[0] != 0 or ordered[-1] != 1:
        raise ValueError("interpolation knots must lie in [0,1]")
    return {"coverages": ordered, "values": [float(table[k]) for k in ordered]}


def export_coverage_thermo(mechanism, phase, output):
    mechanism, output = Path(mechanism), Path(output)
    surface = ct.Interface(str(mechanism), phase)
    if surface.thermo_model != "coverage-dependent-surface":
        raise ValueError("a coverage-dependent-surface phase is required")
    types, coefficients, dependencies = [], [], []
    names = surface.species_names
    for index, species in enumerate(surface.species()):
        thermo_name = type(species.thermo).__name__
        if thermo_name not in ("NasaPoly2", "ConstantCp"):
            raise ValueError(f"unsupported surface species thermo: {thermo_name}")
        types.append(1 if thermo_name == "NasaPoly2" else 2)
        row = species.thermo.coeffs.tolist()
        coefficients.append(row+[0.]*(15-len(row)))
        for influencing, node in species.input_data.get("coverage-dependencies", {}).items():
            model = node["model"]
            if model not in ("linear", "polynomial", "piecewise-linear", "interpolative"):
                raise ValueError(f"unsupported coverage dependency: {model}")
            dependency = {"target": index+1, "influencing": names.index(influencing)+1,
                          "source_model": model, "heat_capacity_a": converted(node, "heat-capacity-a", "J/kmol/K"),
                          "heat_capacity_b": converted(node, "heat-capacity-b", "J/kmol/K")}
            if ("heat-capacity-a" in node) != ("heat-capacity-b" in node):
                raise ValueError("coverage heat-capacity coefficients a and b must occur together")
            if model in ("linear", "polynomial"):
                dependency["kind"] = "polynomial"
                for prop, units in (("enthalpy", "J/kmol"), ("entropy", "J/kmol/K")):
                    if model == "linear":
                        values = [converted(node, prop, units), 0., 0., 0.]
                    else:
                        values = converted(node, prop+"-coefficients", units, [0.]*4)
                    if len(values) != 4:
                        raise ValueError("coverage polynomials require four coefficients")
                    dependency[prop] = list(values)
            else:
                dependency["kind"] = "interpolative"
                dependency["enthalpy"] = interpolation(node, "enthalpy", "J/kmol", model)
                dependency["entropy"] = interpolation(node, "entropy", "J/kmol/K", model)
            dependencies.append(dependency)
    data = {"format": "arrhenius-coverage-thermo-v1", "source_sha256": hashlib.sha256(mechanism.read_bytes()).hexdigest(),
            "cantera_version": ct.__version__, "phase_name": surface.name, "species_names": names,
            "thermo_type": types, "thermo_coefficients": coefficients,
            "reference_state_coverage": float(surface.input_data.get("reference-state-coverage", 1.)),
            "dependencies": dependencies, "units": "K, J/kmol, J/kmol/K"}
    output.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mechanism", type=Path)
    parser.add_argument("phase")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export_coverage_thermo(args.mechanism, args.phase, args.output)
