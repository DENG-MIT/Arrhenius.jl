"""Export a fixed-stoichiometry condensed phase to an Arrhenius.jl SI JSON input.

Usage: python mechanism/export_condensed_phase.py graphite.yaml graphite --output graphite.condensed.json

Cantera is used only while preparing parameters; Julia evaluates all
equilibrium calculations. The phase must be fixed-stoichiometry with exactly
one neutral species, a constant-density (constant-volume) equation of state,
and NASA7 or constant-cp standard thermodynamics. All exported values are SI:
temperatures in K, pressures in Pa, molecular weight in kg/kmol, density in
kg/m^3, molar volume in m^3/kmol, enthalpy in J/kmol, entropy and heat
capacity in J/kmol/K. No evaluated state arrays are exported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import cantera as ct

FORMAT = "arrhenius-condensed-phase-v1"
UNITS = {
    "temperature": "K",
    "pressure": "Pa",
    "molecular_weight": "kg/kmol",
    "density": "kg/m^3",
    "molar_volume": "m^3/kmol",
    "enthalpy": "J/kmol",
    "entropy": "J/kmol/K",
    "heat_capacity": "J/kmol/K",
}
COEFFICIENT_ORDER = (
    "NASA7: [T_mid, a1..a7 high-temperature region, a1..a7 low-temperature region]; "
    "constant-cp: [T0, h0, s0, cp0]"
)


def export_condensed_phase(mechanism, phase, output: Path):
    source = Path(str(mechanism))
    if not source.is_file():
        source = next((Path(directory)/source for directory in ct.get_data_directories()
                       if (Path(directory)/source).is_file()), source)
    if not source.is_file():
        raise FileNotFoundError(f"cannot resolve {mechanism} in Cantera's data directories")
    source = source.resolve()
    condensed = ct.Solution(str(source), phase)
    if condensed.thermo_model != "fixed-stoichiometry":
        raise ValueError("only fixed-stoichiometry condensed phases are supported")
    if condensed.n_species != 1:
        raise ValueError("exactly one condensed species is required")
    species = condensed.species(0)
    if species.charge != 0:
        raise ValueError("charged condensed species are unsupported")
    eos = species.input_data.get("equation-of-state", {})
    eos_model = eos.get("model", "constant-volume") if isinstance(eos, dict) else "constant-volume"
    if eos_model != "constant-volume":
        raise ValueError(
            f"unsupported equation-of-state model {eos_model}: a constant molar volume is required")
    elements = {str(name): float(count) for name, count in species.composition.items()}
    if not elements or not all(math.isfinite(count) for count in elements.values()):
        raise ValueError("element counts must be finite")
    if any(count < 0 for count in elements.values()) or not any(count > 0 for count in elements.values()):
        raise ValueError("element counts must be nonnegative with at least one positive count")
    thermo = species.thermo
    model_name = type(thermo).__name__
    coefficients = [float(value) for value in thermo.coeffs]
    if model_name == "NasaPoly2":
        model = "NASA7"
        if len(coefficients) != 15:
            raise ValueError("NASA7 requires 15 coefficients")
    elif model_name == "ConstantCp":
        model = "constant-cp"
        if len(coefficients) != 4:
            raise ValueError("constant-cp requires 4 coefficients")
        if coefficients[0] <= 0 or coefficients[3] <= 0:
            raise ValueError("constant-cp requires positive T0 and cp0")
    else:
        raise ValueError(f"unsupported species thermo {model_name}: {species.name}")
    if not all(math.isfinite(value) for value in coefficients):
        raise ValueError("thermo coefficients must be finite")
    tmin, tmax = float(thermo.min_temp), float(thermo.max_temp)
    if not (math.isfinite(tmin) and 0 <= tmin < tmax):
        raise ValueError("nonnegative ordered thermodynamic temperature bounds are required")
    if model == "NASA7" and not (math.isfinite(tmax) and 0 < tmin < coefficients[0] < tmax):
        raise ValueError("the NASA7 mid temperature must lie strictly inside the temperature range")
    reference_pressure = float(thermo.reference_pressure)
    molecular_weight = float(condensed.molecular_weights[0])
    density = float(condensed.density)
    if not (math.isfinite(molecular_weight) and molecular_weight > 0):
        raise ValueError("a positive molecular weight is required")
    if not (math.isfinite(density) and density > 0):
        raise ValueError("a positive phase density is required")
    if not (math.isfinite(reference_pressure) and reference_pressure > 0):
        raise ValueError("a positive reference pressure is required")
    payload = {
        "format": FORMAT,
        "generator": "Arrhenius.jl mechanism/export_condensed_phase.py",
        "cantera_version": ct.__version__,
        "source_file": source.name,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "phase_name": str(condensed.name),
        "species_name": str(species.name),
        "elements": elements,
        "charge": 0.0,
        "molecular_weight_kg_per_kmol": molecular_weight,
        "density_kg_per_m3": density,
        "molar_volume_m3_per_kmol": molecular_weight / density,
        "thermo": {
            "model": model,
            "coefficients": coefficients,
            "coefficients_order": COEFFICIENT_ORDER,
            "temperature_range_K": [tmin, tmax if math.isfinite(tmax) else None],
            "reference_pressure_Pa": reference_pressure,
        },
        "units": UNITS,
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    return condensed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mechanism")
    parser.add_argument("phase")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export_condensed_phase(args.mechanism, args.phase, args.output)
    print(f"wrote {args.output}")
