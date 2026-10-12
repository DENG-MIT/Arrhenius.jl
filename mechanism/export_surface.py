"""Export numeric ideal-surface kinetics and a standalone adjacent ideal gas.

Usage: python mechanism/export_surface.py ptcombust.yaml Pt_surf --output pt.surface.npz
Cantera is used only while preparing parameters; Julia evaluates all runtime rates.
Supported thermo: NASA7 and constant-cp; adjacent phases: one ideal gas and
persistent fixed-stoichiometry solids. Electrochemistry and nonideal surfaces
are rejected. Units in the archive are kmol, m, s, K, and J.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
from pathlib import Path

import cantera as ct
import numpy as np


def utf8(text):
    return np.frombuffer(text.encode("utf-8"), dtype=np.uint8)


def export_surface(mechanism, phase, output: Path, *, export_gas=True):
    surface = ct.Interface(str(mechanism), phase)
    if surface.thermo_model != "ideal-surface":
        raise ValueError("only ideal-surface thermodynamics is supported")
    gases = [p for p in surface.adjacent.values() if p.thermo_model == "ideal-gas"]
    if len(gases) != 1:
        raise ValueError("exactly one adjacent ideal gas is required")
    gas = gases[0]
    solids = [p for p in surface.adjacent.values() if p is not gas]
    if any(p.thermo_model != "fixed-stoichiometry" or p.n_species != 1 for p in solids):
        raise ValueError("additional phases must be single-species fixed-stoichiometry solids")
    phases = [surface, gas, *solids]
    names = [name for p in phases for name in p.species_names]
    if len(set(names)) != len(names):
        raise ValueError("species names must be unique across phases")
    indices = {name: k for k, name in enumerate(names)}
    order = [surface.kinetics_species_index(name) for name in names]
    species = [sp for p in phases for sp in p.species()]
    if any(sp.charge != 0 for sp in species):
        raise ValueError("charged species and electrochemical reactions are unsupported")
    ns, ng, nr, nt = surface.n_species, gas.n_species, surface.n_reactions, len(names)
    reactants = np.asarray(surface.reactant_stoich_coeffs)[order, :]
    products = np.asarray(surface.product_stoich_coeffs)[order, :]
    orders = reactants.copy()
    arrhenius = np.zeros((nr, 3))
    cov_a, cov_m, cov_e = np.zeros((nr, ns)), np.zeros((nr, ns)), np.zeros((nr, ns, 4))
    sticking, motz_wise, sticking_order, sticking_factor = (
        np.zeros(nr, dtype=np.int64), np.zeros(nr, dtype=np.bool_), np.zeros(nr), np.zeros(nr))
    reversible = np.zeros(nr, dtype=np.bool_)
    sizes = np.array([sp.size for sp in surface.species()])
    for j, reaction in enumerate(surface.reactions()):
        rate = reaction.rate
        rate_type = type(rate).__name__
        if rate_type not in ("InterfaceArrheniusRate", "StickingArrheniusRate"):
            raise ValueError(f"unsupported surface rate {rate_type}: {reaction.equation}")
        data = rate.input_data
        if any(k in data for k in ("beta", "exchange-current-density-formulation")):
            raise ValueError("electrochemical rate corrections are unsupported")
        base = data["sticking-coefficient" if rate_type == "StickingArrheniusRate" else "rate-constant"]
        arrhenius[j] = [base["A"], base.get("b", 0.0), base.get("Ea", 0.0)]
        if arrhenius[j, 0] < 0:
            raise ValueError("negative pre-exponential factors are unsupported")
        for name, value in reaction.orders.items():
            if value < 0:
                raise ValueError("negative reaction orders are unsupported")
            orders[indices[name], j] = value
        reversible[j] = reaction.reversible
        for name, dependency in rate.coverage_dependencies.items():
            k = indices[name]
            if k >= ns:
                raise ValueError("coverage dependency must name a surface species")
            cov_a[j, k], cov_m[j, k] = dependency["a"], dependency["m"]
            energy = np.atleast_1d(dependency["E"])
            if len(energy) > 4:
                raise ValueError("coverage energy supports at most fourth order")
            cov_e[j, k, :len(energy)] = energy
        if rate_type == "StickingArrheniusRate":
            k = indices[rate.sticking_species]
            if not ns <= k < ns + ng or orders[k, j] != 1:
                raise ValueError("sticking species must be a first-order ideal-gas reactant")
            if any(orders[q, j] != 0 for q in range(ns, nt) if q != k):
                raise ValueError("sticking reactions with additional non-surface reactants are unsupported")
            sticking[j], motz_wise[j] = k + 1, rate.motz_wise_correction
            sticking_order[j] = orders[:ns, j].sum()
            sticking_factor[j] = np.sqrt(ct.gas_constant/(2*np.pi*gas.molecular_weights[k-ns]))
            sticking_factor[j] *= np.prod(sizes ** orders[:ns, j])
    thermo_type = np.zeros(nt, dtype=np.int64)
    thermo_coefficients = np.zeros((nt, 15))
    for k, sp in enumerate(species):
        model = type(sp.thermo).__name__
        if model == "NasaPoly2":
            thermo_type[k] = 1
        elif model == "ConstantCp":
            thermo_type[k] = 2
        else:
            raise ValueError(f"unsupported species thermo {model}: {sp.name}")
        thermo_coefficients[k, :len(sp.thermo.coeffs)] = sp.thermo.coeffs
    elements = sorted({element for sp in species for element in sp.composition})
    elemental = np.array([[sp.composition.get(e, 0.0) for sp in species] for e in elements])
    nu = products - reactants
    if not np.allclose(sizes @ nu[:ns], 0, atol=1e-12, rtol=0):
        raise ValueError("surface reactions do not conserve site occupancy")
    if not np.allclose(elemental @ nu, 0, atol=1e-12, rtol=0):
        raise ValueError("surface reactions do not conserve elements")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    gas_file = output.with_suffix(".gas.yaml")
    if export_gas:
        gas.write_yaml(str(gas_file))
        if gas.n_reactions == 0:
            with gas_file.open("a") as stream:
                stream.write("\nreactions: []\n")
        spec = importlib.util.spec_from_file_location("export_sidecar", Path(__file__).with_name("export_sidecar.py"))
        exporter = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(exporter)
        exporter.export(gas_file, Path(str(gas_file) + ".npz"))
    source_file = output.with_suffix(".source.yaml")
    surface.write_yaml(str(source_file))
    payload = dict(
        surface_format_utf8=utf8("arrhenius-surface-v1"),
        source_name_utf8=utf8(str(mechanism)), phase_name_utf8=utf8(phase),
        cantera_version_utf8=utf8(ct.__version__),
        source_parameters_sha256_utf8=utf8(hashlib.sha256(source_file.read_bytes()).hexdigest()),
        gas_file_utf8=utf8(gas_file.name if export_gas else ""),
        species_names_utf8=utf8("\n".join(names)), element_names_utf8=utf8("\n".join(elements)),
        n_surface=np.array([ns]), n_gas=np.array([ng]),
        site_density=np.array([surface.site_density]), site_sizes=sizes,
        molecular_weights=np.concatenate([p.molecular_weights for p in phases]),
        elemental_matrix=elemental, reactants=reactants, products=products, orders=orders,
        arrhenius=arrhenius, reversible=reversible,
        coverage_a=cov_a, coverage_m=cov_m, coverage_energy=cov_e,
        sticking_species=sticking, sticking_order=sticking_order,
        sticking_factor=sticking_factor, motz_wise=motz_wise,
        thermo_type=thermo_type, thermo_coefficients=thermo_coefficients,
        reference_pressure=np.array([sp.thermo.reference_pressure for sp in species]),
        bulk_molar_volumes=np.array([1/p.density_mole for p in solids]),
        initial_coverages=surface.coverages, initial_mole_fractions=gas.X,
        initial_temperature=np.array([surface.T]), initial_pressure=np.array([surface.P]),
        reaction_equations_utf8=utf8("\n".join(surface.reaction_equations())))
    # NPZ.jl cannot load zero-length NumPy array payloads on all Julia versions.
    np.savez(output, **{key: value for key, value in payload.items() if value.size})
    return surface, gas_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mechanism")
    parser.add_argument("phase")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export_surface(args.mechanism, args.phase, args.output)
