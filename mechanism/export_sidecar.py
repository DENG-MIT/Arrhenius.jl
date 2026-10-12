"""Create an Arrhenius.jl mechanism sidecar with Cantera 3.2+ and NumPy."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import cantera as ct
import numpy as np
from ruamel.yaml import YAML


SUPPORTED_RATE_TYPES = {
    "ArrheniusRate",
    "LindemannRate",
    "PlogRate",
    "TroeRate",
    "BlowersMaselRate",
}
J_PER_KMOL_PER_CAL_PER_MOL = 4184.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_import(mechanism: Path, filename: str) -> Path:
    candidate = Path(filename)
    if candidate.is_absolute():
        candidates = [candidate]
    else:
        candidates = [mechanism.parent / filename] + [
            Path(directory) / filename for directory in ct.get_data_directories()
        ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError(f"cannot resolve imported mechanism file: {filename}")


def phase_selection(mechanism: Path, gas: ct.Solution) -> bytes:
    yaml = YAML(typ="safe")
    with mechanism.open("r", encoding="utf-8") as stream:
        document = yaml.load(stream)
    imports: set[str] = set()
    phases = (document or {}).get("phases") or []
    if isinstance(phases, list) and phases:
        first = phases[0]
        if isinstance(first, dict):
            for field in ("species", "reactions"):
                selector = first.get(field)
                entries = selector if isinstance(selector, list) else [selector]
                for entry in entries:
                    if isinstance(entry, dict):
                        for key in entry:
                            key = str(key)
                            if "/" in key:
                                imports.add(key.rsplit("/", 1)[0])
    dependencies = {
        name: sha256(resolve_import(mechanism, name)) for name in sorted(imports)
    }
    selection = dict(
        species_names=list(gas.species_names),
        n_reactions=int(gas.n_reactions),
        dependencies=dependencies,
    )
    return json.dumps(selection, ensure_ascii=True).encode("utf-8")


def arrhenius_row(data: dict[str, Any]) -> tuple[float, float, float]:
    return (
        float(data["A"]),
        float(data.get("b", 0.0)),
        float(data.get("Ea", 0.0)) / J_PER_KMOL_PER_CAL_PER_MOL,
    )


def export(mechanism: Path, output: Path) -> None:
    gas = ct.Solution(str(mechanism))
    reactions = gas.reactions()
    species_index = {name: index for index, name in enumerate(gas.species_names)}
    unsupported = sorted(
        {type(reaction.rate).__name__ for reaction in reactions}
        - SUPPORTED_RATE_TYPES
    )
    if unsupported:
        raise ValueError(f"unsupported Cantera rate types: {unsupported}")

    reactant_stoich = np.asarray(gas.reactant_stoich_coeffs, dtype=np.float64)
    product_stoich = np.asarray(gas.product_stoich_coeffs, dtype=np.float64)
    reactant_orders = reactant_stoich.copy()
    reversible = np.asarray(
        [reaction.reversible for reaction in reactions], dtype=np.bool_
    )
    arrhenius = np.zeros((gas.n_reactions, 3), dtype=np.float64)
    efficiencies = np.zeros((gas.n_species, gas.n_reactions), dtype=np.float64)
    falloff_low: list[tuple[float, float, float]] = []
    troe_rows: list[tuple[float, float, float, float]] = []
    three_body_indices: list[int] = []
    falloff_indices: list[int] = []
    falloff_troe_indices: list[int] = []
    blowers_indices = []
    blowers_coefficients = []

    plog_reaction_indices: list[int] = []
    plog_group_offsets = [1]
    plog_pressures: list[float] = []
    plog_rate_offsets = [1]
    plog_arrhenius: list[tuple[float, float, float]] = []

    for reaction_index, reaction in enumerate(reactions):
        for name, order in reaction.orders.items():
            reactant_orders[species_index[name], reaction_index] = float(order)

        if reaction.third_body is not None:
            efficiencies[:, reaction_index] = reaction.third_body.default_efficiency
            for name, efficiency in reaction.third_body.efficiencies.items():
                efficiencies[species_index[name], reaction_index] = float(efficiency)

        rate = reaction.rate
        rate_type = type(rate).__name__
        rate_data = rate.input_data
        if rate_type == "BlowersMaselRate":
            coefficients = rate_data["rate-constant"]
            blowers_indices.append(reaction_index + 1)
            blowers_coefficients.append([float(coefficients[key]) for key in ("A","b","Ea0","w")])
            arrhenius[reaction_index,:] = (coefficients["A"],coefficients["b"],0.0)
            if reaction.third_body is not None:
                three_body_indices.append(reaction_index + 1)
        elif rate_type == "ArrheniusRate":
            arrhenius[reaction_index, :] = arrhenius_row(rate_data["rate-constant"])
            if reaction.third_body is not None:
                three_body_indices.append(reaction_index + 1)
        elif rate_type in {"LindemannRate", "TroeRate"}:
            falloff_indices.append(reaction_index + 1)
            arrhenius[reaction_index, :] = arrhenius_row(
                rate_data["high-P-rate-constant"]
            )
            falloff_low.append(arrhenius_row(rate_data["low-P-rate-constant"]))
            if rate_type == "TroeRate":
                troe = rate_data["Troe"]
                troe_rows.append(
                    (
                        float(troe["A"]),
                        float(troe["T1"]),
                        float(troe.get("T2", 1.0e30)),
                        float(troe["T3"]),
                    )
                )
                falloff_troe_indices.append(len(troe_rows))
            else:
                falloff_troe_indices.append(-1)
        elif rate_type == "PlogRate":
            grouped: dict[float, list[tuple[float, float, float]]] = defaultdict(list)
            for item in rate_data["rate-constants"]:
                grouped[float(item["P"])].append(arrhenius_row(item))
            groups = sorted(grouped.items())
            if not groups:
                raise ValueError(f"PLOG reaction {reaction_index + 1} has no rates")
            arrhenius[reaction_index, :] = groups[0][1][0]
            plog_reaction_indices.append(reaction_index + 1)
            for pressure, rows in groups:
                plog_pressures.append(pressure)
                plog_arrhenius.extend(rows)
                plog_rate_offsets.append(len(plog_arrhenius) + 1)
            plog_group_offsets.append(len(plog_pressures) + 1)

    transport = {}
    if gas.transport_model in {"mixture-averaged", "multicomponent", "unity-Lewis-number", "ionized-gas"}:
        # Cantera's native degree-four fits in log(T), in ascending order.
        # Julia evaluates the fits and mixture rules; no Python runtime is used.
        # Keep the active ionized model: its binary fits include ion-neutral
        # collision corrections and the O2/O2- resonant-collision override.
        transport = {
            "transport_model_utf8": np.frombuffer(gas.transport_model.encode(), dtype=np.uint8),
            "species_viscosities_poly": np.array([
                gas.get_viscosity_polynomial(k) for k in range(gas.n_species)
            ]).T,
            "thermal_conductivity_poly": np.array([
                gas.get_thermal_conductivity_polynomial(k) for k in range(gas.n_species)
            ]).T,
            "binary_diff_coeffs_poly": np.array([
                gas.get_binary_diff_coeffs_polynomial(i, j)
                for j in range(gas.n_species) for i in range(gas.n_species)
            ]).T,
        }
    output.parent.mkdir(parents=True, exist_ok=True)
    selection = phase_selection(mechanism, gas)
    payload = dict(
        molecular_weights=np.asarray(gas.molecular_weights, dtype=np.float64),
        reactant_stoich_coeffs=reactant_stoich,
        product_stoich_coeffs=product_stoich,
        reactant_orders=reactant_orders,
        is_reversible=reversible,
        Arrhenius_coeffs=arrhenius,
        efficiencies_coeffs=efficiencies,
        Arrhenius_A0=np.asarray([row[0] for row in falloff_low]),
        Arrhenius_b0=np.asarray([row[1] for row in falloff_low]),
        Arrhenius_Ea0=np.asarray([row[2] for row in falloff_low]),
        Troe_A=np.asarray([row[0] for row in troe_rows]),
        Troe_T1=np.asarray([row[1] for row in troe_rows]),
        Troe_T2=np.asarray([row[2] for row in troe_rows]),
        Troe_T3=np.asarray([row[3] for row in troe_rows]),
        index_three_body=np.asarray(three_body_indices, dtype=np.int64),
        index_falloff=np.asarray(falloff_indices, dtype=np.int64),
        index_falloff_Troe=np.asarray(falloff_troe_indices, dtype=np.int64),
        BlowersMasel_reaction_indices=np.asarray(blowers_indices,dtype=np.int64),
        BlowersMasel_coefficients=np.asarray(blowers_coefficients,dtype=float).reshape((-1,4)),
        Plog_reaction_indices=np.asarray(plog_reaction_indices, dtype=np.int64),
        Plog_group_offsets=np.asarray(plog_group_offsets, dtype=np.int64),
        Plog_pressures=np.asarray(plog_pressures, dtype=np.float64),
        Plog_rate_offsets=np.asarray(plog_rate_offsets, dtype=np.int64),
        Plog_Arrhenius=np.asarray(plog_arrhenius, dtype=np.float64).reshape((-1, 3)),
        sidecar_format_utf8=np.frombuffer(b"arrhenius-sidecar-v4", dtype=np.uint8),
        source_sha256_utf8=np.frombuffer(sha256(mechanism).encode(), dtype=np.uint8),
        phase_selection_utf8=np.frombuffer(selection, dtype=np.uint8),
        cantera_version_utf8=np.frombuffer(ct.__version__.encode(), dtype=np.uint8),
        **transport,
    )
    # NPZ.jl 0.4 / ZipFile cannot read a zero-length member at EOF. Optional
    # reaction families are represented by absent keys, as accepted by the loader.
    np.savez(output, **{key: value for key, value in payload.items() if value.size})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mechanism", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mechanism = args.mechanism.resolve()
    output = (args.output or Path(f"{mechanism}.npz")).resolve()
    export(mechanism, output)
    print(output)


if __name__ == "__main__":
    main()
