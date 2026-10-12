"""Exact pinned Cantera 4 thermo/equivalenceRatio.py reference generation.

https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/thermo/equivalenceRatio.py
Reproduces every calculation in source order with the stateful gas object:
mixtures start at the default 300 K / 1 atm and the burnt temperature persists
for all mixtures created after equilibrate('HP') (called with source defaults).
Every recorded state (T, P, X, Y) and every printed scalar is stored in NPZ,
with metadata and mechanism hashes mirrored in JSON.

Run: python validation/equivalence_ratio_case.py mechanism output_directory

equivalence_ratio_calculation(gas) is importable for paired timing harnesses;
it runs the same source-order calculation on any prepared stateful gas object
and returns (mixture_order, states, scalars) without writing files.
"""
from pathlib import Path
import sys
import hashlib
import json
import cantera as ct
import numpy as np

SOURCE_COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"


def equivalence_ratio_calculation(gas):
    """Every thermo/equivalenceRatio.py calculation in source order.

    The gas must start at the source's default 300 K / 1 atm state. Returns
    (mixture_order, states, scalars) with states[name] = (T, P, X, Y).
    """
    mixture_order, states, scalars = [], {}, {}

    def record(name):
        mixture_order.append(name)
        states[name] = (gas.T, gas.P, gas.X.copy(), gas.Y.copy())

    air = "O2:0.21,N2:0.79"

    # Stoichiometric CH4/air mixture; keeps T = 300 K, P = 1 atm.
    gas.set_equivalence_ratio(phi=1.0, fuel="CH4:1", oxidizer=air)
    record("X_stoich_mole")

    # Mass-basis fuel/oxidizer compositions.
    gas.set_equivalence_ratio(1.0, fuel="CH4:1", oxidizer="O2:0.233,N2:0.767", basis="mass")
    record("X_stoich_mass")

    # Unnormalized amounts (233/767) normalize to the same oxidizer as above.
    scalars["phi_mass_basis"] = gas.equivalence_ratio(fuel="CH4:1", oxidizer="O2:233,N2:767",
                                                      basis="mass")
    scalars["phi_assumed_origin"] = gas.equivalence_ratio()
    scalars["Z_default"] = gas.mixture_fraction(fuel="CH4:1", oxidizer=air)
    scalars["Z_bilger"] = gas.mixture_fraction(fuel="CH4:1", oxidizer=air, element="Bilger")
    scalars["Z_carbon"] = gas.mixture_fraction(fuel="CH4:1", oxidizer=air, element="C")
    scalars["Y_CH4_stoich"] = gas["CH4"].Y[0]

    gas.set_mixture_fraction(0.055, fuel="CH4:1", oxidizer=air)
    record("X_Z055")
    scalars["Y_CH4_Z055"] = gas["CH4"].Y[0]

    fuel = {"CH4": 1}  # fuel composition as dictionary instead of string
    gas.set_equivalence_ratio(1, fuel, air)
    record("X_fresh_burnt_case")
    gas.equilibrate("HP")
    record("X_burnt")  # source T persists after this point
    scalars["phi_burnt"] = gas.equivalence_ratio(fuel, air)
    scalars["Z_burnt"] = gas.mixture_fraction(fuel, air)

    fuel_arbitrary = "CH4:1,O2:0.01,CO:0.05,N2:0.1"
    oxidizer_arbitrary = "O2:0.2,N2:0.8,CO2:0.05,CH4:0.01"
    gas.set_equivalence_ratio(2.5, fuel=fuel_arbitrary, oxidizer=oxidizer_arbitrary)
    record("X_arbitrary")
    scalars["phi_arbitrary"] = gas.equivalence_ratio(fuel=fuel_arbitrary,
                                                     oxidizer=oxidizer_arbitrary)
    scalars["phi_arbitrary_assumed"] = gas.equivalence_ratio()

    gas.set_equivalence_ratio(2.0, "H2:1", "O2:1", diluent="H2O", fraction={"diluent": 0.3})
    record("X_diluted_H2O")
    scalars["X_H2O_diluted"] = gas["H2O"].X[0]
    scalars["H2_O2_mole_ratio"] = gas["H2"].X[0] / gas["O2"].X[0]

    gas.set_equivalence_ratio(2.0, "H2", "O2", diluent="CO2:0.5,H2O:0.5",
                              fraction={"fuel": 0.1}, basis="mass")
    record("X_diluted_mass")
    scalars["Y_H2_diluted"] = gas["H2"].Y[0]  # source prints this mass fraction
    scalars["phi_include_species"] = gas.equivalence_ratio(fuel="H2", oxidizer="O2",
                                                           include_species=["H2", "O2"])

    gas.set_equivalence_ratio(2.0, fuel="H2:0.5,H2O:0.5", oxidizer=air)
    record("X_diluted_fuel")
    scalars["phi_diluted_fuel"] = gas.equivalence_ratio(fuel="H2:0.5,H2O:0.5", oxidizer=air)
    return mixture_order, states, scalars


def generate(mechanisms, output):
    """Write equivalence_ratio-gri30.npz/.json into output; returns the manifest."""
    mechanisms, output = Path(mechanisms), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if not ct.__version__.startswith("4."):
        raise RuntimeError("This reference requires Cantera 4")
    mechanism = mechanisms / "gri30.yaml"
    gas = ct.Solution(str(mechanism))
    if gas.n_species != 53:
        raise RuntimeError("stock 53-species gri30.yaml required")
    mixture_order, states, scalars = equivalence_ratio_calculation(gas)

    utf8 = lambda text: np.frombuffer(text.encode(), dtype=np.uint8)
    mechanism_sha256 = hashlib.sha256(mechanism.read_bytes()).hexdigest()
    arrays = {
        "X": np.stack([states[name][2] for name in mixture_order], axis=1),
        "Y": np.stack([states[name][3] for name in mixture_order], axis=1),
        "T_K": np.array([states[name][0] for name in mixture_order]),
        "P_Pa": np.array([states[name][1] for name in mixture_order]),
        "cantera_version_utf8": utf8(ct.__version__),
        "source_commit_utf8": utf8(SOURCE_COMMIT),
        "mechanism_sha256_utf8": utf8(mechanism_sha256),
        "mixture_names_utf8": utf8("\n".join(mixture_order)),
        "scalar_names_utf8": utf8("\n".join(scalars)),
        "species_names_utf8": utf8("\n".join(gas.species_names)),
    }
    arrays.update({f"scalar_{name}": value for name, value in scalars.items()})
    np.savez(output / "equivalence_ratio-gri30.npz", **arrays)
    manifest = dict(
        cantera=ct.__version__,
        source_commit=SOURCE_COMMIT,
        source="samples/python/thermo/equivalenceRatio.py",
        mechanism=mechanism.name,
        mechanism_sha256=mechanism_sha256,
        species=gas.n_species,
        mixtures={name: {"T_K": states[name][0], "P_Pa": states[name][1]}
                  for name in mixture_order},
        scalars={name: float(value) for name, value in scalars.items()},
    )
    (output / "equivalence_ratio-gri30.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv=None):
    mechanisms, output = map(Path, (sys.argv[1:3] if argv is None else argv))
    manifest = generate(mechanisms, output)
    print(json.dumps({"cantera": manifest["cantera"], "mixtures": list(manifest["mixtures"]),
                      "scalars": manifest["scalars"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
