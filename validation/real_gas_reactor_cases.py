"""Independent RK/IG reaction rates and source shock-tube trajectories.

Usage: python real_gas_reactor_cases.py OUTPUT [--trajectories] [--sweep] [--refine]
References use Cantera 4; the generated ideal-gas YAML/sidecar contains the same
reaction/species data as its RK phase. Julia evaluates all runtime calculations.
Use --sweep --refine to retain the exact published inputs and sampling loops
while adding converged states for strict pointwise trajectory comparisons.
Source: https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/reactors/non_ideal_shock_tube.py
"""
from pathlib import Path
import argparse
import importlib.util
import json
import cantera as ct
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("output", type=Path)
parser.add_argument("--trajectories", action="store_true")
parser.add_argument("--sweep", action="store_true")
parser.add_argument("--refine", action="store_true", help="add tight CT4 states at the original sampling times")
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location("sidecar", Path(__file__).parents[1] / "mechanism" / "export_sidecar.py")
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)


def prepare(name, rk):
    rk_path = args.output / f"{name}_RK.yaml"
    ig_path = args.output / f"{name}_IG.yaml"
    rk.write_yaml(rk_path)
    ig = ct.Solution(thermo="ideal-gas", kinetics="bulk", species=rk.species(), reactions=rk.reactions())
    ig.name = name + "_IG"
    ig.write_yaml(ig_path)
    exporter.export(ig_path, Path(str(ig_path) + ".npz"))
    return ig


def rate_reference(name, gas):
    keys = ["T", "P", "density", "X", "net_production_rates", "forward_rates_of_progress",
            "reverse_rates_of_progress", "equilibrium_constants", "concentrations",
            "activities", "partial_molar_int_energies_TV", "cv_mass"]
    records = {k: [] for k in keys}
    for temperature in (760, 1000, 1500, 2500):
        for pressure in (ct.one_atm, 40*ct.one_atm, 100*ct.one_atm):
            # Positive amounts exercise every reversible/third-body/falloff/PLOG rate.
            gas.TPX = temperature, pressure, np.arange(1, gas.n_species+1)
            for key in keys:
                records[key].append(getattr(gas, key))
    values = {k: np.array(v).T for k, v in records.items()}
    values["version_utf8"] = np.frombuffer(ct.__version__.encode(), dtype=np.uint8)
    np.savez(args.output / f"rates_{name}.npz", **values)


dodecane = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_RK")
dodecane_ig = prepare("dodecane", dodecane)
rate_reference("dodecane_RK", dodecane)
rate_reference("dodecane_IG", dodecane_ig)
h2 = ct.Solution("h2o2.yaml", "ohmech-RK")
# A balanced PLOG fixture tests true EOS pressure independently of concentration.
# It is a validation reaction, not a proposed chemical mechanism.
plog = ct.Reaction(equation="H2 + OH <=> H + H2O", rate=ct.PlogRate([
    (1e5, ct.ArrheniusRate(2e7, 0.3, 1e7)),
    (1e6, ct.ArrheniusRate(4e7, 0.2, 1.1e7)),
    (1e7, ct.ArrheniusRate(8e7, 0.1, 1.2e7))]))
plog.duplicate = True
reactions = h2.reactions()
for reaction in reactions:
    if reaction.equation == plog.equation:
        reaction.duplicate = True
h2 = ct.Solution(thermo="Redlich-Kwong", kinetics="bulk", species=h2.species(), reactions=reactions+[plog])
prepare("h2-plog", h2)
rate_reference("h2-plog_RK", h2)
print("Prepared independent rates including physical third bodies, falloff and PLOG.", flush=True)


def trajectory(temperature, phase):
    gas = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_"+phase)
    gas.TP = temperature, 40*ct.one_atm
    gas.set_equivalence_ratio(1, "c12h26", {"o2": 1, "n2": 3.76})
    reactor = ct.IdealGasMoleReactor(gas, clone=False)
    network = ct.ReactorNet([reactor])
    network.preconditioner = ct.AdaptivePreconditioner()
    # Match the published 0.005 s loop and every-20th-step output sampling.
    history = {key: [] for key in ("time", "state", "P", "rho", "rhs", "energy")}

    def append():
        production = gas.net_production_rates
        dY = production*gas.molecular_weights/gas.density
        dT = -np.dot(production, gas.partial_molar_int_energies_TV)/(gas.density*gas.cv_mass)
        history["time"].append(network.time)
        history["state"].append(np.r_[gas.Y, gas.T])
        history["P"].append(gas.P)
        history["rho"].append(gas.density)
        history["rhs"].append(np.r_[dY, dT])
        history["energy"].append(gas.int_energy_mass)

    append()
    counter = 1
    while network.time < 0.005:
        network.step()
        if counter % 20 == 0:
            append()
        counter += 1
    values = {k: np.array(v).T for k, v in history.items()}
    oh = gas.species_index("oh")
    values["oh_index"] = np.array([oh+1])
    values["ignition_delay"] = np.array([values["time"][np.argmax(values["state"][oh])]])
    values["steps"] = np.array([counter-1])
    values["version_utf8"] = np.frombuffer(ct.__version__.encode(), dtype=np.uint8)
    np.savez(args.output / f"trajectory_{phase}_{temperature}.npz", **values)
    print(f"{phase} {temperature} K: {counter-1} steps, tau={values['ignition_delay'][0]:.9g} s", flush=True)


if args.trajectories or args.sweep:
    temperatures = [1000]
    if args.sweep:
        temperatures += [1250, 1170, 1120, 1080, 1040, 1010, 990, 970, 950, 930, 910, 880, 850, 820, 790, 760]
    for temperature in temperatures:
        for phase in ("RK", "IG"):
            trajectory(temperature, phase)

if args.refine:
    for path in sorted(args.output.glob("trajectory_*.npz")):
        values = dict(np.load(path))
        phase = path.stem.split("_")[1]
        gas = ct.Solution("nDodecane_Reitz.yaml", "nDodecane_"+phase)
        gas.TPY = values["state"][-1, 0], values["P"][0], values["state"][:-1, 0]
        reactor = ct.IdealGasMoleReactor(gas, clone=False)
        network = ct.ReactorNet([reactor])
        network.preconditioner = ct.AdaptivePreconditioner()
        network.rtol, network.atol = 1e-11, 1e-22
        refined = {key: [] for key in ("state", "P", "rho", "rhs")}
        for time in values["time"]:
            network.advance(float(time))
            production = gas.net_production_rates
            refined["state"].append(np.r_[gas.Y, gas.T])
            refined["P"].append(gas.P)
            refined["rho"].append(gas.density)
            refined["rhs"].append(np.r_[production*gas.molecular_weights/gas.density,
                -np.dot(production,gas.partial_molar_int_energies_TV)/(gas.density*gas.cv_mass)])
        for key, records in refined.items():
            values["refined_"+key] = np.array(records).T
        values["reference_rtol"] = np.array([network.rtol])
        values["reference_atol"] = np.array([network.atol])
        np.savez(path, **values)
        delta = np.max(np.abs(values["refined_state"][-1]-values["state"][-1]))
        print(f"{path.name}: source-default vs refined max temperature difference {delta:.6g} K", flush=True)
