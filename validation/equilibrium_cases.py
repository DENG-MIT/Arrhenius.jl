"""Generate independent Cantera references for native ideal-gas equilibrium.

Run with the same mechanism files used by Julia, e.g.:
    python validation/equilibrium_cases.py mechanism output_directory

The mixture cases reproduce the calculations in Cantera's thermo/
equivalenceRatio.py. Equilibrium modes follow ThermoPhase.equilibrate:
https://cantera.org/stable/python/thermo.html#cantera.ThermoPhase.equilibrate
https://cantera.org/stable/examples/python/thermo/equivalenceRatio.html
No reference compositions are used to initialize the native equilibrium solver.
"""
from pathlib import Path
import sys
import cantera as ct
import numpy as np

mechanisms, output = map(Path, sys.argv[1:3])
output.mkdir(parents=True, exist_ok=True)
MODES = ("TP", "TV", "HP", "UV", "SP", "SV")


def equilibrium_cases(name):
    gas = ct.Solution(str(mechanisms / f"{name}.yaml"))
    cases = []
    if name == "h2o2":
        for composition in ("H2:2,O2:1", "H2:1.1,O2:1,AR:5"):
            for mode in MODES:
                # Entropy-conserving reactions can cool below the NASA fit range
                # at room temperature, so use a hotter inlet for SP/SV.
                cases.append((mode, 1200.0 if mode in ("SP", "SV") else 300.0,
                              ct.one_atm, composition))
        for composition in ("H2:1", "O2:1", "AR:1", "H2:1,O2:1e-14"):
            cases.append(("TP", 1000.0, ct.one_atm, composition))
        for temperature in (300., 800., 1500., 2500., 3500.):
            for pressure in (.01, 30.):
                cases.append(("TP", temperature, pressure*ct.one_atm,
                              "H2:2,O2:1,AR:3.76"))
    else:
        for phi in (.5, 1., 2.):
            gas.set_equivalence_ratio(phi, "CH4", "O2:1,N2:3.76")
            composition = gas.X.copy()
            for mode in MODES:
                cases.append((mode, 1500.0 if mode in ("SP", "SV") else 300.0,
                              5*ct.one_atm, composition))
        cases.extend(("TP", t, p*ct.one_atm, "CO:1,CO2:1,N2:2")
                     for t, p in ((300., 1.), (1200., .1), (2500., 30.)))
        cases.append(("TP", 1200., ct.one_atm, "CH4:1,O2:2,N2:7.52,AR:1e-16"))

    data = {key: [] for key in ("mode", "T0", "P0", "X0", "T", "P", "X", "Y", "h", "u", "s", "rho")}
    for mode, temperature, pressure, composition in cases:
        gas.TPX = temperature, pressure, composition
        data["mode"].append(MODES.index(mode)+1)
        data["T0"].append(gas.T)
        data["P0"].append(gas.P)
        data["X0"].append(gas.X.copy())
        gas.equilibrate(mode, rtol=1e-11, max_steps=2000)
        for key, value in (("T", gas.T), ("P", gas.P), ("X", gas.X), ("Y", gas.Y),
                           ("h", gas.enthalpy_mass), ("u", gas.int_energy_mass),
                           ("s", gas.entropy_mass), ("rho", gas.density)):
            data[key].append(value.copy() if isinstance(value, np.ndarray) else value)
    arrays = {key: np.stack(value, axis=1) if key in ("X0", "X", "Y") else np.array(value)
              for key, value in data.items()}
    arrays["cantera_version_utf8"] = np.frombuffer(ct.__version__.encode(), dtype=np.uint8)
    np.savez(output / f"equilibrium-{name}.npz", **arrays)
    print(f"{name}: {len(cases)} equilibrium references ({ct.__version__})", flush=True)


for mechanism in ("h2o2", "gri30"):
    equilibrium_cases(mechanism)

gas = ct.Solution(str(mechanisms / "gri30.yaml"))
mixtures = []
for phi, fuel, oxidizer, kwargs in (
    (1., "CH4", "O2:0.21,N2:0.79", {}),
    (1., "CH4", "O2:0.233,N2:0.767", {"basis": "mass"}),
    (2.5, "CH4:1,O2:0.01,CO:0.05,N2:0.1", "O2:0.2,N2:0.8,CO2:0.05,CH4:0.01", {}),
    (2., "H2", "O2", {"diluent": "H2O", "fraction": {"diluent": .3}}),
    (2., "H2", "O2", {"diluent": "CO2:.5,H2O:.5", "fraction": {"fuel": .1}, "basis": "mass"}),
    (2., "H2:.5,H2O:.5", "O2:.21,N2:.79", {}),
):
    gas.set_equivalence_ratio(phi, fuel, oxidizer, **kwargs)
    mixtures.append(gas.X.copy())
gas.set_mixture_fraction(.055, "CH4", "O2:.21,N2:.79")
mixtures.append(gas.X.copy())
np.savez(output / "composition-gri30.npz", X=np.stack(mixtures, axis=1))
print("gri30: 7 composition references", flush=True)
