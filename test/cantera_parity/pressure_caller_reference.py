"""Print Cantera ROP fixtures for the 10 atm pressure-forwarding tests.

Run from the Arrhenius repository root with Cantera Python installed.
"""
import cantera as ct

mechanism = "mechanism/stagni_2023.yaml"
temperature = 2164.6519745633973
pressure = 10 * ct.one_atm
X = {"NH3": 0.7, "H2": 0.3, "O2": 0.675, "N2": 0.675 * 3.76,
     "HNO": 0.02, "H": 0.01, "NO": 0.02}
gas = ct.Solution(mechanism)
gas.TPX = temperature, pressure, X
initial = gas.Y.copy()
helium = gas.species_index("HE")
nitrogen = gas.species_index("N2")
reaction = 142  # zero-based index for Julia reaction 143
print(f"Cantera kf143={gas.forward_rate_constants[reaction]:.17g}")

for value in (0.0, 1e-25, -1e-25, -1e-7):
    trial = initial.copy()
    trial[helium] = value
    gas.TP = temperature, pressure
    gas.set_unnormalized_mass_fractions(trial)
    gas.TP = temperature, pressure
    print(f"reactor HE={value:.17g} sumY={gas.Y.sum():.17g} "
          f"P={gas.P:.17g} rop143={gas.net_rates_of_progress[reaction]:.17g}")

trial = initial.copy()
trial[helium] = -1e-7
trial[nitrogen] += 1e-7
gas.TP = temperature, pressure
gas.set_unnormalized_mass_fractions(trial)
gas.TP = temperature, pressure
print(f"network HE=-1e-7 sumY={gas.Y.sum():.17g} "
      f"P={gas.P:.17g} rop143={gas.net_rates_of_progress[reaction]:.17g}")
