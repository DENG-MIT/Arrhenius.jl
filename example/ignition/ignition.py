import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from skopt.space import Space
from skopt.sampler import Lhs

import cantera as ct

np.random.seed(1)

gas = ct.Solution('../../mechanism/gri30.yaml')
gas.TPX = 980, 15*ct.one_atm, 'CH4:1.0,O2:2.0,N2:7.52'

r = ct.IdealGasConstPressureReactor(gas, energy='on')
sim = ct.ReactorNet([r])
time = 0.0
states = ct.SolutionArray(gas, extra=['t'])

print('%10s %10s %10s %14s' % ('t [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))
while sim.time < 1.e-1:
    states.append(r.thermo.state, t=sim.time)
    sim.step()

nodedata = np.vstack((states.t, states.T, states.Y.T)).T
np.savetxt('ignition.dat', nodedata)

plt.clf()
plt.subplot(2, 2, 1)
plt.plot(states.t, states.T, 'o')
plt.xlabel('Time (ms)')
plt.ylabel('Temperature (K)')
plt.subplot(2, 2, 2)
plt.plot(states.t, states.Y[:, gas.species_index('H2')])
plt.xlabel('Time (ms)')
plt.ylabel('CH4 Mole Fraction')
plt.subplot(2, 2, 3)
plt.plot(states.t, states.Y[:, gas.species_index('O2')])
plt.xlabel('Time (ms)')
plt.ylabel('O2 Mole Fraction')
plt.subplot(2, 2, 4)
plt.plot(states.t, states.Y[:, gas.species_index('H2O')])
plt.xlabel('Time (ms)')
plt.ylabel('CO2 Mole Fraction')
plt.tight_layout()
plt.show()
