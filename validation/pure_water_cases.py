"""Cantera Reynolds-water oracle for Rankine and the exact vapor-dome grid.

Usage: python validation/pure_water_cases.py OUTPUT_DIRECTORY
Source examples: thermo/rankine.py, rankine_units.py, vapordome.py at Cantera
revision 726522be4e2a13454d8415b7ef799d621f665cf3. Runtime Julia never loads this.
"""
from pathlib import Path
import sys
import numpy as np
import cantera as ct

output = Path(sys.argv[1])
output.mkdir(parents=True, exist_ok=True)
w = ct.Water()
keys = ("T", "P", "Q", "rho", "h", "u", "s", "cp", "cv")


def state(heat_capacity=False):
    return np.array([w.T, w.P, w.Q, w.density, w.h, w.u, w.s,
                     w.cp if heat_capacity else np.nan, w.cv if heat_capacity else np.nan])


degc = np.hstack([np.array([w.min_temp - 273.15, 4, 5, 6, 8]),
                  np.arange(10, 37), np.array([38]), np.arange(40, 100, 5),
                  np.arange(100, 300, 10), np.arange(300, 380, 20),
                  np.arange(370, 374), np.array([w.critical_temperature-273.15])])
dome = []
for temperature in degc+273.15:
    row = []
    for quality in (0., 1.):
        w.TQ = temperature, quality
        row.append(state())
    dome.append(row)
np.savez(output / "water-dome.npz", T=degc+273.15, states=np.array(dome).transpose(2,1,0))
print(f"Vapor dome: {len(dome)} temperatures, critical endpoints {dome[-1]}", flush=True)

cases, results = [], []
for temperature, pressure in ((273.16,1e5),(300.,1e5),(300.,8e5),(373.15,2e5),
                              (500.,1e5),(500.,5e6),(640.,25e6),(650.,1e5),
                              (700.,25e6),(1000.,1e5),(1500.,10e6)):
    w.TP = temperature,pressure
    cases.append([1.,temperature,pressure])
    results.append(state(True))
for temperature in (300.,373.15,500.,640.,647.):
    for quality in (0.,.1,.5,.9,1.):
        w.TQ = temperature,quality
        cases.append([2.,temperature,quality])
        results.append(state(True))
for pressure in (1e4,1e5,8e5,5e6,20e6):
    for quality in (0.,.5,1.):
        w.PQ = pressure,quality
        cases.append([3.,pressure,quality])
        results.append(state())
np.savez(output / "water-states.npz", inputs=np.array(cases).T, states=np.array(results).T)
print(f"State references: {len(cases)}", flush=True)

cycles, metrics = [], []
for initial_temperature,maximum_pressure in ((300.,8e5),((80.33-32)*5/9+273.15,116.03*6894.757293168364)):
    states = []
    w.TQ = initial_temperature,0
    states.append(state())
    h1,s1,p1 = w.h,w.s,w.P
    w.SP = s1,maximum_pressure
    states.append(state())
    pump_work = (w.h-h1)/.6
    w.HP = h1+pump_work,maximum_pressure
    states.append(state())
    h2 = w.h
    w.PQ = maximum_pressure,1
    states.append(state())
    h3,s3 = w.h,w.s
    heat_added = h3-h2
    w.SP = s3,p1
    states.append(state())
    turbine_work = (h3-w.h)*.8
    w.HP = h3-turbine_work,p1
    states.append(state())
    cycles.append(states)
    metrics.append([pump_work,turbine_work,heat_added,(turbine_work-pump_work)/heat_added])
np.savez(output / "water-rankine.npz", states=np.array(cycles).transpose(2,1,0), metrics=np.array(metrics).T)
print(f"Rankine metrics {metrics}; Cantera {ct.__version__}", flush=True)
