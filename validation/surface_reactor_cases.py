"""Selected isothermal Pt/H2 closed-vessel and flow-reactor CT4 references.

These exercise gas/surface coupling at the catalytic_combustion.py initial
temperature/composition; geometry and times are selected regression conditions.
They are zero-dimensional reactors, not an impinging-jet flame calculation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import cantera as ct
import numpy as np


def main(directory):
    manifest = {"cantera_version":ct.__version__,"cases":{}}
    for case in ("closed", "flow"):
        surface = ct.Interface(str(directory/"pt_h2.surface.source.yaml"),"Pt_surf")
        gas = surface.adjacent["gas"]
        gas.TPX = 900,ct.one_atm,"H2:0.05,O2:0.21,N2:0.78,AR:0.01"
        surface.TP = gas.TP
        surface.coverages = {"PT(S)":0.5,"O(S)":0.5}
        reactor = ct.IdealGasReactor(gas,energy="off",volume=1e-6,clone=False)
        catalyst = ct.ReactorSurface(surface,reactor,A=0.001,clone=False)
        if case == "flow":
            inlet, exhaust = ct.Reservoir(gas,clone=True),ct.Reservoir(gas,clone=True)
            feed = ct.MassFlowController(inlet,reactor,mdot=reactor.mass/0.001)
            outlet = ct.PressureController(reactor,exhaust,primary=feed,K=1e-7)
        network = ct.ReactorNet([reactor])
        network.rtol,network.atol = 1e-11,1e-20
        network.max_steps = 100000
        times = np.concatenate(([0.0],np.geomspace(1e-12,0.02,141)))
        mass, Y, theta, pressure = [],[],[],[]
        for t in times:
            if t > 0:
                network.advance(t)
            mass.append(reactor.mass)
            Y.append(reactor.phase.Y.copy())
            theta.append(catalyst.coverages.copy())
            pressure.append(reactor.phase.P)
        np.savez(directory/f"surface_reactor_{case}.reference.npz",times=times,
                 mass=np.array(mass),mass_fractions=np.array(Y).T,
                 coverages=np.array(theta).T,pressure=np.array(pressure))
        manifest["cases"][case] = {"temperature":900,"pressure":ct.one_atm,
            "volume":1e-6,"surface_area":0.001,"sample_points":len(times),
            "residence_time":0.001 if case == "flow" else None}
    # Full diamond_cvd.py continuation at the published 20 atomic-H values.
    surface=ct.Interface("diamond.yaml","diamond_100")
    gas=surface.adjacent["gas"]
    X=gas.X.copy()
    kH=gas.species_index("H")
    hydrogen,growth,coverages=[],[],[]
    for _ in range(20):
        gas.TPX=1200,20*ct.one_atm/760,X
        surface.TP=gas.TP
        surface.advance_coverages_to_steady_state()
        bulk=surface.adjacent["diamond"]
        rate=surface.get_net_production_rates(bulk)[0]/bulk.density_mole*1e6*3600
        hydrogen.append(gas.X[kH]);growth.append(rate);coverages.append(surface.coverages.copy())
        X[kH]/=1.4
    np.savez(directory/"diamond_continuation.reference.npz",hydrogen=hydrogen,
             growth=growth,coverages=np.array(coverages).T)
    (directory/"surface_reactor_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    main(parser.parse_args().directory)
