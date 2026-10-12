"""Exact source settings for Cantera surf_pfr_chain.py and surf_pfr.py.

The chain reports its first reactor outlet at z=0, as the source does. The
direct DAE reference is sampled at exact requested coordinates to avoid the
source CSV loop's pre-step coordinate/post-step state offset.
Both the unmodified chain tolerances and a tighter reference are retained,
because the source defaults accumulate measurable elemental drift over 201 cells.
https://cantera.org/dev/examples/python/reactors/surf_pfr_chain.html
https://cantera.org/dev/examples/python/reactors/surf_pfr.html
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import cantera as ct
import numpy as np


def main(directory):
    directory.mkdir(parents=True,exist_ok=True)
    spec=importlib.util.spec_from_file_location("export_surface",Path(__file__).resolve().parents[1]/"mechanism"/"export_surface.py")
    exporter=importlib.util.module_from_spec(spec);spec.loader.exec_module(exporter)
    exporter.export_surface("methane_pox_on_pt.yaml","Pt_surf",directory/"pox.surface.npz")
    def phase():
        s=ct.Interface("methane_pox_on_pt.yaml","Pt_surf")
        g=s.adjacent["gas"]
        s.TP=1073.15,ct.one_atm
        g.TPX=1073.15,ct.one_atm,"CH4:1,O2:1.5,AR:0.1"
        return s,g
    s,g=phase()
    A,porosity,length=1e-4,0.3,0.003
    mdot=0.4/60*g.density*A*porosity
    n=201;dz=length/(n-1)
    for filename,rtol,atol in (("chain.reference.npz",None,None),("chain.refined.npz",1e-11,1e-22)):
        s,g=phase()
        r=ct.IdealGasReactor(g,energy="off",volume=A*dz*porosity,clone=False)
        upstream=ct.Reservoir(g,clone=True);downstream=ct.Reservoir(g,clone=True)
        rs=ct.ReactorSurface(s,r,A=1e5*r.volume,clone=False)
        feed=ct.MassFlowController(upstream,r,mdot=mdot)
        outlet=ct.PressureController(r,downstream,primary=feed,K=1e-6)
        sim=ct.ReactorNet([r])
        if rtol is not None:sim.rtol,sim.atol=rtol,atol
        X,Y,theta,P=[],[],[],[]
        for k in range(n):
            upstream.phase.TDY=r.phase.TDY
            sim.reinitialize();sim.solve_steady()
            X.append(r.phase.X.copy());Y.append(r.phase.Y.copy());theta.append(rs.coverages.copy());P.append(r.phase.P)
        np.savez(directory/filename,distance=np.arange(n)*dz,mole_fractions=np.array(X).T,
                 mass_fractions=np.array(Y).T,coverages=np.array(theta).T,pressure=P,
                 mass_flow_rate=np.array([mdot]),rtol=np.array([sim.rtol]),atol=np.array([sim.atol]))
    s,g=phase()
    r=ct.FlowReactor(g,clone=False)
    r.area=A;r.mass_flow_rate=mdot;r.energy_enabled=False
    rs=ct.FlowReactorSurface(s,r,clone=False);rs.area=1e5*porosity*A
    sim=ct.ReactorNet([r]);sim.rtol=1e-10;sim.atol=1e-18;sim.max_steps=100000
    sim.initialize()
    distance=np.linspace(0,length,301)
    X,Y,theta,P,rho,speed=[],[],[],[],[],[]
    for z in distance:
        if z>0:sim.advance(z)
        X.append(r.phase.X.copy());Y.append(r.phase.Y.copy());theta.append(rs.coverages.copy())
        P.append(r.phase.P);rho.append(r.phase.density);speed.append(r.speed)
    np.savez(directory/"flow.reference.npz",distance=distance,mole_fractions=np.array(X).T,
             mass_fractions=np.array(Y).T,coverages=np.array(theta).T,pressure=P,density=rho,speed=speed,
             mass_flow_rate=np.array([mdot]))
    manifest={"cantera_version":ct.__version__,"mechanism":"methane_pox_on_pt.yaml",
        "temperature":1073.15,"pressure":ct.one_atm,"composition":"CH4:1,O2:1.5,AR:0.1",
        "length":length,"area":A,"porosity":porosity,"velocity":0.4/60,
        "catalyst_area_per_volume":1e5,"mass_flow_rate":mdot,"chain_reactors":n,
        "flow_sample_points":len(distance),"flow_rtol":sim.rtol,"flow_atol":sim.atol,
        "chain_refined_rtol":1e-11,"chain_refined_atol":1e-22}
    (directory/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("directory",type=Path)
    main(parser.parse_args().directory)
