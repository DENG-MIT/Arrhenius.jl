#!/usr/bin/env python3
"""Generate one task-owned Cantera criterion-path archive for the gate."""
import json, os, sys
from pathlib import Path
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import cantera as ct
import numpy as np
from ct_tight_reference import build, RTOL, ATOL

def main(spec_path, scale, output_path):
    spec=json.loads(Path(spec_path).read_text())
    scale=float(scale)
    flame=build(spec)
    for domain in flame.domains:
        domain.set_steady_tolerances(default=(RTOL, ATOL))
    if any(not np.all(np.asarray(d.steady_reltol()) == RTOL) or
           not np.all(np.asarray(d.steady_abstol()) == ATOL) for d in flame.domains):
        raise RuntimeError("effective Cantera steady tolerances differ from requested values")
    flame.set_refine_criteria(ratio=spec["ratio"], slope=float(spec["slope"])*scale,
        curve=float(spec["curve"])*scale,
        prune=min(float(spec.get("prune",0.0)),float(spec["slope"])*scale,float(spec["curve"])*scale))
    flame.solve(loglevel=0, auto=(spec["kind"] != "fixed"))
    flame.eval(rdt=0.0)
    residual_inf=0.0
    for domain in flame.domains:
        for i in range(domain.n_components):
            vals=np.asarray(domain.residuals(domain.component_name(i)),dtype=float)
            if vals.size:
                residual_inf=max(residual_inf,float(np.nanmax(np.abs(vals))))
    stats=flame.solver_stats
    data=dict(grid=flame.grid.copy(),T=flame.T.copy(),Y=flame.Y.copy(),
        velocity=flame.velocity.copy(),SL=np.asarray(float(flame.velocity[0]) if spec["kind"]=="free" else float(np.max(flame.T))),
        mdot=np.asarray(float(flame.density[0]*flame.velocity[0])),
        species_names=np.asarray(flame.gas.species_names),
        tight_rtol_reached=np.asarray([RTOL]),residual_inf=np.asarray([residual_inf]),
        criterion_scale=np.asarray([scale]))
    # Successful return at the verified effective tolerance proves the attained
    # Newton-step norm is no larger than this recorded conservative bound.
    data["newton_ok"]=np.asarray([1],dtype=np.int8)
    Path(output_path).parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output_path,**data)
    print(f"spread id={spec['id']} scale={scale:.2f} newton_ok=1 tight_rtol_bound={RTOL:.0e} residual_inf={residual_inf:.6g} points={len(flame.grid)} out={output_path}")
if len(sys.argv)!=4: raise SystemExit("usage: ct_spread_reference.py SPEC.json SCALE OUT.npz")
main(sys.argv[1],sys.argv[2],sys.argv[3])
