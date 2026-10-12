"""Compare independently initialized native and Cantera catalytic profiles."""
import argparse,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("native",type=Path)
p.add_argument("reference",type=Path)
p.add_argument("output",type=Path)
p.add_argument("--check",action="store_true")
a=p.parse_args()
results={}
for path in sorted(a.native.glob("*.npz")):
    if not (a.reference/path.name).exists(): continue
    n=np.load(path);r=np.load(a.reference/path.name)
    d={"coverage_max_abs":float(np.max(np.abs(n["coverages"]-r["coverages"])))}
    if "grid" in n:
        d["native_points"]=len(n["grid"]);d["reference_points"]=len(r["grid"])
        for key in ["T","Y","velocity","spread_rate","Lambda"]:
            expected=np.stack([np.interp(n["grid"],r["grid"],row) for row in r[key]]) if r[key].ndim==2 else np.interp(n["grid"],r["grid"],r[key])
            d[key+"_max_abs"]=float(np.max(np.abs(n[key]-expected)))
        for key in ["coverage_rates","wall_species_residual","elemental_production","total_mass_production","flow_residual"]:
            d[key+"_max_abs"]=float(np.max(np.abs(n[key])))
        for key in ["diffusive_mass_flux","surface_production_rates"]:
            if key in r:
                d[key+"_comparison_max_abs"]=float(np.max(np.abs(n[key]-r[key])))
    results[path.stem]=d
if a.check:
    assert len(results)==9, "all initial, inert, six H2, and methane outputs are required"
    for name,d in results.items():
        assert d["coverage_max_abs"]<1e-7,(name,d)
        if "T_max_abs" not in d: continue
        assert d["T_max_abs"]<5e-4,(name,d)
        assert d["Y_max_abs"]<1e-7,(name,d)
        assert d["velocity_max_abs"]<1e-7,(name,d)
        assert d["coverage_rates_max_abs"]<1e-6,(name,d)
        assert d["wall_species_residual_max_abs"]<1e-9,(name,d)
        assert d["elemental_production_max_abs"]<1e-12,(name,d)
        assert d["total_mass_production_max_abs"]<1e-12,(name,d)
        assert d["diffusive_mass_flux_comparison_max_abs"]<1e-9,(name,d)
a.output.write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
