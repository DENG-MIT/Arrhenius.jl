"""Seven independent inert-wall continuation comparisons and boundary checks."""
import json,sys
from pathlib import Path
import numpy as np

root=Path(sys.argv[1])
results=[]
for m in range(6,13):
    n=np.load(root/f"counterflow-wall-{m}-native.npz")
    r=np.load(root/f"counterflow-wall-{m}-reference-samegrid.npz")
    adaptive=np.load(root/f"counterflow-wall-{m}-reference.npz")
    assert np.array_equal(n["grid"],r["grid"])
    item=dict(mdot=m/100,native_points=len(n["grid"]),reference_adaptive_points=len(adaptive["grid"]),
              native_Tmax_K=float(max(n["T"])),reference_adaptive_Tmax_K=float(max(adaptive["T"])),
              residual_max=float(np.max(np.abs(n["residual"]))),
              inlet_velocity=float(n["velocity"][0]),wall_velocity=float(n["velocity"][-1]))
    for key in ["T","Y","velocity","spread_rate","Lambda"]:
        error=float(np.max(np.abs(n[key]-r[key])))
        item[key+"_max_abs_error"]=error
        item[key+"_scaled_error"]=error/max(float(np.max(np.abs(r[key]))),1e-30)
    assert item["residual_max"]<1e-8
    assert item["T_max_abs_error"]<.01
    assert item["Y_max_abs_error"]<1e-6
    assert item["Lambda_scaled_error"]<1e-6
    assert item["inlet_velocity"]>0 and abs(item["wall_velocity"])<1e-10
    assert abs(n["T"][-1]-500.)<1e-6
    assert abs(n["spread_rate"][-1])<1e-10
    results.append(item)
with open(root/"counterflow-wall-sequence-comparison.json","w") as stream:
    json.dump(results,stream,indent=2)
print(json.dumps(results,indent=2))
