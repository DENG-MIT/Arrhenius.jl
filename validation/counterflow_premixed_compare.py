"""Compare independent solutions, preserving mesh and profile error evidence."""
import sys,json
import numpy as np
n,r=map(np.load,sys.argv[1:3])
result={"native_points":len(n["grid"]),"reference_points":len(r["grid"]),
        "same_grid":bool(np.array_equal(n["grid"],r["grid"])),
        "native_Tmax_K":float(max(n["T"])),"reference_Tmax_K":float(max(r["T"])),
        "residual_max":float(np.max(np.abs(n["residual"])))}
for key in ["T","Y","velocity","spread_rate","Lambda"]:
    ref=r[key]
    sampled=np.interp(n["grid"],r["grid"],ref) if ref.ndim==1 else np.array([np.interp(n["grid"],r["grid"],v) for v in ref])
    error=float(np.max(np.abs(n[key]-sampled)))
    result[key+"_max_abs_error"]=error
    result[key+"_scaled_error"]=error/max(float(np.max(np.abs(ref))),1e-30)
print(json.dumps(result,indent=2))
if len(sys.argv)>3:
    with open(sys.argv[3],"w") as out: json.dump(result,out,indent=2)
assert result["residual_max"]<1e-8
assert result["T_max_abs_error"]<1
assert result["Y_max_abs_error"]<1e-4
assert result["Lambda_scaled_error"]<1e-4
