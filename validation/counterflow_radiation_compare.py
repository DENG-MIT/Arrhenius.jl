"""Validate complete nonradiating-to-radiating ethane flame profiles."""
import argparse
import json
from pathlib import Path
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("native",type=Path)
parser.add_argument("reference",type=Path)
parser.add_argument("output",type=Path)
args=parser.parse_args()
a,b=np.load(args.native),np.load(args.reference)
metrics={"native_grid_points":len(a['grid']),"reference_grid_points":len(b['grid']),
    "same_grid":bool(np.array_equal(a['grid'],b['grid'])),
    "native_nonradiating_Tmax_K":float(max(a['nonradiating_T'])),
    "native_radiating_Tmax_K":float(max(a['T'])),
    "reference_nonradiating_Tmax_K":float(max(b['nonradiating_T'])),
    "reference_radiating_Tmax_K":float(max(b['T'])),
    "native_peak_heat_loss_W_m3":float(max(a['radiative_heat_loss'])),
    "scaled_residual_max":float(abs(a['residual']).max())}
assert max(a['T'])<max(a['nonradiating_T'])
assert metrics['scaled_residual_max']<1e-8
for key in ['T','Y','velocity','spread_rate','Lambda','radiative_heat_loss']:
    av,bv=a[key],b[key]
    expected=np.interp(a['grid'],b['grid'],bv) if av.ndim==1 else np.array([
        np.interp(a['grid'],b['grid'],v) for v in bv])
    error=float(abs(av-expected).max())
    metrics[key+'_max_abs_error']=error
    metrics[key+'_scaled_error']=error/max(float(abs(expected).max()),1e-30)
assert metrics['T_max_abs_error']<1
assert metrics['Y_max_abs_error']<1e-4
assert metrics['Lambda_scaled_error']<1e-4
assert metrics['radiative_heat_loss_scaled_error']<.001
args.output.write_text(json.dumps(metrics,indent=2)+"\n")
print(json.dumps(metrics,indent=2))
