"""Compare independently solved counterflow profiles on their native grids."""
import argparse
import json
from pathlib import Path
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("directory",type=Path)
args=parser.parse_args()
summary={}
for kind in ["h2","ethane","methane"]:
    a=np.load(args.directory/f"counterflow-{kind}-native.npz")
    adaptive=np.load(args.directory/f"counterflow-{kind}-reference.npz")
    samegrid=args.directory/f"counterflow-{kind}-reference-samegrid.npz"
    b=np.load(samegrid) if samegrid.is_file() else adaptive
    metrics={"native_points":len(a['grid']),"reference_points":len(b['grid']),
        "same_grid":np.array_equal(a['grid'],b['grid']),
        "native_peak_T_K":float(max(a['T'])),"reference_peak_T_K":float(max(b['T'])),
        "native_Lambda_Pa_per_m2":float(a['Lambda'][0]),"reference_Lambda_Pa_per_m2":float(b['Lambda'][0]),
        "scaled_residual_max":float(abs(a['residual']).max()),
        "minimum_Y":float(a['Y'].min()),"mass_fraction_sum_error":float(abs(a['Y'].sum(axis=0)-1).max())}
    metrics['adaptive_reference_points']=len(adaptive['grid'])
    metrics['adaptive_reference_peak_T_K']=float(max(adaptive['T']))
    metrics['adaptive_reference_Lambda_Pa_per_m2']=float(adaptive['Lambda'][0])
    metrics['adaptive_reference_T_profile_max_abs_K']=float(abs(a['T']-np.interp(a['grid'],adaptive['grid'],adaptive['T'])).max())
    for key in ['T','velocity','spread_rate','Lambda','Y']:
        av,bv=a[key],b[key]
        compare=np.interp(a['grid'],b['grid'],bv) if av.ndim==1 else np.array([
            np.interp(a['grid'],b['grid'],v) for v in bv])
        error=float(abs(av-compare).max())
        metrics[key+"_max_abs_error"]=error
        metrics[key+"_scaled_error"]=error/max(float(abs(compare).max()),1e-30)
    assert metrics['scaled_residual_max']<1e-8
    assert metrics['mass_fraction_sum_error']<1e-10
    assert metrics['T_max_abs_error']<1
    assert metrics['velocity_scaled_error']<.001
    assert metrics['Lambda_scaled_error']<.0001
    assert metrics['Y_max_abs_error']<.001
    summary[kind]=metrics
output=args.directory/'counterflow-comparison.json'
output.write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
