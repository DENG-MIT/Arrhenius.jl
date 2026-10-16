"""Collect native axisymmetric validation, including twin-flame diagnostics."""
import json,sys
from pathlib import Path
import numpy as np

root=Path(sys.argv[1])
n=np.load(root/"counterflow-twin-native-diagnostics.npz")
r=np.load(root/"counterflow-twin-reference-tight-diagnostics.npz")
diagnostics={}
for key in n.files:
    error=float(np.max(np.abs(n[key]-r[key])))
    relative=error/max(float(np.max(np.abs(r[key]))),1e-30)
    diagnostics[key]={"max_abs_error":error,"scaled_error":relative}
    if np.ndim(n[key])==0:
        diagnostics[key].update(native=float(n[key]),reference=float(r[key]))
    assert relative<1e-5,(key,relative)
summary={
    "cantera_source_commit":"726522be4e2a13454d8415b7ef799d621f665cf3",
    "cantera_version":"4.0.0a2",
    "native_initialization":"HP equilibrium and analytic profiles; no Cantera state enters a native solve",
    "reference_comparison":"independent CT initial profiles, common grid, tight CT steady tolerances (1e-9,1e-14)",
    "scope":"neutral ideal-gas mixture-averaged axisymmetric flames; inert wall only",
    "official_examples":[
        "https://cantera.org/dev/examples/python/onedim/premixed_counterflow_flame.html",
        "https://cantera.org/dev/examples/python/onedim/premixed_counterflow_twin_flame.html",
        "https://cantera.org/dev/examples/python/onedim/stagnation_flame.html"],
    "premixed":json.loads((root/"counterflow-premixed-comparison.json").read_text()),
    "twin":json.loads((root/"counterflow-twin-comparison.json").read_text()),
    "wall_sequence":json.loads((root/"counterflow-wall-sequence-comparison.json").read_text()),
    "twin_diagnostics":diagnostics,
}
(root/"counterflow-premixed-summary.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(diagnostics,indent=2))
