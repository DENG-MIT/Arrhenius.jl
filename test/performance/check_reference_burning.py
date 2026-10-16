import csv, json, sys
from pathlib import Path
import numpy as np

manifest, root, tight_dir = map(Path, sys.argv[1:4])
failures = []
rows = list(csv.DictReader(manifest.open(newline="")))
for row in rows:
    case_id = row["id"]
    spec = json.loads((root / row["spec"]).read_text())
    inlet = max(float(spec.get("T", 300.0)), float(spec.get("T_surface", 0.0) or 0.0))
    with np.load(tight_dir / f"{case_id}.ct.npz", allow_pickle=False) as record:
        peak = float(np.max(record["T"]))
    delta = peak - inlet
    print(f"{case_id},inlet_T={inlet:.9g},peak_T={peak:.9g},rise={delta:.9g}")
    if delta < 500.0:
        failures.append(case_id)
print(f"checked={len(rows)} invalid={len(failures)}")
sys.exit(bool(failures))
