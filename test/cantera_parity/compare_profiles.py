"""Compare stored counterflow profiles using the existing validation reader.

Usage: python test/cantera_parity/compare_profiles.py reference.npz native.npz report.json
The comparison uses fixed physical coordinates, with no flame translation.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'validation'))
from flame_benchmarks import profile_errors
ref,native=(dict(np.load(path)) for path in sys.argv[1:3])
names=bytes(ref['species_names_utf8']).decode().split('\n')
ref['species_names']=names
native['species_names']=bytes(native.get('species_names_utf8',ref['species_names_utf8'])).decode().split('\n')
report=profile_errors(ref,native,'burner')
Path(sys.argv[3]).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report['temperature']))
print('all_species_pass:',report['all_species_pass'])
raise SystemExit(0 if report['all_species_pass'] and report['temperature']['normalized_peak_error']<.01 else 1)
