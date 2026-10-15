"""Restricted gate stage of h/validate_heldout.sh using its original helpers.

Inputs are a file of authorized IDs and paired outputs from
h/nh3_native_server.jl and h/nh3_cantera.py. This deliberately supports a
branch and fewer than twelve cases without reading the remaining held-out IDs.
Artifacts are named ID.ct.npz and ID.native.npz. All original accuracy gates
are preserved; this is an accuracy check, not a timing comparison.
"""
from pathlib import Path
import argparse,json,os,subprocess
p=argparse.ArgumentParser()
p.add_argument('root',type=Path);p.add_argument('ids',type=Path)
p.add_argument('artifacts',type=Path);p.add_argument('out',type=Path)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
ids=[v.strip() for v in a.ids.read_text().splitlines() if v.strip()]
assert ids and len(ids)==len(set(ids)), 'empty or duplicate restricted selection'
accuracy_gates=('speed_SL','temperature','species_profiles','reference_converged','elements')
rows=[]
for case_id in ids:
    specs=[a.root/d/(case_id+'.json') for d in ('specs_free_ma','specs_burner_respec','specs_fixed_small','ct_survey/specs_cf2') if (a.root/d/(case_id+'.json')).is_file()]
    refs=[a.root/d/(case_id+'.npz') for d in ('refs16','refs_cf2/counterflow') if (a.root/d/(case_id+'.npz')).is_file()]
    assert len(specs)==len(refs)==1, 'ambiguous or missing case mapping'
    result=a.out/(case_id+'.json')
    subprocess.run([os.environ['TRACE_CAP'],'60',os.environ['PY'],str(a.root/'h/nh3_gate.py'),str(specs[0]),str(refs[0]),str(a.artifacts/(case_id+'.ct.npz')),str(a.artifacts/(case_id+'.native.npz')),str(result)],cwd=a.root,check=True)
    data=json.loads(result.read_text())
    failures=[k for k in accuracy_gates if not data['gates'].get(k,False)]
    rows.append({'id':case_id,'pass':not failures,'failing_gates':failures})
summary={'cases':len(rows),'accuracy_pass':sum(r['pass'] for r in rows),'results':rows}
(a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
raise SystemExit(0 if all(r['pass'] for r in rows) else 1)
