"""Gate all paired branchcheck outputs with the existing harness and compare.

Usage: summarize_branchcheck.py HARNESS_ROOT BRANCHCHECK_OUT
TRACE_CAP and PY name the capped runner and Cantera Python. Inputs are the
manifest, main/branch counts.tsv and native NPZs, and common ct/*.npz outputs.
CounterflowStats counts include the adapter warm and measured solve. Zeros in
other flame families mean those counters are not instrumented there.
"""
import csv,json,os,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1]).resolve(); out=Path(sys.argv[2]).resolve()
rows=list(csv.DictReader(open(out/'manifest.csv')))
summary=[]
for r in rows:
    entry={'id':r['id']}
    for label in ('main','branch'):
        native=out/label/(r['id']+'.npz'); target=out/label/(r['id']+'.json')
        if not native.is_file():
            error=native.with_suffix('.npz.error')
            entry[label]={'pass':False,'error':error.read_text() if error.exists() else 'missing native artifact'}
            continue
        subprocess.run([os.environ['TRACE_CAP'],'60',os.environ['PY'],str(root/'h/nh3_gate.py'),r['spec'],r['reference'],str(out/'ct'/(r['id']+'.npz')),str(native),str(target)],cwd=root,check=True,stdout=subprocess.DEVNULL)
        d=json.loads(target.read_text())
        entry[label]={'pass':d['accuracy_pass'],'failed_gates':[k for k in ('speed_SL','temperature','species_profiles','reference_converged','elements') if not d['gates'][k]]}
    summary.append(entry)
counts={label:{r['id']:r for r in csv.DictReader(open(out/label/'counts.tsv'),delimiter='\t')} for label in ('main','branch')}
regressions=[]
for r in summary:
    a=counts['main'][r['id']];b=counts['branch'][r['id']]
    r['counts']={label:counts[label][r['id']] for label in ('main','branch')}
    if r['main']['pass'] and not r['branch']['pass']: regressions.append(r['id'])
result={'cases':len(summary),'main_pass':sum(r['main']['pass'] for r in summary),'branch_pass':sum(r['branch']['pass'] for r in summary),'accuracy_regressions':regressions,'results':summary}
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['| Case | Main | Branch | Jacobians main/branch | Newton main/branch | Steps main/branch |','|---|---|---|---|---|---|']
for r in summary:
    a=r['counts']['main'];b=r['counts']['branch']
    lines.append(f"| {r['id']} | {r['main']['pass']} | {r['branch']['pass']} | {a['jacobians']}/{b['jacobians']} | {a['newton']}/{b['newton']} | {a['steps']}/{b['steps']} |")
lines+=['',f"Accuracy: main {result['main_pass']}/{len(summary)}, branch {result['branch_pass']}/{len(summary)}; regressions={len(regressions)}.", 'Counts include the adapter warm and measured solves.']
(out/'summary.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='results'}))
for r in summary:
    a=r['counts']['main'];b=r['counts']['branch']
    if a['jacobians']!=b['jacobians'] or a['ok']!=b['ok']:
        print(r['id'],a,b)
raise SystemExit(bool(regressions))
