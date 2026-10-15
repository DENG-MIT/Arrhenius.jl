"""Summarize the A/B policy batch using the existing branch and held-out gates."""
import csv, json, os, subprocess
import numpy as np
from pathlib import Path
root=Path.home()/"work/nh3-parity"
task=Path.home()/"work/hotwin"
out=task/"round2"
tool=task/"repo/test/cantera_parity"
assert "REF_DIR" not in os.environ
for name,target in (("main","A"),("branch","B")):
    p=out/name
    if not p.exists(): p.symlink_to(target,target_is_directory=True)
gate_run=subprocess.run([os.environ["TRACE_CAP"],"900",os.environ["PY"],str(tool/"summarize_branchcheck.py"),str(root),str(out)])
assert gate_run.returncode in (0,1), f"branch gate command failed: {gate_run.returncode}"
d=json.load(open(out/"summary.json"))
baseline=task/"final_branchcheck"
old={r["id"]:r for r in csv.DictReader(open(baseline/"main/counts.tsv"),delimiter="\t")}
baseline_equal=True
for row in d["results"]:
    cid=row["id"]
    baseline_equal &= all(row["counts"]["main"][k]==old[cid][k] for k in ("ok","jacobians","newton","steps","failures","points"))
    with np.load(out/"A"/(cid+".npz")) as a, np.load(baseline/"main"/(cid+".npz")) as b:
        baseline_equal &= all(a[k].shape==b[k].shape and a[k].dtype==b[k].dtype and a[k].tobytes()==b[k].tobytes() for k in ("grid","T","Y","velocity","SL"))
assert baseline_equal, "A no longer reproduces the verified main baseline"
flags=[]
changes=[]
for row in d["results"]:
    for k in ("jacobians","newton","steps","failures","points"):
        a=int(row["counts"]["main"][k]); b=int(row["counts"]["branch"][k])
        if b!=a:
            r={"id":row["id"],"counter":k,"A":a,"B":b,"ratio":b/a if a else None}
            changes.append(r)
            if b>1.25*a: flags.append(r)
for label in ("A","B"):
    folder=out/label
    case_id=(task/"H09_only.txt").read_text().strip()
    for target,source in ((folder/(case_id+".ct.npz"),task/"evidence/H09.ct.npz"),
                          (folder/(case_id+".native.npz"),folder/"H09.npz")):
        if not target.exists(): target.symlink_to(source)
    subprocess.run([os.environ["TRACE_CAP"],"120",os.environ["PY"],str(tool/"validate_heldout_subset.py"),str(root),str(task/"H09_only.txt"),str(folder),str(folder/"heldout_gate")],check=True)
result={"cases":d["cases"],"A_pass":d["main_pass"],"B_pass":d["branch_pass"],
        "A_matches_original_main_counts_and_profiles":bool(baseline_equal),
        "accuracy_regressions":d["accuracy_regressions"],"growth_flags":flags,"count_changes":changes,
        "selected":"A" if flags or d["accuracy_regressions"] else "B",
        "H09":{label:next(r for r in csv.DictReader(open(out/label/"counts.tsv"),delimiter="\t") if r["id"]=="H09") for label in ("A","B")}}
with np.load(out/"B/H09.npz") as b, np.load(task/"final_H09/H09.native.npz") as retained:
    result["H09_B_matches_retained_profiles"]=all(b[k].shape==retained[k].shape and b[k].dtype==retained[k].dtype and b[k].tobytes()==retained[k].tobytes() for k in ("grid","T","Y","velocity","SL"))
    result["H09_sumY_max"]=float(np.max(np.abs(np.sum(b["Y"],axis=0)-1)))
assert result["H09_B_matches_retained_profiles"]
(out/"decision.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result,indent=2))
