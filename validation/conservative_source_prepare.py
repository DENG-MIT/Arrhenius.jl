"""Prepare exact pinned Cantera example parameters and coefficient sidecars only."""
import argparse,ast,hashlib,json,shutil,sys
from pathlib import Path
import numpy as np
import cantera as ct
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("project",type=Path);p.add_argument("cantera_source",type=Path);p.add_argument("output",type=Path)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(a.project/"mechanism"))
from export_sidecar import export
from export_multicomponent import export_multicomponent
metadata=dict(cantera_version=ct.__version__,source_commit="726522be4e2a13454d8415b7ef799d621f665cf3",mechanisms={},examples={})
for name in ["h2o2","gri30"]:
    source=a.cantera_source/"data"/(name+".yaml");target=a.output/(name+".yaml")
    shutil.copyfile(source,target)
    export(target,Path(str(target)+".npz"))
    export_multicomponent(str(target),str(target)+".multicomponent.npz")
    g=ct.Solution(str(target))
    metadata["mechanisms"][name]=dict(sha256=hashlib.sha256(target.read_bytes()).hexdigest(),species=g.n_species,reactions=g.n_reactions)
for name in ["adiabatic_flame","burner_flame","flame_fixed_T"]:
    path=a.cantera_source/"samples/python/onedim"/(name+".py")
    metadata["examples"][name]=dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    if name=="flame_fixed_T":
        values={}
        for node in ast.parse(path.read_text()).body:
            if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
                key=node.targets[0].id
                if key in ("zloc","tvalues"):
                    values[key]=np.asarray(ast.literal_eval(node.value.args[0]),dtype=float)
        np.savez(a.output/"fixed-profile.npz",positions=values["zloc"],temperatures=values["tvalues"])
(a.output/"provenance.json").write_text(json.dumps(metadata,indent=2))
print(json.dumps(metadata,indent=2))
