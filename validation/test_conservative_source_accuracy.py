"""Ensure incomplete or failed reference comparisons cannot qualify a source run."""
import argparse,json,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
import cantera as ct
p=argparse.ArgumentParser(description=__doc__);p.add_argument("mechanism",type=Path);a=p.parse_args()
g=ct.Solution(str(a.mechanism));g.TPX=373.,.05*ct.one_atm,"H2:1.5,O2:1,AR:7"
script=Path(__file__).with_name("conservative_source_accuracy.py")
with tempfile.TemporaryDirectory() as temporary:
    root=Path(temporary);native=root/"native";reference=root/"reference";native.mkdir();reference.mkdir()
    report=root/"report.json"
    def write(path,level):
        z=np.linspace(0.,.5,4*2**level+1)
        np.savez(path,grid=z,T=373.+1000*z,velocity=np.ones(len(z)),Y=np.repeat(g.Y[:,None],len(z),axis=1),inlet_Y=g.Y)
    def run(expected):
        result=subprocess.run([sys.executable,str(script),str(a.mechanism),str(native),str(reference),"burner",str(report)],capture_output=True,text=True)
        data=json.loads(report.read_text())
        assert (result.returncode==0)==expected,(result.returncode,result.stdout,result.stderr)
        assert data["passed"]==expected,data
        assert set(data["results"])=={"mole","multi"},data
    for mode in ["mole","multi"]:
        write(native/f"burner-{mode}-0.npz",0)
        for level in range(3):write(reference/f"burner-{mode}-{level}.npz",level)
    run(True)
    (native/"burner-multi-0.npz").unlink();run(False);write(native/"burner-multi-0.npz",0)
    (reference/"burner-mole-1.npz").unlink();run(False);write(reference/"burner-mole-1.npz",1)
    malformed=reference/"burner-mole-invalid.npz";malformed.write_bytes(b"invalid level");run(False);malformed.unlink()
    for mode in ["mole","multi"]:
        for level in range(3,11):write(reference/f"burner-{mode}-{level}.npz",level)
    run(True) # Level 10 must follow level 9, not level 1.
    file=native/"burner-mole-0.npz";v=np.load(file);bad={key:v[key] for key in v.files};v.close()
    bad["Y"][g.species_index("H2")]+=.001;bad["Y"][g.species_index("AR")] -= .001
    np.savez(file,**bad);run(False)
    write(file,0);bad["Y"]=bad["Y"][:,:-1];np.savez(file,**bad);run(False)
print("Source accuracy validator: 7 positive/negative integrity checks passed")
