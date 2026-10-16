"""Apply existing flame physical gates to source runs and refined CT references."""
import argparse,json,re
from pathlib import Path
import numpy as np
import cantera as ct
from flame_benchmarks import profile_errors,element_conservation,SPEED_TOL,TMAX_TOL,PROFILE_TOL,PROFILE_FLOOR,ELEMENT_TOL
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("mechanism",type=Path);p.add_argument("native",type=Path);p.add_argument("reference",type=Path)
p.add_argument("case",choices=["free","burner","fixed"]);p.add_argument("output",type=Path)
a=p.parse_args();g=ct.Solution(str(a.mechanism))

def load(path):
    v=np.load(path);d={key:v[key] for key in v.files};d["species_names"]=g.species_names
    for key in ["grid","T","Y","velocity","inlet_Y"]:
        if key not in d or not np.all(np.isfinite(d[key])):raise ValueError(f"{path}: missing/nonfinite {key}")
    N=len(d["grid"])
    if N<5 or d["grid"].shape!=(N,) or not np.all(np.diff(d["grid"])>0):raise ValueError(f"{path}: invalid grid")
    if d["T"].shape!=(N,) or d["velocity"].shape!=(N,) or d["Y"].shape!=(g.n_species,N) or d["inlet_Y"].shape!=(g.n_species,):
        raise ValueError(f"{path}: inconsistent state/species dimensions")
    if np.min(d["T"])<200 or np.max(d["T"])>6000 or np.min(d["Y"]) < -1e-7:
        raise ValueError(f"{path}: invalid temperature or mass fractions")
    if not np.allclose(np.sum(d["Y"],axis=0),1.,rtol=0,atol=1e-7) or not np.isclose(sum(d["inlet_Y"]),1.,rtol=0,atol=1e-8):
        raise ValueError(f"{path}: unnormalized mass fractions")
    return d

def extrapolate(previous,last):
    za,zb=previous["grid"].copy(),last["grid"].copy()
    if a.case=="free":
        za-=np.interp(700.,previous["T"],za);zb-=np.interp(700.,last["T"],zb)
    lo=max(za[0],zb[0]);hi=min(za[-1],zb[-1])
    q=np.unique(np.r_[za,zb]);q=q[(q>=lo)&(q<=hi)]
    d=dict(grid=q,T=2*np.interp(q,zb,last["T"])-np.interp(q,za,previous["T"]),
        Y=np.stack([2*np.interp(q,zb,yb)-np.interp(q,za,ya) for ya,yb in zip(previous["Y"],last["Y"])]),
        velocity=2*np.interp(q,zb,last["velocity"])-np.interp(q,za,previous["velocity"]),
        species_names=g.species_names,inlet_Y=last["inlet_Y"])
    return d

results={};failures=[]
required_modes=["mass","mass-soret","multi","multi-soret"] if a.case=="free" else ["mole","multi"]
for mode in required_modes:
    try:
        path=a.native/f"{a.case}-{mode}-0.npz"
        if not path.is_file():raise ValueError(f"missing native source stage: {path}")
        indexed=[]
        for candidate in a.reference.glob(a.case+"-*.npz"):
            stem,suffix=candidate.stem.rsplit("-",1)
            if stem!=f"{a.case}-{mode}":continue
            if not re.fullmatch(r"0|[1-9][0-9]*",suffix):raise ValueError(f"malformed reference level: {candidate}")
            indexed.append((int(suffix),candidate))
        indexed.sort(key=lambda pair:pair[0])
        levels=[level for level,_ in indexed]
        if len(levels)<3 or levels!=list(range(levels[-1]+1)):
            raise ValueError(f"{mode}: at least three contiguous reference levels starting at zero required; got {levels}")
        references=[load(candidate) for _,candidate in indexed]
        if any(len(second["grid"])!=2*len(first["grid"])-1 or not np.allclose(second["grid"][::2],first["grid"],rtol=0,atol=1e-14)
                for first,second in zip(references,references[1:])):
            raise ValueError(f"{mode}: references are not successive whole-grid bisections")
        jul=load(path);old,previous,last=references[-3:]
        ref=extrapolate(previous,last);earlier=extrapolate(old,previous)
        errors=profile_errors(ref,jul,a.case);uncertainty=profile_errors(ref,earlier,a.case)
        speed_error=abs(jul["velocity"][0]-ref["velocity"][0])/abs(ref["velocity"][0])
        Tmax_error=abs(max(jul["T"])-max(ref["T"]))/abs(max(ref["T"]))
        elements=element_conservation(g,jul)
        gates=dict(species_profiles=errors["all_species_pass"],temperature_profile=errors["temperature"]["normalized_peak_error"]<=TMAX_TOL,
            peak_temperature=Tmax_error<=TMAX_TOL,elements=all(e["pass"] for e in elements.values()))
        if a.case=="free":gates["speed"]=speed_error<=SPEED_TOL
        results[mode]=dict(native_points=len(jul["grid"]),native_domain=[float(jul["grid"][0]),float(jul["grid"][-1])],
            reference_levels=levels,reference_points=[len(previous["grid"]),len(last["grid"])],reference_domain=[float(last["grid"][0]),float(last["grid"][-1])],
            profiles=errors,reference_extrapolation_change=uncertainty,speed_relative_error=float(speed_error),Tmax_relative_error=float(Tmax_error),
            elements=elements,gates=gates,passed=all(gates.values()))
        if not all(gates.values()):failures.append(f"{mode}: failed gates {[key for key,value in gates.items() if not value]}")
        print(a.case,mode,"PASS" if all(gates.values()) else "FAIL","worst species",max((v["normalized_peak_error"],k) for k,v in errors["species"].items()),flush=True)
    except Exception as error:
        results[mode]=dict(passed=False,error=str(error));failures.append(f"{mode}: {error}")
        print(a.case,mode,"FAIL",error,flush=True)
a.output.write_text(json.dumps(dict(reference="First-order CT extrapolation from controlled whole-grid bisection, with last-extrapolation change retained.",
    criteria=dict(speed=SPEED_TOL,temperature=TMAX_TOL,species=PROFILE_TOL,species_floor=PROFILE_FLOOR,elements=ELEMENT_TOL),
    required_modes=required_modes,results=results,failures=failures,passed=not failures),indent=2,
    default=lambda obj:obj.item() if isinstance(obj,np.generic) else obj.tolist()))
raise SystemExit(1 if failures else 0)
