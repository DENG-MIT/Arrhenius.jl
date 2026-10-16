"""Check the complete stationary mixer outputs against independent source solves."""
import argparse
import hashlib
import json
from pathlib import Path
import tomllib
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare(native_dir, reference_dir, destination):
    n=np.load(native_dir/"native.npz",allow_pickle=False)
    meta=tomllib.loads((native_dir/"native.toml").read_text())
    original=np.load(reference_dir/"original.npz",allow_pickle=False)
    refined=np.load(reference_dir/"refined.npz",allow_pickle=False)
    if not meta["inputs_unchanged"] or meta["inputs_before"]!=meta["inputs_after"]:
        raise ValueError("native source/input identity check failed")
    if sha(native_dir/"native.npz")!=meta["native_sha256"]:
        raise ValueError("native archive digest mismatch")
    for key in ("runtime_before","runtime_after"):
        if any(value!=1 for value in meta[key].values()):
            raise ValueError("native numerical thread check failed")
    reference=json.loads((reference_dir/"reference.json").read_text())
    # Fixed physical/output tolerances. Near-zero balances use throughput scales.
    limits={name:(2e-7,1e-7) for name in ("pressure","density","mass","mean_molecular_weight",
        "enthalpy","internal_energy","entropy","gibbs","cp","cv","flows")}
    limits.update(Y=(2e-7,1e-10),X=(2e-7,1e-10),temperature=(0.,2e-5))
    comparisons={}
    for label,left,right in (("reference_refinement",original,refined),
                             ("native_original",n,original),("native_refined",n,refined)):
        values={}
        for name,(rtol,atol) in limits.items():
            a,b=left[name],right[name]
            valid=a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
            error=np.abs(a-b)
            values[name]=dict(passed=bool(valid and np.all(error<=atol+rtol*np.abs(b))),
                maximum_absolute_error=float(np.max(error)),rtol=rtol,atol=atol)
        comparisons[label]=values
    mass=float(n["mass"][0]);temperature=float(n["temperature"][0])
    rate=float(max(np.max(np.abs(n["species_rate"]))/mass,abs(n["temperature_rate"][0])/temperature))
    energy_scale=max(abs(float(refined["flows"][2]*refined["enthalpy"][0])),1.)
    energy_residual=abs(float(n["external_energy_rate"][0]))/energy_scale
    element_scale=max(float(refined["flows"][2]/refined["mean_molecular_weight"][0]),1e-300)
    element_residual=float(np.max(np.abs(n["element_rates"]))/element_scale)
    initial=n["initial_state"]
    checks=dict(reference_passed=reference["passed"],converged=meta["converged"],
        chemistry_enabled=meta["chemistry_enabled"],separate_species_sets=meta["species_counts"]==[53,8],
        original_initial_state=bool(initial.shape==(54,) and initial[-1]==300.
            and np.allclose(initial[:-1]/np.sum(initial[:-1]),original["initial_Y"],rtol=2e-14,atol=0.)),
        all_outputs_match=all(row["passed"] for group in comparisons.values() for row in group.values()),
        physical_residual=rate<=1e-9,physical_species_domain=float(np.min(n["Y"]))>=-1e-13,
        mass_fraction_sum=abs(float(np.sum(n["Y"]))-1.)<=1e-12,
        steady_energy_balance=energy_residual<=1e-8,element_balance=element_residual<=1e-10,
        energy_ledger=bool(np.allclose(n["total_energy_rate"],n["external_energy_rate"],rtol=1e-12,atol=1e-12)))
    report=dict(passed=all(checks.values()),checks=checks,comparisons=comparisons,
        physical_residual_per_s=rate,energy_residual_relative_to_enthalpy_throughput=energy_residual,
        element_residual_relative_to_molar_throughput=element_residual,
        hashes={str(p):sha(p) for p in (native_dir/"native.npz",native_dir/"native.toml",
            reference_dir/"original.npz",reference_dir/"refined.npz",reference_dir/"reference.json",Path(__file__))},
        performance_qualified=False)
    destination.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k not in ("comparisons","hashes")},indent=2))
    if not report["passed"]:
        raise ValueError("complete mixer gates failed; report saved")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("native","reference","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();compare(a.native,a.reference,a.output)
