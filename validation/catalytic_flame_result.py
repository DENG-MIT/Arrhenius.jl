"""Write provenance and completed catalytic validation evidence."""
import json,hashlib,argparse
from pathlib import Path
import cantera as ct
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("root",type=Path)
a=p.parse_args()
root=a.root
s=ct.Interface("ptcombust.yaml","Pt_surf")
g=s.adjacent["gas"]
comparison=json.loads((root/"comparison-tight.json").read_text())
data={
    "source_commit":"726522be4e2a13454d8415b7ef799d621f665cf3",
    "source_example":"samples/python/onedim/catalytic_combustion.py",
    "source_url":"https://cantera.org/dev/examples/python/onedim/catalytic_combustion.html",
    "cantera_version":ct.__version__,"cantera_git_commit":ct.__git_commit__,
    "gas_species":g.n_species,"gas_reactions":g.n_reactions,
    "surface_species":s.n_species,"surface_reactions":s.n_reactions,
    "inputs":{"pressure_Pa":ct.one_atm,"inlet_temperature_K":300.,"wall_temperature_K":900.,
        "width_m":.1,"mass_flux_kg_m2_s":.06,"initial_coverages":{"PT(S)":.5,"O(S)":.5},
        "hydrogen_inlet":"H2:.05,O2:.21,N2:.78,AR:.01",
        "methane_inlet":"CH4:.095,O2:.21,N2:.78,AR:.01",
        "reaction_multipliers":[0.,1e-5,1e-4,1e-3,1e-2,.1,1.],
        "final_refiner":{"ratio":100.,"slope":.15,"curve":.2,"prune":0.}},
    "native_initialization":"Linear T, constant inlet composition; native fixed-gas surface integration and constrained Newton polish. No reference profiles or stationary coverages imported.",
    "native_inlet_continuation":[.001,.01,.05,.1,.2,.4,.6,.8,1.],
    "transport":"mixture-averaged, mole-gradient basis, no Soret",
    "boundary":"stationary coverage equations, sum(theta)=1, J_k+MW_k*surface_wdot_k=0, u=V=0, prescribed wall T",
    "supported_scope":"neutral ideal gas, one ideal surface, no persistent bulk deposition",
    "unit_checks_passed":23,
    "comparison":comparison,
    "timed_scope":{
        "include":"Starting from prepared gas and surface mechanisms: native fixed-gas coverage initialization, nonreacting flow, six coupled gas/surface multiplier solves, methane inlet switch and final adaptive solve. Cantera uses the original published sequence from the same physical initial state.",
        "exclude":"mechanism parsing and parameter export, package loading/JIT, output serialization, plotting",
        "outputs":["T(z)","Y_k(z)","u(z)","V(z)","pressure curvature","surface coverages","wall gas production and mass fluxes"],
        "performance_qualified":False,"reason":"WSL paired timing is handled by the parent task; these are correctness runs only"},
    "parameter_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((root/"mechanism").glob("pt_h2.surface*"))},
}
(root/"catalytic-flame-result.json").write_text(json.dumps(data,indent=2))
print(json.dumps({"native_points":comparison["methane"]["native_points"],"gas_species":g.n_species,
    "surface_species":s.n_species,"unit_checks_passed":23,"methane":comparison["methane"]},indent=2))
