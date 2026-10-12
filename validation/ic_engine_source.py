"""Pinned engine source calculations and common comparison primitives.

Extracted unchanged from the existing timing driver so timing and the corrected
accepted-history reference share these operations without circular imports.
"""
from pathlib import Path
import ast
import copy
import hashlib
import numpy as np
import cantera as ct
SOURCE_COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
SOURCE_SHA256 = "43acf803aa5e4589bb1e718a480cbcc7db4558e49cdb97bb5b7ff67d928a68bb"
REFINED_RTOL, REFINED_ATOL = 1e-13, 1e-24
TRAPEZOID = getattr(np,"trapezoid",None) or np.trapz


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def decode_metadata(archive):
    return {key[:-5]:bytes(archive[key]).decode() for key in archive.files if key.endswith("_utf8")}


def source_programs(path,mechanism):
    """Select original calculation statements without rewriting its physics."""
    if sha(path) != SOURCE_SHA256:
        raise ValueError("ic_engine.py does not match the pinned source SHA256")
    tree = ast.parse(Path(path).read_text())
    setup,loop = [],None
    for original in tree.body:
        node = copy.deepcopy(original)
        if isinstance(node,(ast.Import,ast.ImportFrom)) or (isinstance(node,ast.Expr) and
                isinstance(node.value,ast.Constant) and isinstance(node.value.value,str)):
            continue
        if isinstance(node,ast.FunctionDef) and node.name=="ca_ticks":
            break
        if isinstance(node,ast.While):
            loop = node
            continue
        if loop is None:
            if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=="reaction_mechanism" for t in node.targets):
                node.value = ast.Constant(str(mechanism))
            setup.append(node)
    # Integral assignments occur after plots, so find their section separately.
    integrals=[]
    in_integrals=False
    for original in tree.body:
        if isinstance(original,ast.Assign) and any(isinstance(t,ast.Name) and t.id=="Q" for t in original.targets):
            in_integrals=True
        if in_integrals and isinstance(original,(ast.Assign,ast.AugAssign)):
            integrals.append(copy.deepcopy(original))
    if loop is None or not integrals:
        raise ValueError("pinned source calculation/integral sections were not found")
    def code(nodes):
        return compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),str(path),"exec")
    return code(setup),code([loop]),code(ast.parse("t = states.t").body+integrals)


def new_source_state(programs):
    namespace={"ct":ct,"np":np,"trapezoid":TRAPEZOID}
    exec(programs[0],namespace)
    return namespace


def source_calculation(programs):
    namespace=new_source_state(programs)
    exec(programs[1],namespace)
    exec(programs[2],namespace)
    states=namespace["states"]
    return {
        "time":states.t.copy(),"crank_angle":states.ca.copy(),"temperature":states.T.copy(),
        "pressure":states.P.copy(),"volume":states.V.copy(),"mass":states.m.copy(),
        "entropy_mass":states.entropy_mass.copy(),"mean_molecular_weight":states.mean_molecular_weight.copy(),
        "mdot_in":states.mdot_in.copy(),"mdot_out":states.mdot_out.copy(),
        "work_rate":states.dWv_dt.copy(),"heat_release_rate":(states.heat_release_rate*states.V).copy(),
        "X":states.X.copy(),"Y":states.Y.copy(),
        "integrals":dict(heat_J=float(namespace["Q"]),work_J=float(namespace["W"]),
                         efficiency=float(namespace["eta"]),CO_ppm=float(1e6*namespace["CO_emission"])),
        "solver_stats":namespace["sim"].solver_stats,
        "rtol":namespace["sim"].rtol,"atol":namespace["sim"].atol,
    }


def output_replay(actual,expected):
    return all(np.array_equal(actual[key],expected[key]) for key in expected
               if isinstance(expected[key],np.ndarray)) and actual["integrals"]==expected["integrals"] \
        and actual["solver_stats"]==expected["solver_stats"] and actual["rtol"]==expected["rtol"] \
        and actual["atol"]==expected["atol"]


def source_checks(result):
    return (result["rtol"]==1e-12 and result["atol"]==1e-16 and len(result["time"])>2880
        and .16 <= result["time"][-1] < .16+1/(360*50)
        and all(np.isfinite(value).all() for value in result.values() if isinstance(value,np.ndarray))
        and np.all(np.diff(result["time"])>0) and np.min(result["Y"])>=-1e-12)


def source_snapshot(namespace,time_offset=0.):
    gas=namespace["cyl"].phase
    cylinder=namespace["cyl"]
    t=namespace["sim"].time+time_offset
    work=-(gas.P-namespace["ambient_air"].phase.P)*namespace["A_piston"]*namespace["piston_speed"](t)
    return dict(time=t,temperature=gas.T,pressure=gas.P,volume=cylinder.volume,mass=cylinder.mass,
        entropy_mass=gas.entropy_mass,mean_molecular_weight=gas.mean_molecular_weight,
        mdot_in=namespace["inlet_valve"].mass_flow_rate,mdot_out=namespace["outlet_valve"].mass_flow_rate,
        mdot_fuel=namespace["injector_mfc"].mass_flow_rate,work_rate=work,
        heat_release_rate=gas.heat_release_rate*cylinder.volume,CO_X=gas["co"].X[0])


def reference_segments(namespace):
    """Independent tight reference on continuous, local-time source segments."""
    network=namespace["sim"]
    network.rtol,network.atol=REFINED_RTOL,REFINED_ATOL
    network.max_steps=1000000
    network.max_time_step=1/(360*50)
    stops=np.unique(np.r_[0.,[(720*cycle+angle)/(360*50) for cycle in range(4)
        for angle in (18,198,350,365,522,702)],.16])
    for start,stop in zip(stops[:-1],stops[1:]):
        midpoint=(start+stop)/2
        for device,opening,delta in (("inlet_valve","inlet_open","inlet_delta"),
                ("outlet_valve","outlet_open","outlet_delta"),
                ("injector_mfc","injector_open","injector_delta")):
            value=np.mod(namespace["crank_angle"](midpoint)-namespace[opening],4*np.pi)<namespace[delta]
            namespace[device].time_function=lambda t,value=value:value
        namespace["piston"].velocity=lambda t,start=start:namespace["piston_speed"](t+start)
        network.initial_time=0.
        network.initialize()
        yield float(start),float(stop)


def integral_terms(output,indices=None):
    ix=np.arange(len(output["time"])) if indices is None else indices
    t=output["time"][ix]
    heat=TRAPEZOID(output["heat_release_rate"][ix],t)
    work=TRAPEZOID(output["work_rate"][ix],t)
    weights=output["mean_molecular_weight"][ix]*output["mdot_out"][ix]
    return np.array([heat,work,TRAPEZOID(weights*output["CO_X"][ix],t),TRAPEZOID(weights,t)])


def integral_values(terms):
    heat,work,numerator,denominator=terms
    return dict(heat_J=float(heat),work_J=float(work),efficiency=float(work/heat),CO_ppm=float(1e6*numerator/denominator))


def compare_native(native,summary,reference,reference_integrals):
    actual={key[len("output_"):]:native[key] for key in native.files if key.startswith("output_")}
    errors={
        "temperature_K":float(np.max(np.abs(actual["temperature"]-reference["temperature"]))),
        "pressure_relative":float(np.max(np.abs(actual["pressure"]/reference["pressure"]-1))),
        "volume_m3":float(np.max(np.abs(actual["volume"]-reference["volume"]))),
        "mass_relative":float(np.max(np.abs(actual["mass"]/reference["mass"]-1))),
        "mass_fraction":float(np.max(np.abs(actual["Y"]-reference["Y"]))),
        "entropy_J_per_kg_K":float(np.max(np.abs(actual["entropy_mass"]-reference["entropy_mass"]))),
    }
    limits=dict(temperature_K=.5,pressure_relative=2e-4,volume_m3=2e-10,mass_relative=5e-5,
                mass_fraction=1e-4,entropy_J_per_kg_K=.5)
    switches=np.array([(720*cycle+angle)/(360*50) for cycle in range(4)
        for angle in (18,198,350,365,522,702)])
    continuous=np.min(np.abs(actual["time"][:,None]-switches),axis=1)>1e-11
    rate_errors={key:float(np.max(np.abs(actual[key][continuous]-reference[key][continuous]))/
        max(np.max(np.abs(reference[key][continuous])),1e-30))
        for key in ("mdot_in","mdot_out","mdot_fuel","work_rate","heat_release_rate")}
    rate_limits=dict(mdot_in=2e-3,mdot_out=2e-3,mdot_fuel=1e-13,work_rate=2e-3,heat_release_rate=5e-3)
    integral_errors={key:abs(summary["integrals"][key]/value-1) for key,value in reference_integrals.items()}
    integral_limits=dict(heat_J=1e-4,work_J=1e-5,efficiency=1e-4,CO_ppm=1e-4)
    return dict(correctness_pass=all(errors[k]<v for k,v in limits.items()) and
        all(integral_errors[k]<v for k,v in integral_limits.items()) and
        all(rate_errors[k]<v for k,v in rate_limits.items()),trajectory_errors=errors,
        trajectory_limits=limits,rate_peak_scaled_errors=rate_errors,rate_limits=rate_limits,
        integral_relative_errors=integral_errors,integral_limits=integral_limits)


def verify_mechanisms(source,native_path,native_meta):
    if sha(native_path)!=native_meta["mechanism_sha256"] or sha(str(native_path)+".npz")!=native_meta["sidecar_sha256"]:
        raise ValueError("native mechanism/sidecar files differ from the timed Julia artifact")
    sidecar_meta=decode_metadata(np.load(str(native_path)+".npz"))
    if sidecar_meta.get("source_sha256")!=sha(native_path):
        raise ValueError("native sidecar provenance does not match its YAML")
    original,prepared=ct.Solution(str(source),"nDodecane_IG"),ct.Solution(str(native_path))
    if original.species_names!=prepared.species_names or original.n_species!=100 or not original.n_reactions==prepared.n_reactions==553:
        raise ValueError("source/prepared phases do not contain the required 100 species and 553 reactions")
    for T in (300.,1000.,2500.):
        original.TPX=prepared.TPX=T,1.3e5,"o2:1,n2:3.76"
        for key in ("molecular_weights","standard_enthalpies_RT","standard_cp_R","standard_entropies_R",
                    "forward_rate_constants","reverse_rate_constants"):
            np.testing.assert_allclose(getattr(original,key),getattr(prepared,key),rtol=2e-12,atol=1e-20,
                                       err_msg=f"prepared phase differs in {key} at {T} K")
    return sidecar_meta
