"""Pair the full native shock-tube calculation with Cantera 4 on one host.

Run real_gas_timing.jl first, then this file with --julia-result and --output.
Mechanism loading, imports, file I/O, plots and independent refined-reference
checks are outside measured loops. Both timed calculations initialize each
reactor, take adaptive steps until time >= 0.005 s, save every 20th step and
extract the OH peak. Seven warm repetitions are required. --validate-only
executes one untimed-for-qualification calculation and the correctness checks.
Controlled timing also requires the pinned --source-example, a
--cantera-build-record and --native-mechanism-directory. The native file hashes,
sidecar association, and complete species/reaction data are checked before timing.

The published default Cantera tolerances are retained in its timed loop. A
separate tight reference checks native trajectories: at 760 K the published
default and converged Cantera temperatures differ by about 26 K. This is
reported explicitly rather than weakening the native trajectory criterion.

Source: https://github.com/Cantera/cantera/blob/726522be4e2a13454d8415b7ef799d621f665cf3/samples/python/reactors/non_ideal_shock_tube.py
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import os
import statistics
import time
import tomllib
from collections.abc import Mapping
from numbers import Real
from benchmark_environment import host_metadata, matches_target, loaded_library_paths, verify_numerical_threads

THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
for name in THREAD_VARIABLES:
    os.environ[name] = "1"
started = time.perf_counter()
import numpy as np
import cantera as ct
import cantera._cantera as compiled
IMPORT_SECONDS = time.perf_counter()-started
SOURCE_COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
# Byte hashes of the files at SOURCE_COMMIT (LF checkout).
SOURCE_EXAMPLE_SHA256 = "88f73fd9a1acf122ba1f71789d974d221b133f82c18381615c8099ff84b2f1de"
SOURCE_MECHANISM_SHA256 = "3d3b59ed91dec0d0bcbac2fa2ef2cba13fbd565bfff8f266ca847ed6aa92f7f1"
TEMPERATURES = [1250,1170,1120,1080,1040,1010,990,970,950,930,910,880,850,820,790,760]
EXAMPLE_HELPERS = ("non_ideal_shock_tube.jl","real_gas_ode_solver.jl","real_gas_ad_jacobian.jl",
                   "real_gas_trial_states.jl","real_gas_qndf_solver.jl")
NATIVE_MECHANISMS = ("dodecane_RK.yaml","dodecane_IG.yaml","dodecane_IG.yaml.npz")


def case_order(smoke=False):
    return [("RK",1000),("IG",1000)] + ([] if smoke else
        [(phase,temperature) for phase in ("RK","IG") for temperature in TEMPERATURES])


def capture_phase_states(gases):
    """Save the freshly constructed phases before any whole-example call."""
    return {phase:gas.state.copy() for phase,gas in gases.items()}


def restore_phase_states(gases,states):
    """Restore whole-example inlet history outside the calculation timer."""
    if gases.keys()!=states.keys():
        raise ValueError("prepared phase membership differs from the calculation")
    for phase,gas in gases.items():
        gas.state = states[phase]


def make_network(gas,temperature,refined=False):
    gas.TP = temperature,40*ct.one_atm
    gas.set_equivalence_ratio(1.,"c12h26",{"o2":1.,"n2":3.76})
    reactor = ct.IdealGasMoleReactor(gas,clone=False)
    network = ct.ReactorNet([reactor])
    network.preconditioner = ct.AdaptivePreconditioner()
    if refined:
        network.rtol,network.atol = 1e-11,1e-22
    return reactor,network


def calculate(gases,smoke=False,refined=False):
    results = {}
    for phase,temperature in case_order(smoke):
        gas = gases[phase]
        reactor,network = make_network(gas,temperature,refined=refined)
        history = ct.SolutionArray(gas,extra=["t"])
        counter = 0
        while network.time < .005:
            network.step()
            counter += 1
            if counter % 20 == 0:
                history.append(reactor.phase.state,t=network.time)
        index = int(history("oh").Y.argmax())
        results[f"{phase}_{temperature}"] = {
            "delay":float(history.t[index]),"time":history.t.copy(),
            "state":np.vstack([history.Y.T,history.T]),
            "steps":counter,"final_time":network.time,"final_state":np.r_[gas.Y,gas.T],
            "rtol":network.rtol,"atol":network.atol,"solver_stats":network.solver_stats}
    return results


def validate_case(gas,temperature,key,native,published):
    _,network = make_network(gas,temperature,refined=True)
    saved_times = native[key+"_time"]
    # Check the final overshooting endpoint as well as the source's saved grid.
    times = np.r_[saved_times,native[key+"_final_time"]]
    states = np.column_stack([native[key+"_state"],native[key+"_final_state"]])
    initial_energy = gas.int_energy_mass
    element_matrix = np.array([[gas.n_atoms(k,e) for k in range(gas.n_species)]
        for e in range(gas.n_elements)]) / gas.molecular_weights
    initial_elements = element_matrix @ gas.Y
    all_times = np.unique(np.r_[times,published["time"]])
    all_reference = np.empty((states.shape[0],len(all_times)))
    all_pressure = np.empty_like(all_times)
    for i,t in enumerate(all_times):
        network.advance(float(t))
        all_reference[:,i] = np.r_[gas.Y,gas.T]
        all_pressure[i] = gas.P
    indices = np.searchsorted(all_times,times)
    reference,pressure = all_reference[:,indices],all_pressure[indices]
    published_reference = all_reference[:,np.searchsorted(all_times,published["time"])]
    dT = float(np.max(np.abs(states[-1]-reference[-1])))
    dY = float(np.max(np.abs(states[:-1]-reference[:-1])))
    oh = gas.species_index("oh")
    expected_peak = int(np.argmax(reference[oh,:len(saved_times)]))
    actual_peak = int(np.argmax(states[oh,:len(saved_times)]))
    density = float(native[key+"_density"][0])
    native_pressure = np.empty_like(times)
    energy = np.empty_like(times)
    for i in range(len(times)):
        gas.TDY = float(states[-1,i]),density,np.maximum(states[:-1,i],0.)
        native_pressure[i],energy[i] = gas.P,gas.int_energy_mass
    dP = float(np.max(np.abs(native_pressure/pressure-1)))
    mass_drift = float(np.max(np.abs(np.sum(states[:-1],axis=0)-1)))
    minimum_mass_fraction = float(np.min(states[:-1]))
    energy_drift = float(np.max(np.abs(energy-initial_energy))/max(abs(initial_energy),1e6))
    element_drift = float(np.max(np.abs(element_matrix @ states[:-1]-initial_elements[:,None])))
    # The independent-reference OH peak must occupy the same native sample or
    # one adjacent sample (a peak may fall almost halfway between sample times).
    peak_pass = abs(expected_peak-actual_peak)<=1
    source_delay = published["delay"]
    delay = float(native[key+"_delay"][0])
    report = {"correctness_pass":bool(dT<.3 and dY<2e-5 and dP<1e-4
                and mass_drift<5e-10 and element_drift<5e-11 and energy_drift<1e-6
                and minimum_mass_fraction>-1e-12 and peak_pass),
            "temperature_error_K":dT,"mass_fraction_error":dY,
            "pressure_relative_error":dP,"mass_drift":mass_drift,
            "minimum_mass_fraction":minimum_mass_fraction,
            "energy_relative_drift":energy_drift,"element_drift":element_drift,
            "source_default_temperature_error_K":float(np.max(np.abs(published["state"][-1]-published_reference[-1]))),
            "source_default_mass_fraction_error":float(np.max(np.abs(published["state"][:-1]-published_reference[:-1]))),
            "native_ignition_delay_s":delay,"published_ignition_delay_s":source_delay,
            "published_delay_relative_difference":abs(delay/source_delay-1),
            "refined_peak_on_native_grid_s":float(times[expected_peak]),
            "refined_peak_sample_distance":abs(expected_peak-actual_peak),
            "native_steps":int(native[key+"_steps"][0]),"cantera_steps":published["steps"],
            "native_final_time_s":float(native[key+"_final_time"][0]),
            "native_saved_samples":len(saved_times),"endpoint_checked":True,
            "cantera_final_time_s":published["final_time"],
            "cantera_solver_stats":published["solver_stats"],
            "native_rhs_evaluations":int(native[key+"_rhs_evaluations"][0]),
            "native_jacobian_evaluations":int(native[key+"_jacobian_evaluations"][0]),
            "reference_rtol":network.rtol,"reference_atol":network.atol}
    report["native_solver_stats"] = {field:int(native[key+"_"+field][0])
        for field in ("accepted_steps","rejected_steps","linear_solves","matrix_updates",
            "nonlinear_iterations","nonlinear_convergence_failures") if key+"_"+field in native}
    return report,reference


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def current_source_hashes(project):
    project = Path(project)
    paths = [p for p in (project/"src").rglob("*") if p.is_file()]
    paths += [project/"example/reactors"/name for name in EXAMPLE_HELPERS]
    paths.append(project/"validation/real_gas_timing.jl")
    return {p.relative_to(project).as_posix():sha(p) for p in sorted(paths)}


def verify_native_source_inventory(project,meta):
    before = tomllib.loads(meta.get("source_hashes_toml",""))
    after = tomllib.loads(meta.get("source_hashes_after_toml",""))
    current = current_source_hashes(project)
    if not before or before!=after or before!=current:
        raise ValueError("native source inventory or bytes differ before/after/current checkout")
    return current


def require_unchanged(expected):
    if any(not Path(path).is_file() or sha(path)!=value for path,value in expected.items()):
        raise ValueError("benchmark input or driver bytes changed during calculations")


def mapped_cantera_hashes(record=None):
    paths = {p.resolve() for p in loaded_library_paths() if "cantera" in p.name.lower()}
    extension = Path(compiled.__file__).resolve()
    if extension not in paths or not any("libcantera" in p.name for p in paths):
        raise ValueError("the Cantera extension and shared library must both be actually mapped")
    hashes = {str(p):sha(p) for p in sorted(paths)}
    if record is not None:
        recorded = {Path(k).name:v for k,v in record.get("library_hashes",{}).items()}
        if any(recorded.get(Path(p).name)!=value for p,value in hashes.items()):
            raise ValueError("actually mapped Cantera library/extension bytes differ from build record")
    return hashes


def same_trajectories(actual,expected):
    """Check complete deterministic replay outside the timed calculation."""
    return actual.keys()==expected.keys() and all(
        all(np.array_equal(actual[key][field],expected[key][field])
            for field in ("time","state","final_state"))
        and all(actual[key][field]==expected[key][field]
            for field in ("delay","steps","final_time"))
        for key in expected)


def require_replay(actual,expected,label):
    if not same_trajectories(actual,expected):
        raise ValueError(f"{label} differs from the checked first trajectory")
    return True


def native_replay_checks(native):
    samples = native["warm_seconds"].tolist()
    matches = native["warm_matches_first"].tolist() if "warm_matches_first" in native else []
    if len(samples)!=len(matches) or any(value!=1 for value in matches):
        raise ValueError("native warm repetitions lack successful full-trajectory replay checks")
    if any(not np.isfinite(value) or value<=0 for value in samples):
        raise ValueError("native warm durations must be finite and positive")
    return samples,[bool(value) for value in matches]


def same_input_data(actual,expected):
    """Allow only YAML writer roundoff when comparing complete species/rate data."""
    if isinstance(expected,Mapping):
        return (isinstance(actual,Mapping) and actual.keys()==expected.keys()
            and all(same_input_data(actual[key],value) for key,value in expected.items()))
    if isinstance(expected,(list,tuple,np.ndarray)):
        return (isinstance(actual,(list,tuple,np.ndarray)) and len(actual)==len(expected)
            and all(same_input_data(a,b) for a,b in zip(actual,expected)))
    if isinstance(expected,Real) and isinstance(actual,Real):
        return bool(np.isclose(actual,expected,rtol=1e-13,atol=0.))
    return actual==expected


def verify_native_mechanisms(directory,meta,gases):
    if directory is None:
        return {"verified":False,"reason":"native mechanism directory was not supplied"}
    names = NATIVE_MECHANISMS
    recorded = tomllib.loads(meta["mechanism_hashes_toml"])
    after = tomllib.loads(meta.get("mechanism_hashes_after_toml",""))
    current = {name:sha(directory/name) for name in names}
    if recorded!=after or recorded!=current:
        raise ValueError("native mechanism files differ from the files hashed by the Julia run")
    with np.load(directory/names[2]) as sidecar:
        sidecar_source = bytes(sidecar["source_sha256_utf8"]).decode()
    if sidecar_source!=current["dodecane_IG.yaml"]:
        raise ValueError("native kinetic sidecar does not match its ideal-gas YAML")
    for phase in ("RK","IG"):
        exported = ct.Solution(str(directory/f"dodecane_{phase}.yaml"))
        reference = gases[phase]
        if (exported.thermo_model!=reference.thermo_model
            or exported.species_names!=reference.species_names
            or exported.n_reactions!=reference.n_reactions):
            raise ValueError(f"native {phase} phase differs from the pinned source phase")
        for kind,actual,expected in (("species",exported.species(),reference.species()),
                                     ("reaction",exported.reactions(),reference.reactions())):
            for i,(a,b) in enumerate(zip(actual,expected)):
                if not same_input_data(a.input_data,b.input_data):
                    raise ValueError(f"native {phase} {kind} {i} differs from the pinned source data")
    return {"verified":True,"files_sha256":current,"sidecar_source_sha256":sidecar_source,
        "reference_comparison":"complete species thermo/EOS and reaction input data; YAML roundoff rtol=1e-13, atol=0"}


def verify_threads(enforce=True):
    result = verify_numerical_threads(set_accelerate=enforce)
    result["environment"] = {name:os.environ[name] for name in THREAD_VARIABLES}
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--julia-result",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--repetitions",type=int,default=7)
    parser.add_argument("--qualification",choices=("informational","controlled"),default="informational")
    parser.add_argument("--benchmark-target",choices=("wsl","apple-m4"),default="wsl")
    parser.add_argument("--cantera-build-record",type=Path)
    parser.add_argument("--source-example",type=Path)
    parser.add_argument("--mechanism",type=Path,help="defaults to Cantera's installed nDodecane_Reitz.yaml")
    parser.add_argument("--native-mechanism-directory",type=Path,
        help="prepared native YAML/sidecar directory; required for controlled timing")
    parser.add_argument("--native-project",type=Path,default=Path(__file__).resolve().parent.parent,
        help="current native source checkout; must exactly match the Julia before/after inventory")
    parser.add_argument("--validate-only",action="store_true")
    args = parser.parse_args()
    if args.repetitions<7:
        parser.error("at least seven warm repetitions are required")
    native = np.load(args.julia_result)
    meta = {k[:-5]:bytes(native[k]).decode() for k in native.files if k.endswith("_utf8")}
    if meta.get("scope") not in ("1000_K_pair","full_34_trajectories"):
        raise ValueError("unknown native shock-tube scope")
    smoke = meta["scope"]=="1000_K_pair"
    native_samples,native_warm_matches_first = native_replay_checks(native)
    source_hashes = verify_native_source_inventory(args.native_project,meta)
    if args.source_example is not None and sha(args.source_example)!=SOURCE_EXAMPLE_SHA256:
        raise ValueError("source example hash differs from the pinned Cantera source")
    if args.qualification=="controlled":
        if (args.source_example is None or args.native_mechanism_directory is None
            or args.cantera_build_record is None):
            parser.error("controlled timing requires --source-example, --native-mechanism-directory and --cantera-build-record")
        if smoke or args.validate_only or meta.get("qualification")!="controlled" or len(native_samples)<7:
            raise ValueError("controlled timing requires a full native controlled run with at least seven checked repetitions")
    order = [f"{phase}_{T}" for phase,T in case_order(smoke)]
    if meta["case_order"].split(",")!=order:
        raise ValueError("native case ordering differs from the source")
    mechanism = args.mechanism or next((Path(folder)/"nDodecane_Reitz.yaml"
        for folder in ct.get_data_directories() if (Path(folder)/"nDodecane_Reitz.yaml").is_file()),None)
    if mechanism is None:
        raise FileNotFoundError("supply the source nDodecane_Reitz.yaml with --mechanism")
    if sha(mechanism)!=SOURCE_MECHANISM_SHA256:
        raise ValueError("reference mechanism hash differs from the pinned Cantera source")
    gases = {phase:ct.Solution(str(mechanism),"nDodecane_"+phase) for phase in ("RK","IG")}
    prepared_states = capture_phase_states(gases)
    native_provenance = verify_native_mechanisms(args.native_mechanism_directory,meta,gases)
    record = json.loads(args.cantera_build_record.read_text()) if args.cantera_build_record else None
    libraries_before = mapped_cantera_hashes(record)
    input_paths = [Path(__file__),Path(__file__).with_name("benchmark_environment.py"),args.julia_result,mechanism]
    input_paths += [p for p in (args.source_example,args.cantera_build_record) if p is not None]
    if args.native_mechanism_directory is not None:
        input_paths += [args.native_mechanism_directory/name for name in NATIVE_MECHANISMS]
    input_hashes = {str(p.resolve()):sha(p) for p in input_paths}
    hardware = host_metadata()
    thread_checks = {"before":verify_threads()}
    hardware["load_average_start"] = os.getloadavg()
    restore_phase_states(gases,prepared_states)
    start = time.perf_counter()
    first = calculate(gases,smoke)
    first_seconds = time.perf_counter()-start
    print(f"Cantera first calculation: {len(first)} trajectories, {first_seconds:.6g} s",flush=True)
    samples,warm_delays,warm_matches_first = [],[],[]
    if not args.validate_only:
        gc.collect()
        for repetition in range(args.repetitions):
            restore_phase_states(gases,prepared_states)
            start = time.perf_counter()
            result = calculate(gases,smoke)
            samples.append(time.perf_counter()-start)
            warm_delays.append([result[key]["delay"] for key in order])
            warm_matches_first.append(require_replay(result,first,f"Cantera warm repetition {repetition+1}"))
            print(f"Cantera warm repetition {repetition+1}: {samples[-1]:.6g} s",flush=True)
    thread_checks["after_source_default"] = verify_threads(enforce=False)
    # Keep source-default timings above as the qualification baseline. This
    # second source-style loop quantifies the independent accuracy reference's
    # cost using the same output sampling rather than the validation grid.
    restore_phase_states(gases,prepared_states)
    start = time.perf_counter()
    refined_first = calculate(gases,smoke,refined=True)
    refined_first_seconds = time.perf_counter()-start
    refined_samples,refined_warm_matches_first = [],[]
    if not args.validate_only:
        gc.collect()
        for repetition in range(args.repetitions):
            restore_phase_states(gases,prepared_states)
            start = time.perf_counter()
            result = calculate(gases,smoke,refined=True)
            refined_samples.append(time.perf_counter()-start)
            refined_warm_matches_first.append(require_replay(result,refined_first,f"refined Cantera warm repetition {repetition+1}"))
            print(f"Refined Cantera warm repetition {repetition+1}: {refined_samples[-1]:.6g} s",flush=True)
    validation_start = time.perf_counter()
    cases,reference_arrays = {},{}
    for phase,temperature in case_order(smoke):
        key = f"{phase}_{temperature}"
        cases[key],reference_arrays[key+"_refined_state_on_native_grid"] = validate_case(
            gases[phase],temperature,key,native,first[key])
        reference_arrays[key+"_native_time"] = np.r_[native[key+"_time"],native[key+"_final_time"]]
        reference_arrays[key+"_published_time"] = first[key]["time"]
        reference_arrays[key+"_published_state"] = first[key]["state"]
        print(key,"correctness",cases[key]["correctness_pass"],"max dT",cases[key]["temperature_error_K"],flush=True)
    validation_seconds = time.perf_counter()-validation_start
    thread_checks["after"] = verify_threads(enforce=False)
    reference_path = args.output.with_suffix(".cantera.npz")
    np.savez(reference_path,**reference_arrays)
    libraries_after = mapped_cantera_hashes(record)
    if libraries_after!=libraries_before:
        raise ValueError("actually mapped Cantera library inventory or bytes changed during calculations")
    require_unchanged(input_hashes)
    if current_source_hashes(args.native_project)!=source_hashes:
        raise ValueError("native source inventory or bytes changed during Cantera calculations")
    libraries = {Path(p).name:value for p,value in libraries_after.items() if "libcantera" in Path(p).name}
    source_sha = record["source"]["commit"] if record else getattr(ct,"__git_commit__","unknown")
    recorded_libraries = {Path(k).name:v for k,v in (record or {}).get("library_hashes",{}).items()}
    extension_sha = libraries_after[str(Path(compiled.__file__).resolve())]
    build_matches = bool(libraries and recorded_libraries
        and all(recorded_libraries.get(k)==v for k,v in libraries.items())
        and recorded_libraries.get(Path(compiled.__file__).name)==extension_sha)
    ratio = statistics.median(samples)/statistics.median(native_samples) if samples and native_samples else None
    native_host = {k:meta.get(k,"") for k in ("cpu","system","kernel_release")}
    correct = all(case["correctness_pass"] for case in cases.values())
    all_repetitions_checked = (len(warm_matches_first)==len(samples)>=7
        and len(native_warm_matches_first)==len(native_samples)>=7
        and all(warm_matches_first) and all(native_warm_matches_first))
    native_thread_checks = tomllib.loads(meta.get("thread_checks_toml",""))
    native_thread_environment = tomllib.loads(meta.get("thread_environment_toml",""))
    native_threads_pass = (bool(native_thread_checks)
        and all(native_thread_environment.get(name)=="1" for name in THREAD_VARIABLES)
        and all(check.get("julia_threads")==check.get("blas_threads")==1
            and check.get("mkl_threads",1)==1
            and (args.benchmark_target!="apple-m4" or check.get("accelerate_threading_mode")==1)
            for check in native_thread_checks.values()))
    controlled = bool(not smoke and not args.validate_only and args.qualification=="controlled"
        and meta["qualification"]=="controlled" and len(samples)>=7 and len(native_samples)>=7
        and matches_target(hardware,args.benchmark_target) and matches_target(native_host,args.benchmark_target)
        and all(hardware[key]==native_host[key] for key in ("cpu","kernel_release"))
        and ct.__version__.startswith("4.0") and source_sha==SOURCE_COMMIT and build_matches
        and all_repetitions_checked and native_threads_pass and native_provenance["verified"])
    report = {"example":"reactors/non_ideal_shock_tube.py","scope":meta["scope"],
        "case_order":order,"temperature_sweep_K":TEMPERATURES,"pressure_Pa":40*ct.one_atm,
        "composition_moles":{"c12h26":1.,"o2":18.5,"n2":69.56},
        "end_time_s":.005,"save_stride":20,"endpoint":"first accepted step at or beyond end time",
        "mechanism_species":gases["RK"].n_species,"mechanism_reactions":gases["RK"].n_reactions,
        "cantera_mechanism_file":mechanism.name,"cantera_mechanism_sha256":sha(mechanism),
        "native_mechanism_provenance":native_provenance,
        "cantera_version":ct.__version__,"cantera_source_sha":source_sha,
        "cantera_build_record_sha256":sha(args.cantera_build_record) if record else None,
        "cantera_shared_libraries_sha256":libraries,"cantera_extension_sha256":extension_sha,
        "mapped_cantera_hashes_before":libraries_before,"mapped_cantera_hashes_after":libraries_after,
        "input_hashes_before_and_after":input_hashes,"input_source_bytes_unchanged":True,
        "cantera_build_matches_loaded_libraries":build_matches,
        "published_source_sha256":sha(args.source_example) if args.source_example else None,
        "harness_sha256":sha(__file__),"native_artifact_sha256":sha(args.julia_result),
        "environment_helper_sha256":sha(Path(__file__).with_name("benchmark_environment.py")),
        "cantera_trajectory_artifact":reference_path.name,"cantera_trajectory_artifact_sha256":sha(reference_path),
        "native_metadata":meta,"native_source_hashes":tomllib.loads(meta["source_hashes_toml"]),
        "hardware":dict(hardware,load_average_end=os.getloadavg()),
        "thread_environment":{key:os.environ.get(key) for key in THREAD_VARIABLES},
        "thread_checks":thread_checks,"native_thread_checks":native_thread_checks,
        "native_thread_environment":native_thread_environment,"native_threads_pass":native_threads_pass,
        "timestamp_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
        "benchmark_target":args.benchmark_target,"timed_scope":"preloaded mechanisms; fresh reactors, source stepping, output-state storage and OH peak; no imports, model files, plots or reference validation",
        "cantera_import_seconds":IMPORT_SECONDS,"julia_import_seconds":float(native["import_seconds"][0]),
        "cantera_first_seconds":first_seconds,"julia_first_seconds":float(native["first_seconds"][0]),
        "refined_cantera_first_seconds":refined_first_seconds,
        "refined_cantera_warm_seconds":refined_samples,
        "refined_cantera_warm_matches_first":refined_warm_matches_first,
        "refined_cantera_source_style_tolerances":{"rtol":1e-11,"atol":1e-22},
        "refined_cantera_solver_stats":{key:refined_first[key]["solver_stats"] for key in order},
        "refined_cantera_speed_ratio":statistics.median(refined_samples)/statistics.median(native_samples) if refined_samples and native_samples else None,
        "qualification_timing_baseline":"published Cantera default tolerances",
        "whole_call_reset":"restore pristine phase states outside every source-default and refined calculation timer; published per-case TP-then-composition operations are unchanged",
        "first_call_note":"one full first invocation per process; Julia includes JIT; warm repetitions reuse code and preloaded static models",
        "cantera_warm_seconds":samples,"julia_warm_seconds":native_samples,
        "cantera_warm_delays":warm_delays,"julia_warm_delays":native["warm_delays"].T.tolist(),
        "cantera_warm_matches_first":warm_matches_first,"julia_warm_matches_first":native_warm_matches_first,
        "all_repetitions_checked":all_repetitions_checked,
        "repetition_check":"complete saved trajectories and endpoints exactly match the independently checked first invocation",
        "cantera_tolerances":{"rtol":first[order[0]]["rtol"],"atol":first[order[0]]["atol"]},
        "native_solver":meta.get("solver","OrdinaryDiffEqSDIRK.KenCarp4"),
        "native_tolerances":{"rtol":float(native["native_rtol"][0]),"atol":float(native["native_atol"][0])},
        "native_temperature_atol_K":float(native["native_temperature_atol_K"][0]) if "native_temperature_atol_K" in native else None,
        "native_trial_policy":meta.get("trial_policy","clipped concentrations"),
        "native_absolute_tolerance_formula":meta.get("absolute_tolerance_formula","scalar tolerance for every state component"),
        "accuracy_note":"Native pointwise checks use independently refined Cantera tolerances; published default-tolerance results are retained, including known low-temperature convergence error.",
        "correctness_pass":correct,"reference_validation_seconds":validation_seconds,"cases":cases,
        "speed_ratio":ratio,"minimum_speed_ratio":.95,"speed_ratio_definition":"median Cantera / median Julia",
        "qualification":"controlled" if controlled else "not_qualified",
        "performance_pass":bool(controlled and correct and ratio>=.95)}
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
    print("Full-example correctness",correct,"speed ratio",ratio,"qualified",controlled,flush=True)
    if not correct:
        raise SystemExit("native shock-tube trajectories failed refined-reference checks")
    if args.qualification=="controlled" and not controlled:
        raise SystemExit("requested controlled timing failed its provenance or environment guards")


if __name__=="__main__":
    main()
