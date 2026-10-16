"""Pair the public complete QNDF engine with the unchanged pinned CT source.

smoke = first + one warm on each side; controlled = first + nine warm.
Run the Julia driver first. Source construction/solve/required output integrals
are timed; replay, provenance, artifact I/O and independent refined references
are outside timers. No prior trajectory is ever a solver initial condition.
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import os
import statistics
import sys
import time
import tomllib
from types import SimpleNamespace

_imports_started=time.perf_counter()
import numpy as np
import cantera as ct
import cantera._cantera as compiled
from ic_engine_source import (SOURCE_COMMIT,SOURCE_SHA256,source_programs,source_calculation,
    output_replay,source_checks,verify_mechanisms,sha)
import ic_engine_qndf_case as reference
from benchmark_environment import host_metadata,matches_target,loaded_library_paths,verify_numerical_threads
IMPORT_SECONDS=time.perf_counter()-_imports_started

SOURCE_FILES=("Project.toml","example/reactors/ic_engine.jl","example/reactors/ic_engine_setup.jl",
    "example/reactors/ic_engine_qndf_solver.jl","validation/ic_engine_timing.jl",
    "validation/ic_engine_timing.py","validation/ic_engine_source.py","validation/ic_engine_qndf_case.jl",
    "validation/ic_engine_qndf_case.py","validation/numerical_threads.jl","validation/benchmark_environment.py")

def source_inventory(root):
    files={p.relative_to(root).as_posix():sha(p) for folder in ("src","example/reactors/engine_qndf")
           for p in (root/folder).rglob("*") if p.is_file()}
    for name in SOURCE_FILES:files[name]=sha(root/name)
    return files

def require_inventory(actual,expected,label):
    if actual!=expected:
        missing=sorted(set(expected)-set(actual));added=sorted(set(actual)-set(expected))
        changed=sorted(k for k in actual.keys()&expected.keys() if actual[k]!=expected[k])
        raise ValueError(f"{label} membership/bytes changed: missing={missing}, added={added}, changed={changed}")
    return True

def package_sources(module):
    directory=Path(module.__file__).resolve().parent
    return {p.relative_to(directory).as_posix():sha(p) for p in directory.rglob("*")
            if p.is_file() and "__pycache__" not in p.parts and p.suffix not in (".pyc",".pyo")}

def prepare_reference_imports():
    # NumPy unique() checks np.ma.is_masked in the independent comparator.
    # Load these known modules before freezing dependency membership/timers.
    import numpy.ma
    import numpy.ma.core
    import numpy.ma.extras

def dependency_snapshot():
    # Pin installed NumPy/Cantera package membership, interpreter, and loaded
    # source modules. Package inventories also cover their extension modules.
    modules={}
    for name,module in list(sys.modules.items()):
        path=getattr(module,"__file__",None)
        if path and Path(path).is_file():modules[name]={"path":str(Path(path).resolve()),"sha256":sha(path)}
    return dict(python_executable_sha256=sha(sys.executable),python_version=sys.version,
                numpy_version=np.__version__,cantera_version=ct.__version__,
                numpy_files=package_sources(np),cantera_files=package_sources(ct),loaded_module_files=modules)

def mapped_snapshot():
    paths=loaded_library_paths()
    return {str(p):sha(p) for p in paths if p.is_file()}

def verify_cantera_build(mapped,record):
    if record.get("source",{}).get("commit")!=SOURCE_COMMIT or not ct.__version__.startswith("4.0."):
        raise ValueError("pinned Cantera 4 source/build required")
    extension=Path(compiled.__file__).resolve()
    selected={p:h for p,h in mapped.items() if "libcantera" in Path(p).name or Path(p)==extension}
    if str(extension) not in selected or not any("libcantera" in Path(p).name for p in selected):
        raise ValueError("actual mapped Cantera extension/shared library absent")
    expected={Path(p).name:h for p,h in record.get("library_hashes",{}).items()}
    if not all(expected.get(Path(p).name)==h for p,h in selected.items()):
        raise ValueError("mapped Cantera bytes differ from the build record")
    return selected

def full_mechanism_match(source,prepared):
    original=ct.Solution(str(source),"nDodecane_IG");native=ct.Solution(str(prepared))
    def equal(a,b,path):
        if isinstance(a,dict) and isinstance(b,dict):
            if set(a)!=set(b):raise ValueError(f"mechanism keys differ at {path}")
            for key in a:equal(a[key],b[key],f"{path}.{key}")
        elif isinstance(a,(list,tuple,np.ndarray)) and isinstance(b,(list,tuple,np.ndarray)):
            if len(a)!=len(b):raise ValueError(f"mechanism dimensions differ at {path}")
            for i,(x,y) in enumerate(zip(a,b)):equal(x,y,f"{path}[{i}]")
        elif isinstance(a,(float,int,np.number)) and isinstance(b,(float,int,np.number)):
            # Match the existing preparation tolerance, across every numerical
            # species/reaction datum, allowing YAML's decimal serialization.
            if not np.isclose(a,b,rtol=2e-12,atol=0.,equal_nan=False):
                raise ValueError(f"mechanism number differs at {path}")
        elif a!=b:raise ValueError(f"mechanism metadata differs at {path}")
    if original.species_names!=native.species_names or original.n_reactions!=native.n_reactions:
        raise ValueError("source/prepared mechanism membership differs")
    for i,(a,b) in enumerate(zip(original.species(),native.species())):equal(dict(a.input_data),dict(b.input_data),f"species[{i}]")
    for i,(a,b) in enumerate(zip(original.reactions(),native.reactions())):equal(dict(a.input_data),dict(b.input_data),f"reaction[{i}]")
    return True

def native_dependency_files(directory):
    directory=Path(directory)
    files={p.relative_to(directory).as_posix():sha(p) for folder in ("src","ext","deps","lib")
           for p in (directory/folder).rglob("*") if p.is_file()}
    for name in ("Project.toml","Artifacts.toml"):
        if (directory/name).is_file():files[name]=sha(directory/name)
    return files

def require_native_files(record):
    inputs=record["input_before"]
    for name,path in inputs["file_paths"].items():
        if sha(path)!=inputs[name+"_sha256"]:raise ValueError(f"current native {name} bytes changed")
    for identity,package in inputs["dependencies"].items():
        require_inventory(native_dependency_files(package["path"]),package["files"],f"native dependency {identity}")
    pins=record["library_pins"]
    require_inventory({path:sha(path) for path in pins},pins,"native pre-pinned libraries")
    prior=record["runtime_before"]
    for index,run in enumerate(record["runs"]):
        before,after=run["runtime_before"],run["runtime_after"]
        require_inventory(before,prior,"native between-call runtime")
        if before["environment"]!=after["environment"] or not all(v==1 for v in after["settings"].values()):
            raise ValueError("native actual thread settings changed")
        initial,final=set(before["mapped"]),set(after["mapped"])
        if not (initial<=final if index==0 else initial==final):raise ValueError("native mapped membership changed")
        for snapshot in (before,after):
            if set(snapshot["mapped_sha256"])!=set(snapshot["mapped"]):raise ValueError("native mapped hash membership differs")
            if any(pins.get(path)!=value for path,value in snapshot["mapped_sha256"].items()):
                raise ValueError("native mapped bytes differ from pre-run pins")
        prior=after
    require_inventory(record["runtime_after"],prior,"native final runtime")
    return True

def require_native_timing(record,mode,repeats,root,native):
    if not record.get("complete") or record.get("mode")!=mode or record.get("warm_repetitions")!=repeats:
        raise ValueError("native first/warm sequence is incomplete or uses a different mode")
    if record.get("public_callable")!="NativeEngineQNDF.solve_ic_engine_qndf":
        raise ValueError("validated public QNDF callable required")
    if record.get("wall_clock_exclusions")!=["seconds_including_first_specialization"]:
        raise ValueError("unexpected exact-replay exclusion")
    if len(record.get("runs",[]))!=repeats+1 or len(record.get("warm_seconds",[]))!=repeats:
        raise ValueError("missing native first/warm calls")
    if [r["label"] for r in record["runs"]]!=["first"]+[f"warm-{i}" for i in range(1,repeats+1)]:
        raise ValueError("native repetition order differs")
    durations=[record["first_seconds"]]+record["warm_seconds"]
    if not all(np.isfinite(t) and t>0 for t in durations) or durations!=[r["seconds"] for r in record["runs"]]:
        raise ValueError("native durations are invalid or differ from individual calls")
    current=source_inventory(root)
    require_inventory(record["input_before"],record["input_after"],"native before/after inputs")
    require_inventory(record["input_before"]["source_inventory"],current,"native/current source")
    for run in record["runs"]:
        if not run.get("checks_pass") or not run.get("replay_pass"):
            raise ValueError("native physical/replay check failed")
        require_inventory(run["inputs"],record["input_before"],"native repetition inputs")
        if sha(run["archive"])!=run["archive_sha256"] or sha(run["archive"]+".replay.jls")!=run["payload_sha256"]:
            raise ValueError("native saved repetition bytes changed")
        if sha(Path(run["archive"]).with_suffix(".toml"))!=run["metadata_sha256"]:
            raise ValueError("native saved repetition metadata changed")
    if Path(record["runs"][0]["archive"]).resolve()!=native.resolve():raise ValueError("wrong native first artifact")
    require_native_files(record)
    return True

def save_source(result,directory,label):
    directory.mkdir(parents=True,exist_ok=True)
    path=directory/(label+".npz")
    if path.exists():raise ValueError("preserve previous source artifact")
    np.savez(path,**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
    metadata={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
    metadata.update(archive_sha256=sha(path))
    path.with_suffix(".json").write_text(json.dumps(metadata,indent=2,allow_nan=False)+"\n")
    return {"array_path":str(path.resolve()),"array_sha256":sha(path),"metadata_sha256":sha(path.with_suffix(".json"))}

def exact_source_equal(actual,expected):
    if type(actual) is not type(expected):return False
    if isinstance(actual,np.ndarray):
        return actual.shape==expected.shape and actual.dtype==expected.dtype and actual.tobytes()==expected.tobytes()
    if isinstance(actual,dict):
        return set(actual)==set(expected) and all(exact_source_equal(actual[k],expected[k]) for k in actual)
    if isinstance(actual,(list,tuple)):
        return len(actual)==len(expected) and all(exact_source_equal(a,b) for a,b in zip(actual,expected))
    return actual==expected

def require_source_result(actual,expected=None):
    if not source_checks(actual):raise ValueError("source physical/output checks failed")
    if expected is not None and not (output_replay(actual,expected) and exact_source_equal(actual,expected)):
        raise ValueError("source warm output/counter replay failed")
    return True

def require_saved_source(saved):
    path=Path(saved["array_path"])
    if sha(path)!=saved["array_sha256"] or sha(path.with_suffix(".json"))!=saved["metadata_sha256"]:
        raise ValueError("saved source output bytes changed")
    return True

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("native","native-mechanism","source-example","source-mechanism","cantera-build-record","output"):
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--mode",choices=("smoke","controlled"),default="smoke")
    parser.add_argument("--target",choices=("wsl","apple-m4"),required=True)
    args=parser.parse_args();repeats=1 if args.mode=="smoke" else 9
    if args.output.exists():raise ValueError("preserve earlier paired receipt")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    artifact_dir=args.output.with_suffix(".artifacts")
    if artifact_dir.exists():raise ValueError("preserve earlier source artifacts")
    artifact_dir.mkdir()
    native_record_path=Path(str(args.native)+".timing.toml")
    record=tomllib.loads(native_record_path.read_text())
    root=Path(__file__).resolve().parents[1]
    require_native_timing(record,args.mode,repeats,root,args.native)
    if sha(args.source_example)!=SOURCE_SHA256 or sha(args.source_mechanism)!=reference.SOURCE_MECHANISM_SHA256:
        raise ValueError("pinned source example and mechanism required")
    native=np.load(args.native,allow_pickle=False);native_meta=tomllib.loads(args.native.with_suffix(".toml").read_text())
    reference.require_native(native,native_meta)
    verify_mechanisms(args.source_mechanism,args.native_mechanism,native_meta)
    full_mechanism_match(args.source_mechanism,args.native_mechanism)
    programs=source_programs(args.source_example,args.source_mechanism)
    hardware=host_metadata()
    if not matches_target(hardware,args.target):raise ValueError("wrong physical benchmark target")
    if hardware["cpu"]!=record["cpu"] or hardware["kernel_release"]!=record["kernel_release"]:
        raise ValueError("native and CT must run on the same host/kernel")
    thread_keys=("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")+(("VECLIB_MAXIMUM_THREADS",) if hardware["system"]=="Darwin" else ())
    if any(os.environ.get(k)!="1" for k in thread_keys):raise ValueError("one thread must be requested on every backend")
    # Exercise serialization/import paths before freezing loaded-module and
    # library membership, without constructing or stepping an engine.
    prepare_reference_imports()
    np.savez(artifact_dir/"serialization-preflight.npz",empty=np.zeros(0))
    thread_before=verify_numerical_threads(set_accelerate=True)
    build=json.loads(args.cantera_build_record.read_text())
    dependencies=dependency_snapshot();libraries=mapped_snapshot();verify_cantera_build(libraries,build)
    named={"native":args.native,"native_record":native_record_path,"native_metadata":args.native.with_suffix(".toml"),
           "native_mechanism":args.native_mechanism,"native_sidecar":Path(str(args.native_mechanism)+".npz"),
           "source_example":args.source_example,"source_mechanism":args.source_mechanism,"build_record":args.cantera_build_record}
    inputs={name:sha(path) for name,path in named.items()};sources=source_inventory(root)
    report=dict(complete=False,performance_qualified=False,mode=args.mode,target=args.target,hardware=hardware,
        source_commit=SOURCE_COMMIT,source_example_sha256=SOURCE_SHA256,source_mechanism_sha256=reference.SOURCE_MECHANISM_SHA256,
        inputs_before=inputs,source_before=sources,dependencies_before=dependencies,libraries_before=libraries,
        thread_before=thread_before,cantera_import_seconds=IMPORT_SECONDS,runs=[],cantera_warm_seconds=[],
        julia_first_seconds=record["first_seconds"],julia_warm_seconds=record["warm_seconds"],
        timed_scope="fresh source construction, eight-revolution solve, original output properties and sampled integrals",
        native_scope=record["scope"],accuracy_reference="independent full accepted-history 25-segment endpoint/quadrature oracle; excluded from timing ratio")
    def checkpoint():args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    def guard():
        require_inventory({name:sha(path) for name,path in named.items()},inputs,"input files")
        require_inventory(source_inventory(root),sources,"source inventory")
        require_inventory(dependency_snapshot(),dependencies,"Python dependency inventory")
        observed=mapped_snapshot();require_inventory(observed,libraries,"mapped-library inventory")
        verify_cantera_build(observed,build)
        return verify_numerical_threads()
    checkpoint();first=None
    try:
        gc.collect()
        for index in range(repeats+1):
            label="first" if index==0 else f"warm-{index}"
            before_threads=guard()
            started=time.perf_counter();result=source_calculation(programs);seconds=time.perf_counter()-started
            saved=save_source(result,artifact_dir,label)
            run=dict(label=label,seconds=seconds,saved=saved,checks_pass=False,replay_pass=False,threads_before=before_threads)
            report["runs"].append(run);checkpoint()
            run["threads_after"]=guard()
            require_source_result(result,first)
            run["checks_pass"]=run["replay_pass"]=True
            if index==0:first=result;report["cantera_first_seconds"]=seconds
            else:report["cantera_warm_seconds"].append(seconds)
            checkpoint();print(f"CT engine {label}: {seconds} s; full source replay passes",flush=True)
        # Reuse the exact corrected reference implementation, without import
        # cycles or copied integrator/reference equations.
        comparison={"checks_pass":False}
        comparison_path=artifact_dir/"independent.json"
        started=time.perf_counter()
        try:
            reference.compare(SimpleNamespace(native=args.native,native_mechanism=args.native_mechanism,
                source_example=args.source_example,source_mechanism=args.source_mechanism,output=comparison_path),comparison)
        finally:
            report["reference_seconds"]=time.perf_counter()-started
            comparison_path.write_text(json.dumps(comparison,indent=2,allow_nan=False)+"\n")
        if not comparison["checks_pass"]:raise ValueError("independent accepted-history/integral gates failed")
        report["thread_after"]=guard()
        require_native_timing(tomllib.loads(native_record_path.read_text()),args.mode,repeats,root,args.native)
        for run in report["runs"]:
            require_saved_source(run["saved"])
            if not np.isfinite(run["seconds"]) or run["seconds"]<=0:raise ValueError("invalid source duration")
        report["inputs_after"]={name:sha(path) for name,path in named.items()}
        report["source_after"]=source_inventory(root);report["libraries_after"]=mapped_snapshot()
        report["dependencies_after"]=dependency_snapshot();report["accuracy"]=comparison
        ratio=statistics.median(report["cantera_warm_seconds"])/statistics.median(record["warm_seconds"])
        report.update(complete=True,correctness_pass=True,speed_ratio=ratio,minimum_speed_ratio=.95,
                      performance_qualified=args.mode=="controlled" and ratio>=.95,
                      speed_ratio_definition="median source-default CT / median public QNDF complete-call seconds")
        checkpoint()
    except BaseException as error:
        report["error"]=str(error);checkpoint();raise
    print("Engine paired checks complete; ratio",report["speed_ratio"],"qualified",report["performance_qualified"])

if __name__=="__main__":main()
