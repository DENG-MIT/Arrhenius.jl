"""Time complete source-default flame sequences, with separate physical references."""
import argparse,ast,gc,hashlib,json,os,platform,re,statistics,subprocess,sys,time,tomllib
from pathlib import Path
THREAD_KEYS=["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","JULIA_NUM_THREADS"]
for key in THREAD_KEYS:os.environ[key]="1"
import numpy as np
import cantera as ct
from benchmark_environment import host_metadata,matches_target,loaded_library_paths,verify_numerical_threads

PINNED_COMMIT="726522be4e2a13454d8415b7ef799d621f665cf3"
SOURCE_EXAMPLES={
    "free":("adiabatic_flame","f2a845e0e2b0c06d466eadeba9b1be1d9b8124aadf94e83230527d4124b5840d"),
    "burner":("burner_flame","8c95f440748a2559b4f2a940b6b1864144e0668704a8d894e5923663278632c2"),
    "fixed":("flame_fixed_T","102af72f0349116ce7a1258bb1ce0b06cc04ebe7359358b5c4d635799bbd01be"),
    "ion-free":("ion_free_flame","70d0b6d647ce98b7ca7eaa9b749be7d67fe8dccbad2b560b4a2225499ec47878"),
    "ion-burner":("ion_burner_flame","69f4d02bc8482f51c439a879fda59a947ba652a3735cdc4315b3f81d7f47481b")}
MECHANISM_HASHES={"h2o2":"0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4",
    "gri30":"06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345",
    "gri30_ion":"1409cc65488ea885118257c65d0a0b72fbcf5968f7b4a0572f637f9bdec9a66a"}
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

if platform.system()=="Linux":
    if not hasattr(time,"CLOCK_MONOTONIC_RAW") or not hasattr(time,"clock_gettime_ns"):
        raise RuntimeError("CLOCK_MONOTONIC_RAW nanosecond clock is unavailable on this Linux host")
    def benchmark_elapsed_ns():return time.clock_gettime_ns(time.CLOCK_MONOTONIC_RAW)
    ELAPSED_CLOCK="CLOCK_MONOTONIC_RAW"
else:
    def benchmark_elapsed_ns():return time.perf_counter_ns()
    ELAPSED_CLOCK="monotonic"

def checked_elapsed_clock(data):
    try:
        value=bytes(np.asarray(data["elapsed_clock_utf8"],dtype=np.uint8).ravel()).decode("utf-8")
    except (KeyError,ValueError,UnicodeError) as exc:
        raise RuntimeError("missing or invalid native elapsed clock") from exc
    if value!=ELAPSED_CLOCK:
        raise RuntimeError("Julia and Cantera elapsed clocks differ")
    return value


def current_source_hashes(project):
    return {str(path.relative_to(project)):digest(path) for path in sorted((project/"src").rglob("*.jl"))}

def require_unchanged(expected):
    changed=[path for path,value in expected.items() if not Path(path).is_file() or digest(path)!=value]
    if changed:raise RuntimeError(f"benchmark input/source bytes changed: {changed}")

def checked_times(values,repetitions,stages):
    values=np.asarray(values,dtype=float)
    if values.shape!=(repetitions+1,stages) or not np.all(np.isfinite(values)) or not np.all(values>0):
        raise RuntimeError("missing or invalid first/warm stage timings")
    return np.sum(values,axis=1)

def checked_clock_observations(values,repetitions):
    values=np.asarray(values,dtype=float)
    if (values.shape!=(repetitions+1,4) or not np.all(np.isfinite(values)) or
            np.any(values[:,0]<=0) or np.any(values[:,1]<0) or
            np.any(values[:,2:]<0) or np.any(values[:,2:]!=np.floor(values[:,2:]))):
        raise RuntimeError("invalid process-clock observations")
    return values.tolist()


def make_clock_observer():
    # Diagnostic only: CPU IDs are WSL/Linux logical IDs, not host core types.
    import ctypes
    libc=ctypes.CDLL(None)
    getcpu=libc.sched_getcpu;getcpu.argtypes=[];getcpu.restype=ctypes.c_int
    def observe():
        cpu=getcpu()
        if cpu<0:raise RuntimeError("sched_getcpu failed")
        return benchmark_elapsed_ns()/1e9,time.process_time(),cpu
    return observe


def verified_cantera_libraries(build):
    import cantera._cantera as compiled
    paths={path.resolve() for path in loaded_library_paths() if "cantera" in path.name.lower() and path.is_file()}
    paths.add(Path(compiled.__file__).resolve())
    hashes={str(path):digest(path) for path in sorted(paths)}
    if len(hashes)<2 or any(value not in build.get("library_hashes",{}).values() for value in hashes.values()):
        raise RuntimeError(f"loaded Cantera library/build-record mismatch: {hashes}")
    return hashes

def verify_source_inputs(case,source,mechanism,provenance,profile):
    name,expected=SOURCE_EXAMPLES[case]
    if digest(source)!=expected or provenance.get("source_commit")!=PINNED_COMMIT or provenance.get("examples",{}).get(name,{}).get("sha256")!=expected:
        raise ValueError("example source/provenance does not match pinned Cantera")
    if digest(mechanism)!=MECHANISM_HASHES[mechanism.stem]:raise ValueError("mechanism bytes differ from pinned stock source")
    if case=="fixed":
        values={}
        for node in ast.parse(source.read_text()).body:
            if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ("zloc","tvalues"):
                values[node.targets[0].id]=np.asarray(ast.literal_eval(node.value.args[0]),dtype=float)
        if not np.array_equal(values.get("zloc"),profile["positions"]) or not np.array_equal(values.get("tvalues"),profile["temperatures"]):
            raise ValueError("prescribed temperature profile differs from exact source values")

def loaded_libraries():
    return [str(path) for path in loaded_library_paths()]


def verify_threads(*,set_accelerate=False):
    result=verify_numerical_threads(set_accelerate=set_accelerate)
    result["environment"]={key:os.environ[key] for key in THREAD_KEYS}
    return result

def cantera_sequence(gas,case,profile,save_profiles=False):
    if case.startswith("ion-"):
        from ion_source_sequence import run_ion_source_sequence
        source_case=case.removeprefix("ion-")
        if gas.transport_model!="ionized-gas":raise RuntimeError("ion mechanism transport model mismatch")
        seconds,points,snapshots=run_ion_source_sequence(gas,source_case,clock_ns=benchmark_elapsed_ns,loglevel=0)
        return seconds,points,snapshots if save_profiles else {}
    free=case=="free";fixed=case=="fixed"
    modes=["mass","mass-soret","multi","multi-soret"] if free else ["mole","multi"]
    seconds=[];points=[];snapshots={}
    # Restore the preloaded gas's initial transport manager, corresponding to
    # fresh mechanism loading outside the timed source calculation.
    gas.transport_model="mixture-averaged"
    for stage,mode in enumerate(modes):
        multi=mode.startswith("multi");start=benchmark_elapsed_ns()
        if stage==0:
            gas.TPX=(300.,ct.one_atm,"H2:1.1,O2:1,AR:5") if free else (373.7,ct.one_atm,"CH4:.65,O2:1,N2:3.76") if fixed else (373.,.05*ct.one_atm,"H2:1.5,O2:1,AR:7")
            flame=ct.FreeFlame(gas,width=.03) if free else ct.BurnerFlame(gas,width=.01 if fixed else .5)
            if not free:flame.burner.mdot=.04 if fixed else .06
            if fixed:
                flame.flame.set_fixed_temp_profile(profile["positions"]/.01,profile["temperatures"])
                flame.energy_enabled=False
        flame.transport_model="multicomponent" if multi else "mixture-averaged"
        flame.flux_gradient_basis="mass" if free else "molar"
        flame.soret_enabled=mode.endswith("soret")
        slope=.06 if free else (.1 if multi else .3) if fixed else .05
        curve=.12 if free else (.2 if multi else 1.) if fixed else .1
        flame.set_refine_criteria(ratio=3.,slope=slope,curve=curve)
        # Preserve the published default tolerances and Jacobian mode here.
        # Refined accuracy-reference solves are entirely outside this runner.
        flame.solve(loglevel=0,auto=stage==0 and not fixed)
        if save_profiles:
            snapshots[mode]=dict(grid=flame.grid.copy(),T=flame.T.copy(),Y=flame.Y.copy(),velocity=flame.velocity.copy(),
                inlet_Y=flame.inlet.Y.copy() if free else flame.burner.Y.copy(),P=[flame.P])
        seconds.append((benchmark_elapsed_ns()-start)/1e9);points.append(len(flame.grid))
    return seconds,points,snapshots

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("parameters",type=Path);p.add_argument("references",type=Path);p.add_argument("output",type=Path)
    p.add_argument("case",choices=["free","burner","fixed","ion-free","ion-burner"])
    p.add_argument("--project",required=True,type=Path);p.add_argument("--julia",required=True)
    p.add_argument("--build-record",required=True,type=Path)
    p.add_argument("--native-record",required=True,type=Path)
    p.add_argument("--source-examples",required=True,type=Path,help="pinned Cantera samples/python/onedim directory")
    p.add_argument("--target",required=True,choices=["wsl","apple-m4"])
    p.add_argument("--reps",type=int,default=9);p.add_argument("--formal",action="store_true")
    p.add_argument("--order",choices=["cantera-first","julia-first"],default="cantera-first")
    p.add_argument("--clock-diagnostics",action="store_true",help="Linux process clocks and observed CPU IDs; incompatible with --formal")
    a=p.parse_args()
    if a.clock_diagnostics and (a.formal or a.target!="wsl" or platform.system()!="Linux"):
        p.error("clock diagnostics require Linux/WSL and cannot qualify formal performance")
    if a.formal and a.target=="wsl" and ELAPSED_CLOCK!="CLOCK_MONOTONIC_RAW":
        p.error("formal WSL qualification requires the CLOCK_MONOTONIC_RAW elapsed clock")
    observe_clock=make_clock_observer() if a.clock_diagnostics else None
    if a.reps<5 or (a.formal and a.reps<9):p.error("at least five warm repetitions required; formal qualification requires nine")
    if not ct.__version__.startswith("4.0"):p.error("Cantera4 required")
    build=json.loads(a.build_record.read_text())
    if build.get("source",{}).get("commit")!=PINNED_COMMIT:p.error("build record does not identify pinned pristine source")
    native_record=json.loads(a.native_record.read_text())
    if not re.fullmatch(r"[0-9a-f]{40}",native_record.get("commit","")):p.error("native record requires a full checkpoint commit")
    source_hashes=current_source_hashes(a.project)
    if source_hashes!=native_record.get("source_hashes"):p.error("native source files do not match the pinned checkpoint record")
    host=host_metadata()
    if not matches_target(host,a.target):p.error(f"host mismatch: {host}")
    if a.output.exists() and any(a.output.iterdir()):p.error("output directory must be empty to exclude stale artifacts")
    a.output.mkdir(parents=True,exist_ok=True)
    native=a.output/"native";native.mkdir(exist_ok=True)
    cantera=a.output/"cantera";cantera.mkdir(exist_ok=True)
    ion=a.case.startswith("ion-"); source_case=a.case.removeprefix("ion-") if ion else a.case
    mechanism=a.parameters/("gri30_ion.yaml" if ion else "gri30.yaml" if a.case=="fixed" else "h2o2.yaml")
    gas=ct.Solution(str(mechanism))
    expected=(56,331) if ion else (53,325) if a.case=="fixed" else (10,29)
    if (gas.n_species,gas.n_reactions)!=expected or (ion and gas.transport_model!="ionized-gas"):p.error("exact stock mechanism required")
    profile=np.load(a.parameters/"fixed-profile.npz") if a.case=="fixed" else None
    provenance=json.loads((a.parameters/"provenance.json").read_text())
    example=a.source_examples/(SOURCE_EXAMPLES[a.case][0]+".py")
    verify_source_inputs(a.case,example,mechanism,provenance,profile)
    threads=verify_threads(set_accelerate=True)
    libraries_before=verified_cantera_libraries(build)
    driver_paths=[Path(__file__),Path(__file__).with_suffix(".jl"),*[Path(__file__).with_name(name) for name in
        ("conservative_source_accuracy.py","ionized_source_accuracy.py","ion_source_sequence.py","flame_benchmarks.py","benchmark_environment.py","numerical_threads.jl")]]
    example_paths=[a.project/"example/flames"/name for name in
        ("source_flame_sequence.jl","source_ion_flame_sequence.jl","adiabatic_flame.jl","burner_flame.jl","flame_fixed_T.jl")]
    input_paths=[*driver_paths,*example_paths,example,mechanism,Path(str(mechanism)+".npz"),
        a.parameters/"provenance.json",a.build_record,a.native_record,a.project/"Project.toml",a.project/"Manifest.toml",
        *sorted(a.references.glob(source_case+"-*.npz"))]
    if not ion: input_paths.append(Path(str(mechanism)+".multicomponent.npz"))
    if a.case=="fixed":input_paths.append(a.parameters/"fixed-profile.npz")
    input_hashes={str(path.resolve()):digest(path) for path in input_paths}
    report=dict(case=a.case,formal=a.formal,passed=False,performance_pass=False,date=time.strftime("%Y-%m-%d %H:%M:%S %z"),host=host,target=a.target,
        scope="Sum of construction/initialization/adaptive-solve stages and required numerical profiles in the complete published stage sequence. Mechanism/sidecar loading, validation and file output excluded. No refined-reference solve is timed.",
        compilation_scope="Julia runtime startup and using/imports are outside timers. Specialization of run_source_sequence before entry to its internal stage timers is also excluded; the first measured sequence is not whole-program cold latency. Only JIT triggered after a stage timer starts can enter its measurement. First repetition is recorded separately and excluded from warm medians.",
        warmup_gc_policy="One explicit collection after the excluded first complete sequence on each runtime; automatic GC remains enabled and timed during measured calculations.",
        clock_diagnostics=a.clock_diagnostics,elapsed_clock=ELAPSED_CLOCK,
        clock_columns=["wall_seconds","process_cpu_seconds","cpu_before","cpu_after"],
        clock_scope="Outer complete calculation call; process CPU includes all process threads. CPU IDs are observations at call boundaries and cannot rule out migrations within a call.",
        timing_order=a.order,threads=threads,cantera_version=ct.__version__,cantera_build_record=build,
        cantera_build_record_sha256=digest(a.build_record),mechanism_sha256=digest(mechanism),
        sidecar_sha256=digest(str(mechanism)+".npz"),multicomponent_sha256=None if ion else digest(str(mechanism)+".multicomponent.npz"),
        fixed_temperature_profile_sha256=digest(a.parameters/"fixed-profile.npz") if a.case=="fixed" else None,
        parameter_provenance=provenance,original_source_example=dict(path=str(example),sha256=digest(example)),
        reference_hashes={path.name:digest(path) for path in sorted(a.references.glob(source_case+"-*.npz"))},
        benchmark_driver_hashes={path.name:digest(path) for path in driver_paths},input_hashes=input_hashes,
        native_example_hashes={path.name:digest(path) for path in example_paths},
        actual_loaded_cantera_hashes_before=libraries_before,python_executable_sha256=digest(sys.executable),
        native_source_hashes=source_hashes,native_commit=native_record.get("commit"),native_record_sha256=digest(a.native_record))
    def run_ct():
        times=[];nodes=[];baseline=None;replay=[];thread_checks=[];clock_rows=[]
        for repetition in range(a.reps+1):
            clock_before=observe_clock() if observe_clock is not None else None
            elapsed,points,snapshots=cantera_sequence(gas,a.case,profile,True)
            if observe_clock is not None:
                after=observe_clock()
                clock_rows.append([after[0]-clock_before[0],after[1]-clock_before[1],clock_before[2],after[2]])
            thread_checks.append(verify_threads())
            times.append(elapsed);nodes.append(points)
            if baseline is None:baseline=snapshots
            for mode,data in snapshots.items():
                errors={};hasher=hashlib.sha256()
                for key in sorted(data):
                    value=np.asarray(data[key]);initial=np.asarray(baseline[mode][key])
                    if value.shape!=initial.shape or not np.all(np.isfinite(value)) or not np.allclose(value,initial,rtol=1e-12,atol=1e-14):
                        raise RuntimeError(f"Cantera numerical replay failed: {repetition} {mode} {key}")
                    errors[key]=float(np.max(np.abs(value-initial)))
                    hasher.update(key.encode());hasher.update(str(value.shape).encode());hasher.update(value.tobytes())
                replay.append(dict(repetition=repetition,mode=mode,sha256=hasher.hexdigest(),max_abs=errors))
                if repetition==0:np.savez(cantera/f"{a.case}-{mode}-first.npz",**data)
            print("cantera",a.case,repetition,sum(elapsed),elapsed,points,flush=True)
            if repetition==0:
                gc_start=benchmark_elapsed_ns();gc.collect()
                report["cantera_warmup_gc_seconds"]=(benchmark_elapsed_ns()-gc_start)/1e9
            if repetition==a.reps:
                for mode,data in snapshots.items():np.savez(cantera/f"{a.case}-{mode}-0.npz",**data)
        if a.clock_diagnostics:report["cantera_clock_observations"]=checked_clock_observations(clock_rows,a.reps)
        report["cantera_stage_seconds"]=times;report["cantera_stage_points"]=nodes
        report["cantera_repetition_thread_checks"]=thread_checks
        report["cantera_replay"]=dict(passed=True,checked_stages=len(replay),relative_tolerance=1e-12,absolute_tolerance=1e-14,records=replay)
    def run_julia():
        command=[a.julia,f"--project={a.project}",str(Path(__file__).with_suffix(".jl")),str(a.parameters),str(native),a.case,str(a.reps)]
        if a.clock_diagnostics:command.append("clock-diagnostics")
        result=subprocess.run(command,capture_output=True,text=True,env=dict(os.environ))
        (a.output/"julia.log").write_text(result.stdout+result.stderr)
        print(result.stdout,end="",flush=True)
        if result.returncode:raise RuntimeError(f"Julia failed: {result.stderr[-3000:]}")
        data=np.load(native/"timings.npz")
        report["julia_elapsed_clock"]=checked_elapsed_clock(data)
        if a.clock_diagnostics:
            if not bool(data["clock_diagnostics"][0]):raise RuntimeError("native clock diagnostic was not enabled")
            report["julia_clock_observations"]=checked_clock_observations(data["clock_observations"],a.reps)
        if int(data["julia_threads"][0])!=1 or int(data["blas_threads"][0])!=1:raise RuntimeError("Julia thread guard failed")
        if platform.system()=="Darwin" and int(data["accelerate_threading"][0])!=1:raise RuntimeError("Julia Accelerate thread guard failed")
        text=lambda key:bytes(data[key+"_utf8"]).decode()
        if text("kernel")!=host["kernel_release"]:raise RuntimeError("Julia and Cantera host/kernel mismatch")
        if Path(text("package_path")).resolve()!=(a.project/"src/Arrhenius.jl").resolve():raise RuntimeError("Julia loaded a different Arrhenius checkout")
        helper=a.project/"example/flames/"/("source_ion_flame_sequence.jl" if ion else "source_flame_sequence.jl")
        if Path(text("shared_calculation_path")).resolve()!=helper.resolve() or text("shared_calculation_sha256")!=digest(helper):
            raise RuntimeError("Julia source calculation does not match the public example helper")
        report["shared_public_calculation_verified"]=True
        report["julia_stage_seconds"]=data["stage_seconds"].tolist();report["julia_stage_points"]=data["stage_points"].tolist()
        gc_seconds=np.asarray(data["automatic_gc_seconds"],dtype=float)
        if gc_seconds.shape!=(a.reps+1,) or not np.all(np.isfinite(gc_seconds)) or np.any(gc_seconds<0):raise RuntimeError("invalid native automatic-GC observations")
        report["julia_automatic_gc_seconds"]=gc_seconds.tolist()
        report["julia_warmup_gc_seconds"]=float(data["warmup_gc_seconds"][0])
        report["julia_runtime"]={key:text(key) for key in ["julia_version","kernel","machine","package_path","blas_config","loaded_libraries"]}
        report["julia_runtime"]["accelerate_threading"]=int(data["accelerate_threading"][0])
        checks=tomllib.loads(text("thread_checks_toml"))
        if len(checks)!=a.reps+2 or any(value!=1 for check in checks.values() for value in check.values()):
            raise RuntimeError("incomplete or failed Julia repetition thread checks")
        actual_sources=tomllib.loads(text("source_hashes_toml"))
        if not data["source_hashes_unchanged"][0] or actual_sources!=source_hashes:
            raise RuntimeError("Julia source bytes changed or differ from recorded checkpoint")
        report["julia_repetition_thread_checks"]=checks
        report["julia_source_bytes_verified"]=True
        report["julia_runtime"]["julia_executable_sha256"]=text("julia_executable_sha256")
        report["julia_library_hashes"]=tomllib.loads(text("library_hashes_toml"))
        expected=(a.reps+1)*(4 if a.case=="free" else 2)
        if int(data["replay_checked_stages"][0])!=expected:raise RuntimeError("incomplete Julia replay checks")
        report["julia_replay"]=dict(passed=True,checked_stages=expected,relative_tolerance=1e-12,absolute_tolerance=1e-14,
            hashes=text("replay_hashes").splitlines(),fields=text("replay_fields").splitlines(),max_abs=data["replay_max_abs"].tolist())
        report["julia_blas_library_hashes"]={path:digest(path) for path in text("loaded_libraries").splitlines()
            if any(word in Path(path).name.lower() for word in ["blas","lapack"]) and Path(path).is_file()}
    try:
        for run in ([run_ct,run_julia] if a.order=="cantera-first" else [run_julia,run_ct]):run()
        report["threads_after"]=verify_threads()
        hashes=verified_cantera_libraries(build)
        if hashes!=libraries_before:raise RuntimeError("loaded Cantera libraries changed during calculation")
        report["actual_loaded_cantera_hashes"]=hashes
        checker=Path(__file__).with_name("ionized_source_accuracy.py" if ion else "conservative_source_accuracy.py")
        accuracy=a.output/"accuracy.json"
        result=subprocess.run([sys.executable,str(checker),str(mechanism),str(native),str(a.references),a.case,str(accuracy)],capture_output=True,text=True)
        print(result.stdout,end="",flush=True)
        report["accuracy"]=json.loads(accuracy.read_text()) if accuracy.is_file() else dict(passed=False,error=result.stderr)
        require_unchanged(input_hashes)
        if current_source_hashes(a.project)!=source_hashes:raise RuntimeError("native source files changed during benchmark")
        if {path.name:digest(path) for path in a.references.glob(source_case+"-*.npz")}!=report["reference_hashes"]:
            raise RuntimeError("reference file set changed during benchmark")
        report["input_and_source_bytes_unchanged"]=True
        stages=2 if ion else 4 if a.case=="free" else 2
        ct_times=checked_times(report["cantera_stage_seconds"],a.reps,stages)
        jl_times=checked_times(report["julia_stage_seconds"],a.reps,stages)
        report["cantera_seconds"]=ct_times.tolist();report["julia_seconds"]=jl_times.tolist()
        ratio=float(statistics.median(ct_times[1:])/statistics.median(jl_times[1:]))
        report["median_speed_ratio_cantera_over_julia"]=ratio
        physical=result.returncode==0 and report["accuracy"]["passed"]
        replay=report["cantera_replay"]["passed"] and report["julia_replay"]["passed"]
        report["ratio_threshold_met"]=ratio>=.95
        report["physical_gate_application"]="Last repetition of every source stage checked against refined references; strict numerical replay establishes the same outputs for all first/warm repetitions."
        report["performance_pass"]=bool(a.formal and physical and replay and ratio>=.95)
        report["passed"]=bool(physical and replay and (not a.formal or report["performance_pass"]))
    except Exception as error:
        report["failure"]=str(error)
    (a.output/"report.json").write_text(json.dumps(report,indent=2))
    print(json.dumps({key:report.get(key) for key in ["case","passed","median_speed_ratio_cantera_over_julia","failure"]}),flush=True)
    return 0 if report["passed"] else 1
if __name__=="__main__":raise SystemExit(main())
