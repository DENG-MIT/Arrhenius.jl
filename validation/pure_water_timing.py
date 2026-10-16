"""Paired native-Julia/Cantera water-example timing and output comparison.

Run this with Cantera 4 after pure_water_timing.jl writes its NPZ result. Both
programs exclude plotting, rendering, and file I/O from measured calculations.
Use --validate-only while other CPU-heavy work is active. A timed run requires
at least seven warm repetitions and is informational unless --qualification
controlled is set after coordinating an otherwise idle machine.

python validation/pure_water_timing.py --julia-result result-julia.npz \
    --output result-cantera.json --repetitions 9 --qualification informational
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import os
import platform
import statistics
import subprocess
import time
from benchmark_environment import host_metadata, matches_target, cantera_library_hashes

import_start = time.perf_counter()
import numpy as np
import cantera as ct
import cantera._cantera as compiled
import_seconds = time.perf_counter()-import_start

parser = argparse.ArgumentParser()
parser.add_argument("--julia-result", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--repetitions", type=int, default=9)
parser.add_argument("--qualification", choices=("informational","controlled"), default="informational")
parser.add_argument("--target", choices=("wsl", "apple-m4"), default="wsl")
parser.add_argument("--cantera-build-record", type=Path)
parser.add_argument("--validate-only", action="store_true")
args = parser.parse_args()
if args.repetitions < 7:
    parser.error("at least seven warm repetitions are required")


def rankine(temperature=300., pressure=8e5):
    w = ct.Water()
    w.TQ = temperature,0.
    h1,s1,p1 = w.h,w.s,w.P
    w.SP = s1,pressure
    pump_work = (w.h-h1)/.6
    w.HP = h1+pump_work,pressure
    h2 = w.h
    w.PQ = pressure,1.
    h3,s3 = w.h,w.s
    heat = h3-h2
    w.SP = s3,p1
    turbine_work = (h3-w.h)*.8
    w.HP = h3-turbine_work,p1
    return (pump_work,turbine_work,heat,(turbine_work-pump_work)/heat)


def vapordome():
    w = ct.Water()
    degc = np.hstack([np.array([w.min_temp-273.15,4,5,6,8]),np.arange(10,37),[38],
                      np.arange(40,100,5),np.arange(100,300,10),np.arange(300,380,20),
                      np.arange(370,374),[w.critical_temperature-273.15]])
    table = np.zeros((len(degc),14))
    table[:,0] = degc
    states = ct.SolutionArray(w,len(degc))
    # Preserve the source's vectorized Cantera calculation, not Python loops.
    states.TQ = degc+273.15,1.
    table[:,1] = states.P_sat/1e5
    table[:,4] = states.v
    table[:,7] = states.int_energy_mass/1e3
    table[:,10] = states.enthalpy_mass/1e3
    table[:,13] = states.entropy_mass/1e3
    states.TQ = degc+273.15,0.
    table[:,2] = states.v
    table[:,5] = states.int_energy_mass/1e3
    table[:,8] = states.enthalpy_mass/1e3
    table[:,11] = states.entropy_mass/1e3
    for difference,gas,liquid in ((3,4,2),(6,7,5),(9,10,8),(12,13,11)):
        table[:,difference] = table[:,gas]-table[:,liquid]
    w.TQ = w.min_temp,0.
    u0,h0,s0,pv0 = w.u/1e3,w.h/1e3,w.s/1e3,w.P*w.v/1e3
    table[:,(5,7)] -= u0
    table[:,(8,10)] -= h0-pv0
    table[:,(11,13)] -= s0
    return table


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


hardware = host_metadata()

native = np.load(args.julia_result)
native_meta = {key[:-5]: bytes(native[key]).decode() for key in native.files if key.endswith("_utf8")}
build_record = json.loads(args.cantera_build_record.read_text()) if args.cantera_build_record else None
source_sha = build_record["source"]["commit"] if build_record else getattr(ct,"__git_commit__","unknown")
libraries = cantera_library_hashes(ct.__file__)
native_host = {key: native_meta.get(key, "") for key in ("cpu", "system", "kernel_release")}
report = {
    "cantera_version":ct.__version__, "cantera_source_sha":source_sha,
    "cantera_build_record_sha256":sha(args.cantera_build_record) if args.cantera_build_record else None,
    "cantera_shared_libraries_sha256":libraries,
    "cantera_extension_sha256":sha(compiled.__file__), "harness_sha256":sha(__file__),
    "hardware":dict(hardware, load_average_start=os.getloadavg()),
    "benchmark_target": args.target,
    "julia_metadata":native_meta, "qualification": "validation_only" if args.validate_only else args.qualification,
    "import_seconds":import_seconds, "julia_import_seconds":float(native["import_seconds"][0]),
    "timer":"time.perf_counter / Julia time_ns; process import and first calls separate",
    "first_call_note":"first invocation of each case in process order; later cases may reuse compilation from earlier cases",
    "timed_scope":"fresh water model + source calculation; excludes plots, reports, dataframe presentation, file I/O",
    "minimum_speed_ratio":.95, "speed_ratio_definition":"median Cantera seconds / median Julia seconds",
    "timestamp_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
    "inputs":{"rankine":{"inlet_K":300.,"boiler_Pa":8e5,"pump_eta":.6,"turbine_eta":.8},
              "rankine_units":{"inlet_degF":80.33,"boiler_psi":116.03,"psi_to_Pa":6894.757293168364},
              "vapordome":{"temperature_count":74,"saturation_state_assignments":148,
                           "reference_state_assignments":1,"columns":["T","P","vf","vfg","vg","uf","ufg","ug","hf","hfg","hg","sf","sfg","sg"]}},
    "unit_conversion_limitation":"rankine_units uses SI-equivalent calculations; Julia Unitful and Python Pint wrapper overhead are excluded",
    "cases":{},
}
for name,run in (("rankine",rankine),("rankine_units",lambda:rankine((80.33-32)*5/9+273.15,116.03*6894.757293168364)),
                 ("vapordome",vapordome)):
    started = time.perf_counter()
    output = np.asarray(run())
    first_seconds = time.perf_counter()-started
    samples = []
    if not args.validate_only:
        run()
        gc.collect()
        for _ in range(args.repetitions):
            started = time.perf_counter()
            run()
            samples.append(time.perf_counter()-started)
    expected = native[f"{name}_output"]
    tolerance = 2e-6 if name != "vapordome" else 2e-7
    absolute = 1e-3 if name != "vapordome" else 1e-4
    matches = bool(np.allclose(output,expected,rtol=tolerance,atol=absolute))
    native_samples = native[f"{name}_seconds"].tolist()
    ratio = statistics.median(samples)/statistics.median(native_samples) if samples and native_samples else None
    controlled = (args.qualification == "controlled" and native_meta.get("qualification") == "controlled"
                  and not args.validate_only and len(samples)>=7 and len(native_samples)>=7
                  and matches_target(hardware,args.target) and matches_target(native_host,args.target)
                  and hardware["cpu"] == native_host["cpu"]
                  and hardware["kernel_release"] == native_host["kernel_release"]
                  and ct.__version__.startswith("4.0")
                  and source_sha == "726522be4e2a13454d8415b7ef799d621f665cf3")
    report["cases"][name] = {
        "correctness_pass":matches,"max_absolute_error":float(np.max(np.abs(output-expected))),
        "cantera_output":output.tolist(),"native_output":expected.tolist(),
        "cantera_first_seconds":first_seconds,"julia_first_seconds":float(native[f"{name}_first_seconds"][0]),
        "cantera_warm_seconds":samples,"julia_warm_seconds":native_samples,
        "speed_ratio":ratio,"performance_pass": bool(controlled and matches and ratio >= .95),
        "qualification": "controlled" if controlled else "not_qualified",
    }
    print(name,"correctness",matches,"speed ratio",ratio,"qualified",controlled,flush=True)
report["hardware"]["load_average_end"] = os.getloadavg()
args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
if not all(case["correctness_pass"] for case in report["cases"].values()):
    raise SystemExit("native and Cantera water example outputs differ")
