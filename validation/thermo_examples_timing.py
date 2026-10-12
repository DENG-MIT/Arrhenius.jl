"""Paired native-Julia/Cantera timing and output comparison for the complete
calculations in the Cantera examples thermo/equivalenceRatio.py,
thermo/isentropic.py, thermo/isentropic_units.py, thermo/sound_speed.py and
thermo/sound_speed_units.py (pinned source commit
726522be4e2a13454d8415b7ef799d621f665cf3).

Run this with Cantera 4 after validation/thermo_examples_timing.jl writes its
NPZ result. Both programs time only the complete example calculations on
prepared phases; mechanism loading, state reset, printing, plotting and file
I/O are excluded, and first calls are recorded separately from at least nine
warm batches. Samples are per-calculation means of batches targeting 20 ms,
capped at 1000 complete calls; reset and output checks are outside each call's
timer. The Cantera state is reset (TPX with the correct fresh
composition) before EVERY repetition; per-source state changes happen inside
the measured region. The units examples are recomputed with mathematically
equivalent explicit SI conversions (degF->K, m/s->ft/s) instead of the
source's pint wrappers, so unit-wrapper overhead is not measured; the
conversions are identical to the pinned sources. Use --validate-only while
other CPU-heavy work is active. A controlled result additionally requires an
otherwise idle target machine and the exact Cantera build record.

The sound_speed references keep the published source tolerances (rtol=1e-8 /
1e-6) as the timed baseline; those tolerances are NOT physical truth (the
default-tolerance SP setter leaves ~0.019 m/s and ~3.94 ft/s of numerical
error). Native correctness is therefore gated against independently refined
references (frozen SP entropy residual <1e-10 J/kg/K, equilibrate rtol=1e-13,
with the rtol=1e-11 refinement kept as a consistency check), and the
original-source outputs are reported against the refined reference separately.

python validation/thermo_examples_timing.py --julia-result result-julia.npz \
    --output result-cantera.json --repetitions 9 --qualification informational \
    --mechanism-dir /path/to/mechanisms \
    --cantera-build-record /path/to/wsl-cantera-build-record.json [--time-refined]
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import math
import os
import statistics
import time
from benchmark_environment import host_metadata, matches_target, cantera_library_hashes, verify_numerical_threads

import_start = time.perf_counter()
import numpy as np
import cantera as ct
import cantera._cantera as compiled
from isentropic_cases import nozzle_calculation, acoustics_calculation
from equivalence_ratio_case import equivalence_ratio_calculation
import_seconds = time.perf_counter()-import_start

CANTERA_COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
# sha256 of the pinned example source files at CANTERA_COMMIT; the workloads
# below were confirmed against these files, not reimplemented from memory.
EXAMPLE_SOURCE_SHA256 = {
    "thermo/equivalenceRatio.py": "613cdc6325f4a2ac2c411e0880500c47163628c19ecbf4ad77ca1fd537a5248e",
    "thermo/isentropic.py": "96529ca5ee8a1925d181026e7a4c3d2fc6f03832b500bdb715cf6411ab02e845",
    "thermo/isentropic_units.py": "a2a45fd4201d2592739347e9e5cf62ce388e57bd61c6509bcbdd758675accf81",
    "thermo/sound_speed.py": "1ee4e0c2370d66fe602b7c24917905d02abcf1ed108eb31896c872bcfb7d0133",
    "thermo/sound_speed_units.py": "44ff536a83c72c61e8e27e44c180d631e6e6c2637b004ebc5c11298aaa7a9275",
}
MIN_REPETITIONS = 9

parser = argparse.ArgumentParser()
parser.add_argument("--julia-result", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--repetitions", type=int, default=9)
parser.add_argument("--qualification", choices=("informational","controlled"), default="informational")
parser.add_argument("--target", choices=("wsl", "apple-m4"), default="wsl")
parser.add_argument("--cantera-build-record", type=Path)
parser.add_argument("--mechanism-dir", type=Path, required=True,
                    help="directory with h2o2.yaml, gri30.yaml, gri30_highT.yaml and their "
                         ".yaml.npz sidecars (as produced by validation/isentropic_cases.py)")
parser.add_argument("--validate-only", action="store_true")
parser.add_argument("--time-refined", action="store_true",
                    help="additionally time the refined acoustic workload (frozen entropy "
                         "refinement, equilibrate rtol=1e-13); recorded separately and never "
                         "used for the baseline speed ratio")
args = parser.parse_args()
if args.repetitions < MIN_REPETITIONS:
    parser.error(f"at least {MIN_REPETITIONS} warm repetitions are required")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utf8_meta(archive):
    return {key[:-5]: bytes(archive[key]).decode() for key in archive.files if key.endswith("_utf8")}


native = np.load(args.julia_result)
native_meta = utf8_meta(native)
mechanisms = {name: args.mechanism_dir/name
              for name in ("h2o2.yaml", "gri30.yaml", "gri30_highT.yaml")}
sidecars = {name: Path(str(path)+".npz") for name, path in mechanisms.items()}
local_hashes = {}
for name, path in mechanisms.items():
    if not path.is_file():
        raise SystemExit(f"missing mechanism {path}")
    if not sidecars[name].is_file():
        raise SystemExit(f"missing native sidecar {sidecars[name]} "
                         "(export it with mechanism/export_sidecar.py)")
    tag = {"h2o2.yaml": "h2o2", "gri30.yaml": "gri", "gri30_highT.yaml": "gri_highT"}[name]
    local_hashes[f"{tag}_mechanism_sha256"] = sha(path)
    local_hashes[f"{tag}_sidecar_sha256"] = sha(sidecars[name])
for key, digest in local_hashes.items():
    if native_meta.get(key) != digest:
        raise SystemExit(f"{key}: Julia-side record {native_meta.get(key)!r} != local file {digest}; "
                         "both harnesses must consume identical mechanism/sidecar files")
sidecar_meta = {}
for name, path in sidecars.items():
    meta = utf8_meta(np.load(path))
    tag = {"h2o2.yaml": "h2o2", "gri30.yaml": "gri", "gri30_highT.yaml": "gri_highT"}[name]
    if meta.get("source_sha256") != local_hashes[f"{tag}_mechanism_sha256"]:
        raise SystemExit(f"{path} was not exported from {mechanisms[name]} "
                         f"(sidecar source_sha256 {meta.get('source_sha256')!r})")
    sidecar_meta[name] = meta

gases = {name: ct.Solution(str(path)) for name, path in mechanisms.items()}
if gases["h2o2.yaml"].n_species != 10 or "N2" not in gases["h2o2.yaml"].species_names:
    raise SystemExit("the nozzle source uses the stock 10-species h2o2.yaml (with N2); "
                     "the bundled 9-species copy is not accepted")
for name in ("gri30.yaml", "gri30_highT.yaml"):
    if gases[name].n_species != 53:
        raise SystemExit(f"expected the 53-species {name} mechanism, got {gases[name].n_species}")
    if sidecar_meta[name].get("species_names") is not None and \
            sidecar_meta[name]["species_names"] != "\n".join(gases[name].species_names):
        raise SystemExit(f"{sidecars[name]} species do not match {mechanisms[name]}")
if sidecar_meta["h2o2.yaml"].get("species_names") is not None and \
        sidecar_meta["h2o2.yaml"]["species_names"] != "\n".join(gases["h2o2.yaml"].species_names):
    raise SystemExit("h2o2.yaml.npz species do not match h2o2.yaml")

# --- workload definitions (confirmed against the pinned sources) -------------
ACOUSTIC_COMPOSITION = "CH4:1,O2:2,N2:7.52"
NOZZLE_STATE = (1200., 10*ct.one_atm, "H2:1,N2:0.1")
SOUND_SPEED_TEMPERATURES = np.arange(300., 5001., 200.)  # 24 points, source rtol=1e-8
SOUND_SPEED_UNITS_TEMPERATURES = np.linspace(80., 4880., 25)  # degF, source rtol=1e-6
FRESH_MIXTURES = ("X_stoich_mole", "X_stoich_mass", "X_Z055", "X_fresh_burnt_case")
GRI_MW = gases["gri30.yaml"].molecular_weights


def reset_equivalence():
    # Source starts from the default 300 K / 1 atm state; every mixture call
    # fully defines the composition, and the burnt temperature persists after
    # the HP equilibration, so a full TPX reset is required before each run.
    gases["gri30.yaml"].TPX = 300., ct.one_atm, "CH4:1"


def reset_nozzle(mechanism):
    # The source's stagnation-state assignment is preparation; the timed
    # region contains the per-pressure SP state changes.
    gases[mechanism].TPX = NOZZLE_STATE


def reset_acoustics(mechanism):
    # The published sound_speed source carries each equilibrium composition
    # into the next TP state, so the fresh composition must be restored before
    # every repetition.
    gases[mechanism].TPX = 300., ct.one_atm, ACOUSTIC_COMPOSITION


def run_equivalence():
    mixture_order, states, scalars = equivalence_ratio_calculation(gases["gri30.yaml"])
    return {"mixture_order": mixture_order, "states": states, "scalars": scalars}


CASES = {
    "equivalence_ratio": {
        "setup": reset_equivalence,
        "run": run_equivalence,
        "timed_scope": ("9 set_equivalence_ratio/set_mixture_fraction mixture assignments with "
                        "phi/Z/Y scalar evaluations and one HP equilibration (source defaults) on "
                        "53 GRI species, in source order. Excludes phase construction, the 300 K / "
                        "1 atm TPX reset, printing, I/O."),
    },
    "isentropic": {
        "setup": lambda: reset_nozzle("h2o2.yaml"),
        "run": lambda: nozzle_calculation(gases["h2o2.yaml"], 200),
        "timed_scope": ("200 SP state assignments from the stagnation entropy with velocity, area "
                        "and Mach evaluation on stock 10-species h2o2. Excludes phase construction "
                        "and the stagnation TPX assignment, printing, plotting, I/O."),
    },
    "isentropic_units": {
        "setup": lambda: reset_nozzle("gri30.yaml"),
        "run": lambda: nozzle_calculation(gases["gri30.yaml"], 10, units=True),
        "timed_scope": ("10 SP state assignments on 53 GRI species with explicit SI unit handling "
                        "(mathematically equivalent to the source's pint wrappers; pint overhead "
                        "not measured). Excludes phase construction and the stagnation TPX "
                        "assignment, printing, I/O."),
    },
    "sound_speed": {
        "setup": lambda: reset_acoustics("gri30_highT.yaml"),
        "run": lambda: acoustics_calculation(gases["gri30_highT.yaml"],
                                             SOUND_SPEED_TEMPERATURES, 1e-8),
        "timed_scope": ("24 TP assignments with TP equilibration (source rtol=1e-8), 1.0001 "
                        "pressure perturbation, frozen SP state, SP equilibration and the built-in "
                        "sound speed on 53-species gri30_highT; equilibrium composition carried "
                        "between states per the source. Excludes phase construction and the fresh "
                        "TPX reset, printing, plotting, I/O."),
    },
    "sound_speed_units": {
        "setup": lambda: reset_acoustics("gri30.yaml"),
        "run": lambda: acoustics_calculation(gases["gri30.yaml"],
                                             SOUND_SPEED_UNITS_TEMPERATURES, 1e-6, units=True),
        "timed_scope": ("25 TP assignments (80..4880 degF via explicit degF->K conversion) with "
                        "TP equilibration (source rtol=1e-6), pressure perturbation, frozen SP "
                        "state, SP equilibration, built-in sound speed, explicit m/s->ft/s "
                        "conversion on 53 GRI species. Excludes phase construction and the fresh "
                        "TPX reset, printing, I/O."),
    },
}


# --- comparison gates (identical to the existing independent validators) -----
def entry(got, want, rtol, atol, note, extra=None):
    got, want = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    if got.shape != want.shape:
        raise SystemExit(f"shape {got.shape} != native {want.shape}")
    error = np.abs(got-want)
    record = {"correctness_pass": bool(np.allclose(got, want, rtol=rtol, atol=atol)),
              "rtol": rtol, "atol": atol, "tolerance_note": note,
              "max_absolute_error": float(error.max()),
              "max_relative_error": float((error/np.maximum(np.abs(want), max(atol, 1e-300))).max()),
              "cantera_output": got.tolist(), "native_output": want.tolist()}
    record.update(extra or {})
    return record


def compare_equivalence(output):
    native_T = {"fresh": float(native["equivalence_ratio_T_fresh_K"][0]),
                "burnt": float(native["equivalence_ratio_T_post_burnt_K"][0])}
    outputs, case_ok = {}, True
    for name in output["mixture_order"]:
        T_c, P_c, X_c, Y_c = output["states"][name]
        X_n = np.asarray(native[f"equivalence_ratio_{name}"], dtype=float)
        if name == "X_burnt":
            # Native equilibrium tolerances; see validation/equilibrium_cases.jl.
            # Cantera's default HP solve drifts by 8.5e-6 Pa here: the native
            # pressure is exactly 1 atm and the reference is required only
            # within its 1e-9 solver tolerance.
            Y_n = np.asarray(native["equivalence_ratio_Y_burnt"], dtype=float)
            T_n = float(native["equivalence_ratio_T_burnt_K"][0])
            P_n = float(native["equivalence_ratio_P_burnt_Pa"][0])
            ok = (abs(T_n-T_c) <= 2e-7*T_c+2e-5 and P_n == ct.one_atm
                  and abs(P_n-P_c) <= 1e-9*P_c
                  and np.abs(X_n-X_c).max() < 1e-7 and np.abs(Y_n-Y_c).max() < 1e-7)
            outputs[name] = entry(X_c, X_n, 0., 1e-7,
                                  "HP-equilibrated burnt state: X/Y max abs < 1e-7, T within "
                                  "rtol 2e-7/atol 2e-5 K, native P exactly 1 atm within the "
                                  "reference's 1e-9 solver tolerance",
                                  {"temperature_pass": bool(abs(T_n-T_c) <= 2e-7*T_c+2e-5),
                                   "pressure_pass": bool(P_n == ct.one_atm
                                                         and abs(P_n-P_c) <= 1e-9*P_c),
                                   "max_mass_fraction_error": float(np.abs(Y_n-Y_c).max())})
            case_ok &= ok
        else:
            Y_n = X_n*GRI_MW/np.dot(X_n, GRI_MW)
            T_n = native_T["fresh"] if name in FRESH_MIXTURES else native_T["burnt"]
            pressure_rtol = 1e-13 if name in FRESH_MIXTURES else 1e-9
            ok = (np.allclose(X_n, X_c, rtol=1e-12, atol=1e-14)
                  and np.allclose(Y_n, Y_c, rtol=1e-12, atol=1e-14)
                  and abs(T_n-T_c) <= 2e-7*T_c+2e-5
                  and abs(ct.one_atm-P_c) <= pressure_rtol*P_c)
            outputs[name] = entry(X_c, X_n, 1e-12, 1e-14,
                                  "mixture composition from set_equivalence_ratio/"
                                  "set_mixture_fraction; round-off agreement",
                                  {"max_mass_fraction_error": float(np.abs(Y_n-Y_c).max()),
                                   "temperature_pass": bool(abs(T_n-T_c) <= 2e-7*T_c+2e-5),
                                   "pressure_pass": bool(abs(ct.one_atm-P_c) <= pressure_rtol*P_c)})
            case_ok &= ok
    for name, value in output["scalars"].items():
        tolerance = 1e-8 if name in ("phi_burnt", "Z_burnt") else 1e-9
        got = float(native[f"equivalence_ratio_{name}"][0])
        record = entry(value, got, tolerance, 1e-12,
                       "printed scalar; equivalence-ratio/mixture-fraction round-off agreement"
                       + (" (burnt composition, relaxed for the HP solve)" if tolerance == 1e-8
                          else ""))
        outputs[f"scalar_{name}"] = record
        case_ok &= record["correctness_pass"]
    return outputs, case_ok


NOZZLE_GATES = {
    "data": (2e-8, 1e-8, "area ratio, Mach number, temperature ratio (K for the units case) and "
                        "pressure ratio; SP-state round-off agreement"),
    "states": (2e-8, 1e-8, "T, P, density, enthalpy, entropy rows (SI); SP-state round-off agreement"),
    "throat_area": (2e-8, 1e-8, "minimum nozzle area (m^2)"),
}


def compare_nozzle(case, output):
    outputs, case_ok = {}, True
    for key, (rtol, atol, note) in NOZZLE_GATES.items():
        record = entry(output[key], native[f"{case}_{key}"], rtol, atol, note)
        outputs[key] = record
        case_ok &= record["correctness_pass"]
    return outputs, case_ok


def acoustic_reference(mechanism, temperatures, rtol, units):
    # Independent correctness reference; computed outside any timed region on
    # its own phase. refine_entropy tightens the frozen isentrope to an
    # entropy residual <1e-10 J/kg/K before the density difference is taken.
    gas = ct.Solution(str(mechanisms[mechanism]))
    gas.TPX = 300., ct.one_atm, ACOUSTIC_COMPOSITION
    return acoustics_calculation(gas, temperatures, rtol, units=units, refine_entropy=True)


def compare_acoustics(case, output, refined, refined2):
    # Correctness is gated against the refined2 reference (frozen SP entropy
    # residual <1e-10 J/kg/K, equilibrate rtol=1e-13); the published source
    # tolerance is reported separately and is not treated as physical truth.
    outputs, case_ok = {}, True
    data_n = np.asarray(native[f"{case}_data"], dtype=float)
    states_n = np.asarray(native[f"{case}_final_states"], dtype=float)
    record = entry(refined2["data"], data_n, 1e-7, 1e-9,
                   "native vs refined reference (frozen entropy residual <1e-10 J/kg/K, "
                   "equilibrate rtol=1e-13); pressure finite differences amplify the "
                   "equilibrium tolerance",
                   {"reference": "refined2",
                    "cantera_output": np.asarray(refined2["data"]).tolist(),
                    "original_source_output": np.asarray(output["data"]).tolist(),
                    "original_source_max_error_vs_refined":
                        float(np.abs(output["data"]-refined2["data"]).max()),
                    "refinement_rtol_1e_11_vs_1e_13_max_change":
                        float(np.abs(refined["data"]-refined2["data"]).max()),
                    "refinement_consistency_pass":
                        bool(np.allclose(refined["data"], refined2["data"], rtol=1e-6, atol=1e-9)),
                    "published_agreement_at_1e_4":
                        bool(np.allclose(output["data"], data_n, rtol=1e-4, atol=1e-9))})
    outputs["data"] = record
    case_ok &= record["correctness_pass"] and record["refinement_consistency_pass"]
    top = entry(refined2["final_states"][:2], states_n[:2], 2e-8, 1e-8,
                "perturbed equilibrium T, P rows of the carried final states",
                {"original_source_max_error_vs_refined":
                     float(np.abs(output["final_states"][:2]-refined2["final_states"][:2]).max())})
    composition = entry(refined2["final_states"][2:], states_n[2:], 2e-7, 5e-9,
                        "carried equilibrium compositions of the final states",
                        {"original_source_max_error_vs_refined":
                             float(np.abs(output["final_states"][2:]
                                          -refined2["final_states"][2:]).max())})
    outputs["final_states_T_P"] = top
    outputs["final_states_X"] = composition
    case_ok &= top["correctness_pass"] and composition["correctness_pass"]
    return outputs, case_ok


def assert_same_equivalence(first, repeated):
    for name in first["mixture_order"]:
        T1, P1, X1, Y1 = first["states"][name]
        T2, P2, X2, Y2 = repeated["states"][name]
        assert T1 == T2 and P1 == P2, f"equivalence_ratio reset changed the {name} state"
        np.testing.assert_array_equal(X1, X2, err_msg=f"equivalence_ratio reset changed {name}")
        np.testing.assert_array_equal(Y1, Y2, err_msg=f"equivalence_ratio reset changed {name}")
    for name, value in first["scalars"].items():
        assert repeated["scalars"][name] == value, f"equivalence_ratio reset changed {name}"


def assert_same_output(case, first, repeated):
    if case == "equivalence_ratio":
        assert_same_equivalence(first, repeated)
        return
    for key in ("data", "states", "throat_area", "final_states"):
        if key in first:
            np.testing.assert_allclose(repeated[key], first[key], rtol=1e-9, atol=1e-12,
                                       err_msg=f"{case} reset changed {key}")


hardware = host_metadata()
numerical_threads = verify_numerical_threads(set_accelerate=True)
build_record = json.loads(args.cantera_build_record.read_text()) if args.cantera_build_record else None
source_sha = build_record["source"]["commit"] if build_record else getattr(ct, "__git_commit__", "unknown")
loaded_libraries = cantera_library_hashes(ct.__file__)
extension_sha = sha(compiled.__file__)
recorded_libraries = {Path(k).name: v for k, v in (build_record or {}).get("library_hashes", {}).items()}
build_matches = bool(loaded_libraries and recorded_libraries
                     and all(recorded_libraries.get(k) == v for k, v in loaded_libraries.items())
                     and recorded_libraries.get(Path(compiled.__file__).name) == extension_sha)
thread_settings = {key: os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")}
if args.target == "apple-m4":
    thread_settings["VECLIB_MAXIMUM_THREADS"] = os.environ.get("VECLIB_MAXIMUM_THREADS")
native_host = {key: native_meta.get(key, "") for key in ("cpu", "system", "kernel_release")}
report = {
    "cantera_version": ct.__version__, "cantera_source_sha": source_sha,
    "cantera_build_record_sha256": sha(args.cantera_build_record) if args.cantera_build_record else None,
    "cantera_shared_libraries_sha256": loaded_libraries,
    "cantera_extension_sha256": extension_sha, "harness_sha256": sha(__file__),
    "calculation_helpers_sha256": {name: sha(Path(__file__).with_name(name))
                                   for name in ("isentropic_cases.py", "equivalence_ratio_case.py", "benchmark_environment.py")},
    "loaded_libraries_match_build_record": build_matches,
    "thread_settings": dict(thread_settings, julia_threads=int(native["julia_threads"][0]),
                            julia_blas_threads=int(native["blas_threads"][0])),
    "numerical_threads": numerical_threads,
    "cantera_examples": {"commit": CANTERA_COMMIT, "source_sha256": EXAMPLE_SOURCE_SHA256},
    "hardware": dict(hardware, load_average_start=os.getloadavg()),
    "benchmark_target": args.target,
    "julia_metadata": native_meta,
    "mechanisms": {name: {"file": mechanisms[name].name,
                          "sha256": local_hashes[f"{tag}_mechanism_sha256"],
                          "sidecar_sha256": local_hashes[f"{tag}_sidecar_sha256"],
                          "sidecar_exported_from_sha256": sidecar_meta[name].get("source_sha256"),
                          "species": gases[name].n_species}
                   for name, tag in (("h2o2.yaml", "h2o2"), ("gri30.yaml", "gri"),
                                     ("gri30_highT.yaml", "gri_highT"))},
    "qualification": "validation_only" if args.validate_only else args.qualification,
    "import_seconds": import_seconds, "julia_import_seconds": float(native["import_seconds"][0]),
    "timer": "time.perf_counter / Julia time_ns; process import and first calls separate",
    "first_call_note": "first invocation of each case in process order; the Julia first call includes JIT compilation",
    "minimum_warm_repetitions": MIN_REPETITIONS,
    "minimum_speed_ratio": .95, "speed_ratio_definition": "median Cantera seconds / median Julia seconds",
    "units_note": ("units examples use mathematically equivalent explicit SI conversions "
                   "(degF->K, m/s->ft/s) instead of the source's pint wrappers; the timed "
                   "physics is identical and pint overhead is not measured"),
    "reset_policy": ("Cantera: TPX reset with the correct fresh composition before EVERY "
                     "repetition (300 K / 1 atm / CH4:1 for equivalence_ratio, the 1200 K / "
                     "10 atm / H2:1,N2:0.1 stagnation state for the nozzles, 300 K / 1 atm / "
                     "CH4:1,O2:2,N2:7.52 for the acoustic cases, whose source carries the "
                     "equilibrium composition forward); state changes prescribed by the sources "
                     "happen inside the measured region. Native: stateless inputs, nothing to "
                     "reset; native acoustic SP equilibration uses property_rtol=1e-13."),
    "correctness_reference_note": ("Nozzle and equivalence_ratio outputs are compared with the "
                                   "original-source Cantera calculation at source tolerances. "
                                   "Acoustic outputs are gated against the independently refined "
                                   "reference (frozen SP entropy residual <1e-10 J/kg/K, "
                                   "equilibrate rtol=1e-13, with rtol=1e-11 as a consistency "
                                   "check); the published source tolerance is not physical "
                                   "truth, so original-source vs refined errors are reported "
                                   "separately and never used as a pass gate for native output."),
    "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "cases": {},
}

controlled_base = (args.qualification == "controlled" and native_meta.get("qualification") == "controlled"
                   and not args.validate_only
                   and matches_target(hardware, args.target) and matches_target(native_host, args.target)
                   and hardware["cpu"] == native_host["cpu"]
                   and hardware["kernel_release"] == native_host["kernel_release"]
                   and ct.__version__.startswith("4.0")
                   and build_record is not None
                   and source_sha == CANTERA_COMMIT and build_matches
                   and int(native["julia_threads"][0]) == 1 and int(native["blas_threads"][0]) == 1
                   and all(value == "1" for value in thread_settings.values()))

# Refined acoustic correctness references (computed once, outside timing).
references = {}
for case, mechanism, temperatures, units in (
        ("sound_speed", "gri30_highT.yaml", SOUND_SPEED_TEMPERATURES, False),
        ("sound_speed_units", "gri30.yaml", SOUND_SPEED_UNITS_TEMPERATURES, True)):
    references[case] = (acoustic_reference(mechanism, temperatures, 1e-11, units),
                        acoustic_reference(mechanism, temperatures, 1e-13, units))

any_failure = False
for name, spec in CASES.items():
    spec["setup"]()
    started = time.perf_counter()
    output = spec["run"]()
    first_seconds = time.perf_counter()-started
    samples = []
    batch_size = 0
    if not args.validate_only:
        spec["setup"]()
        started = time.perf_counter()
        warmup = spec["run"]()
        warmup_seconds = time.perf_counter()-started
        assert_same_output(name, output, warmup)
        batch_size = max(1, min(1000, math.ceil(.020/max(warmup_seconds, 1e-9))))
        gc.collect()
        for _ in range(args.repetitions):
            elapsed = 0.0
            for _ in range(batch_size):
                spec["setup"]()
                started = time.perf_counter()
                repeated_output = spec["run"]()
                elapsed += time.perf_counter()-started
                # Check outside the measured region: a fast, mutated workload
                # is not a valid repetition of the original calculation.
                assert_same_output(name, output, repeated_output)
            samples.append(elapsed/batch_size)
    if name == "equivalence_ratio":
        outputs, case_ok = compare_equivalence(output)
    elif name in ("isentropic", "isentropic_units"):
        outputs, case_ok = compare_nozzle(name, output)
    else:
        refined, refined2 = references[name]
        outputs, case_ok = compare_acoustics(name, output, refined, refined2)
    native_samples = native[f"{name}_seconds"].tolist()
    native_batch_size = int(native[f"{name}_batch_size"][0])
    native_checked = int(native[f"{name}_warm_outputs_checked"][0])
    if native_checked != len(native_samples)*native_batch_size:
        raise SystemExit(f"{name}: inconsistent native output-check count")
    ratio = statistics.median(samples)/statistics.median(native_samples) if samples and native_samples else None
    controlled = bool(controlled_base and len(samples) >= MIN_REPETITIONS
                      and len(native_samples) >= MIN_REPETITIONS
                      and 1 <= native_batch_size <= 1000)
    case_report = {
        "correctness_pass": bool(case_ok), "outputs": outputs, "timed_scope": spec["timed_scope"],
        "cantera_first_seconds": first_seconds, "julia_first_seconds": float(native[f"{name}_first_seconds"][0]),
        "cantera_warm_seconds": samples, "julia_warm_seconds": native_samples,
        "sample_definition": "mean timed seconds per complete calculation within each warm batch; resets and output checks excluded",
        "batch_policy": "ceil(0.020 / warmup seconds), clamped to 1..1000 complete calculations",
        "cantera_batch_size": batch_size, "julia_batch_size": native_batch_size,
        "warm_outputs_checked": len(samples)*batch_size, "native_warm_outputs_checked": native_checked,
        "speed_ratio": ratio, "performance_pass": bool(controlled and case_ok and ratio is not None and ratio >= .95),
        "qualification": "controlled" if controlled else "not_qualified",
    }
    if args.time_refined and name in references and not args.validate_only:
        # Clearly distinct refined-workload timing; never part of the baseline
        # ratio. Same reset protocol, refined entropy + rtol=1e-13 workload.
        mechanism = "gri30_highT.yaml" if name == "sound_speed" else "gri30.yaml"
        temperatures = SOUND_SPEED_TEMPERATURES if name == "sound_speed" else SOUND_SPEED_UNITS_TEMPERATURES
        units = name != "sound_speed"
        run_refined = lambda: acoustics_calculation(gases[mechanism], temperatures, 1e-13,
                                                    units=units, refine_entropy=True)
        reset = lambda: reset_acoustics(mechanism)
        reset()
        refined_output = run_refined()  # untimed warmup on the timed phase
        refined_samples = []
        for _ in range(args.repetitions):
            reset()
            started = time.perf_counter()
            repeated = run_refined()
            refined_samples.append(time.perf_counter()-started)
            np.testing.assert_allclose(repeated["data"], refined_output["data"],
                                       rtol=1e-9, atol=1e-12,
                                       err_msg=f"{name} refined reset changed data")
        case_report["refined_workload"] = {
            "description": "frozen entropy refinement (<1e-10 J/kg/K) with equilibrate "
                           "rtol=1e-13; distinct from the original-source baseline",
            "cantera_warm_seconds": refined_samples}
    report["cases"][name] = case_report
    any_failure |= not case_ok
    print(name, "correctness", case_ok, "speed ratio", ratio, "qualified", controlled, flush=True)
report["hardware"]["load_average_end"] = os.getloadavg()
report["numerical_threads_after"] = verify_numerical_threads()
args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
if any_failure:
    raise SystemExit("native and Cantera thermo example outputs differ")
