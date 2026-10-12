"""Paired native-Julia/Cantera timing and output comparison for the complete
calculations in the Cantera examples kinetics/blowers_masel.py and
transport/dusty_gas.py (pinned source commit
726522be4e2a13454d8415b7ef799d621f665cf3).

Run this with Cantera 4 after validation/property_examples_timing.jl writes its
NPZ result. Both programs time only the complete example calculations on
prepared phases/workspaces; mechanism loading, workspace construction, H-thermo
restoration, printing, plotting, and file I/O are excluded, and first calls are
recorded separately from at least nine warm repetitions. Use --validate-only
while other CPU-heavy work is active. A controlled result additionally requires
an otherwise idle target machine and the exact Cantera build record.

python validation/property_examples_timing.py --julia-result result-julia.npz \
    --output result-cantera.json --repetitions 9 --qualification informational \
    --h2o2-mechanism /path/to/h2o2.yaml --h2o2-sidecar /path/to/h2o2.multicomponent.npz \
    --cantera-build-record /path/to/wsl-cantera-build-record.json
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import os
import statistics
import time
from benchmark_environment import host_metadata, matches_target, cantera_library_hashes

import_start = time.perf_counter()
import numpy as np
import cantera as ct
import cantera._cantera as compiled
import_seconds = time.perf_counter()-import_start

CANTERA_COMMIT = "726522be4e2a13454d8415b7ef799d621f665cf3"
# sha256 of the pinned example source files at CANTERA_COMMIT; the workloads
# below were confirmed against these files, not reimplemented from memory.
EXAMPLE_SOURCE_SHA256 = {
    "kinetics/blowers_masel.py": "b79f3f8f7c5c0e59454ab2fc8139cc468c5b97b144e55b28e3c0c0b14b6f5714",
    "transport/dusty_gas.py": "0d3c75ba6debd34eafcb2609dac7da1224a97c13694ebeb6bd6ae84589bd487f",
}
MIN_REPETITIONS = 9

parser = argparse.ArgumentParser()
parser.add_argument("--julia-result", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--repetitions", type=int, default=9)
parser.add_argument("--qualification", choices=("informational","controlled"), default="informational")
parser.add_argument("--target", choices=("wsl", "apple-m4"), default="wsl")
parser.add_argument("--cantera-build-record", type=Path)
parser.add_argument("--gri-mechanism", type=Path,
                    default=Path(__file__).parents[1]/"mechanism"/"gri30.yaml")
parser.add_argument("--h2o2-mechanism", type=Path, required=True,
                    help="stock 10-species h2o2.yaml (with N2); the bundled 9-species copy is rejected")
parser.add_argument("--h2o2-sidecar", type=Path, required=True,
                    help="multicomponent sidecar exported from --h2o2-mechanism with mechanism/export_multicomponent.py")
parser.add_argument("--validate-only", action="store_true")
args = parser.parse_args()
if args.repetitions < MIN_REPETITIONS:
    parser.error(f"at least {MIN_REPETITIONS} warm repetitions are required")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utf8_meta(archive):
    return {key[:-5]: bytes(archive[key]).decode() for key in archive.files if key.endswith("_utf8")}


native = np.load(args.julia_result)
native_meta = utf8_meta(native)
local_hashes = {"gri_mechanism_sha256": sha(args.gri_mechanism),
                "h2o2_mechanism_sha256": sha(args.h2o2_mechanism),
                "h2o2_sidecar_sha256": sha(args.h2o2_sidecar)}
for key, digest in local_hashes.items():
    if native_meta.get(key) != digest:
        raise SystemExit(f"{key}: Julia-side record {native_meta.get(key)!r} != local file {digest}; "
                         "both harnesses must consume identical mechanism/sidecar files")
sidecar = np.load(args.h2o2_sidecar)
sidecar_meta = utf8_meta(sidecar)
if sidecar_meta.get("source_sha256") != local_hashes["h2o2_mechanism_sha256"]:
    raise SystemExit("--h2o2-sidecar was not exported from --h2o2-mechanism "
                     f"(sidecar source_sha256 {sidecar_meta.get('source_sha256')!r})")

gri = ct.Solution(str(args.gri_mechanism))
if gri.n_species != 53:
    raise SystemExit(f"expected the 53-species GRI-Mech 3.0 mechanism, got {gri.n_species}")
h2o2 = ct.Solution(str(args.h2o2_mechanism))
if h2o2.n_species != 10 or "N2" not in h2o2.species_names:
    raise SystemExit("dusty_gas.py uses the stock 10-species h2o2.yaml (with N2); "
                     "the bundled 9-species copy is not accepted")
if sidecar_meta.get("species_names") != "\n".join(h2o2.species_names):
    raise SystemExit("--h2o2-sidecar species do not match --h2o2-mechanism")

# --- blowers_masel.py workload (confirmed against the pinned source) ---------
BM_PARAMETERS = (3.87e1, 2.7, 6260*1000*4.184)
bm_gas = ct.Solution(thermo="ideal-gas", kinetics="gas", species=gri.species(),
                     reactions=[ct.Reaction(equation="O + H2 <=> H + OH", rate=ct.ArrheniusRate(*BM_PARAMETERS)),
                                ct.Reaction(equation="O + H2 <=> H + OH", rate=ct.BlowersMaselRate(*BM_PARAMETERS, 1e9)),
                                ct.Reaction(equation="H + CH4 <=> CH3 + H2", rate=ct.BlowersMaselRate(*BM_PARAMETERS, 1e9))])
bm_gas.TP = 300., ct.one_atm
bm_h_index = bm_gas.species_index("H")
bm_h_thermo = bm_gas.species(bm_h_index).thermo
bm_h_coefficients0 = bm_h_thermo.coeffs.copy()
bm_h_bounds = (bm_h_thermo.min_temp, bm_h_thermo.max_temp, bm_h_thermo.reference_pressure)


def reset_bm():
    # The enthalpy sweep rewrites the H NASA thermo, so restore the initial
    # thermo and the source's 300 K state before EVERY repetition.
    # Species wrappers alias the phase's mutable species. Keep a value snapshot
    # of the coefficients and reconstruct thermo instead of saving a wrapper.
    species = bm_gas.species(bm_h_index)
    species.thermo = ct.NasaPoly2(*bm_h_bounds, bm_h_coefficients0.copy())
    bm_gas.modify_species(bm_h_index, species)
    bm_gas.reaction(1).rate.delta_enthalpy = 0.
    bm_gas.TP = 300., ct.one_atm


def bm_calculation(gas):
    temperatures = np.arange(300., 3500., 100.)
    rates = np.zeros((3, temperatures.size))
    for j, T in enumerate(temperatures):
        gas.TP = T, ct.one_atm
        rates[:, j] = gas.forward_rate_constants
    # The standalone rate's delta_enthalpy is initially zero; forward-rate
    # evaluation does not change it. Thus this is the intrinsic barrier.
    E0 = gas.reaction(1).rate.activation_energy
    enthalpies = np.linspace(-5*E0, 5*E0, 100)
    barriers = np.zeros(enthalpies.size)
    for k, desired in enumerate(enthalpies):
        species = gas.species(bm_h_index)
        coefficients = species.thermo.coeffs.copy()
        change = (desired-gas.delta_enthalpy[1])/ct.gas_constant
        coefficients[6] += change
        coefficients[13] += change
        species.thermo = ct.NasaPoly2(species.thermo.min_temp, species.thermo.max_temp,
                                      species.thermo.reference_pressure, coefficients)
        gas.modify_species(bm_h_index, species)
        rate = gas.reaction(1).rate
        rate.delta_enthalpy = gas.delta_enthalpy[1]
        barriers[k] = rate.activation_energy
    return {"temperatures": temperatures, "rates": rates,
            "enthalpies": enthalpies, "barriers": barriers}


# --- dusty_gas.py workload (confirmed against the pinned source) -------------
DUSTY_COMPOSITION = {"OH": 1., "H": 2., "O2": 3., "O": 1e-8, "H2": 1e-8,
                     "H2O": 1e-8, "H2O2": 1e-8, "HO2": 1e-8, "AR": 1e-8}
dusty_holder = {}


def build_dusty():
    # Fresh transport object per repetition (untimed preparation) so every
    # timed call runs the complete calculation with cold transport caches; see
    # validation/dusty_gas_cache_probe.py for the CT4 cache-history issue.
    g = ct.DustyGas(str(args.h2o2_mechanism))
    g.porosity = .2
    g.tortuosity = 4.
    g.mean_pore_radius = 1.5e-7
    g.mean_particle_diameter = 1.5e-6
    dusty_holder["g"] = g


def dusty_calculation():
    g = dusty_holder["g"]
    g.TPX = 500., ct.one_atm, DUSTY_COMPOSITION
    diffusion = g.multi_diff_coeffs.copy()
    conductivity = g.thermal_conductivity
    T1, rho1, Y1 = g.TDY
    g.TP = T1, 1.2*ct.one_atm
    T2, rho2, Y2 = g.TDY
    delta = .001
    zero = g.molar_fluxes(T1, T1, rho1, rho1, Y1, Y1, delta).copy()
    pressure = g.molar_fluxes(T1, T2, rho1, rho2, Y1, Y2, delta).copy()
    return {"diffusion": diffusion, "conductivity": np.atleast_1d(conductivity),
            "zero_flux": zero, "pressure_flux": pressure}


COMPARISONS = {
    "blowers_masel": (
        ("temperatures", 1e-12, 1e-12,
         "identical 300:100:3400 K grid built by arange/range; agree to round-off"),
        ("rates", 2e-12, 1e-9,
         "bimolecular forward rate constants (m^3/kmol/s); tolerance matches the independent rate validation"),
        ("enthalpies", 2e-12, 1e-6,
         "prescribed reaction enthalpies (J/kmol) spanning +/-5 times the intrinsic barrier"),
        ("barriers", 2e-12, 1e-6,
         "activation energies (J/kmol); absolute floor covers cancellation near the zero barrier"),
    ),
    "dusty_gas": (
        ("diffusion", 1e-10, 1e-30,
         "full diffusion matrix (m^2/s), including trace entries near 1e-26; all entries are checked"),
        ("conductivity", 1e-10, 1e-14,
         "thermal conductivity (W/m/K); round-off agreement with independently exported transport fits"),
        ("zero_flux", 0., 0.,
         "zero-gradient molar fluxes (kmol/m^2/s), identically zero by construction on both sides"),
        ("pressure_flux", 1e-10, 1e-30,
         "pressure-gradient molar fluxes (kmol/m^2/s), including the trace N2 floor"),
    ),
}

TIMED_SCOPE = {
    "blowers_masel": ("32 gas.TP state assignments with forward rate constants for the 3 example "
                      "reactions on 53 GRI species, intrinsic barrier, then the 100-step "
                      "enthalpy sweep rewriting H NASA thermo (source behavior, not a Python-loop "
                      "replacement of a vectorized API). Excludes phase/mechanism construction, "
                      "H-thermo/state restoration, printing, plotting, I/O."),
    "dusty_gas": ("TPX state assignment, full multi_diff_coeffs matrix, thermal_conductivity, TDY "
                  "quantities, 1.2-pressure state assignment with TDY, zero-gradient and "
                  "pressure-gradient molar_fluxes. Excludes DustyGas object construction and porous "
                  "medium parameter assignment (fresh per repetition), printing, I/O."),
}

hardware = host_metadata()
build_record = json.loads(args.cantera_build_record.read_text()) if args.cantera_build_record else None
source_sha = build_record["source"]["commit"] if build_record else getattr(ct, "__git_commit__", "unknown")
loaded_libraries = cantera_library_hashes(ct.__file__)
extension_sha = sha(compiled.__file__)
recorded_libraries = {Path(k).name: v for k, v in (build_record or {}).get("library_hashes", {}).items()}
build_matches = bool(loaded_libraries and recorded_libraries
                     and all(recorded_libraries.get(k) == v for k, v in loaded_libraries.items())
                     and recorded_libraries.get(Path(compiled.__file__).name) == extension_sha)
thread_settings = {key: os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")}
if hardware["system"] == "Darwin":
    thread_settings["VECLIB_MAXIMUM_THREADS"] = os.environ.get("VECLIB_MAXIMUM_THREADS")
native_host = {key: native_meta.get(key, "") for key in ("cpu", "system", "kernel_release")}
report = {
    "cantera_version": ct.__version__, "cantera_source_sha": source_sha,
    "cantera_build_record_sha256": sha(args.cantera_build_record) if args.cantera_build_record else None,
    "cantera_shared_libraries_sha256": loaded_libraries,
    "cantera_extension_sha256": extension_sha, "harness_sha256": sha(__file__),
    "loaded_libraries_match_build_record": build_matches,
    "thread_settings": dict(thread_settings, julia_threads=int(native["julia_threads"][0]),
                            julia_blas_threads=int(native["blas_threads"][0])),
    "cantera_examples": {"commit": CANTERA_COMMIT, "source_sha256": EXAMPLE_SOURCE_SHA256},
    "hardware": dict(hardware, load_average_start=os.getloadavg()),
    "benchmark_target": args.target,
    "julia_metadata": native_meta,
    "mechanisms": {
        "gri": {"file": args.gri_mechanism.name, "sha256": local_hashes["gri_mechanism_sha256"], "species": 53},
        "h2o2": {"file": args.h2o2_mechanism.name, "sha256": local_hashes["h2o2_mechanism_sha256"], "species": 10},
        "h2o2_sidecar": {"file": args.h2o2_sidecar.name, "sha256": local_hashes["h2o2_sidecar_sha256"],
                         "exported_from_sha256": sidecar_meta.get("source_sha256"),
                         "export_cantera_version": sidecar_meta.get("cantera_version")},
    },
    "qualification": "validation_only" if args.validate_only else args.qualification,
    "import_seconds": import_seconds, "julia_import_seconds": float(native["import_seconds"][0]),
    "timer": "time.perf_counter / Julia time_ns; process import and first calls separate",
    "first_call_note": "first invocation of each case in process order; the Julia first call includes JIT compilation",
    "minimum_warm_repetitions": MIN_REPETITIONS,
    "minimum_speed_ratio": .95, "speed_ratio_definition": "median Cantera seconds / median Julia seconds",
    "reset_policy": ("Cantera blowers_masel: initial H NASA thermo and 300 K/1 atm state restored before "
                     "every repetition because the sweep mutates species thermo. Cantera dusty_gas: fresh "
                     "DustyGas object per repetition (cold caches). Native: blowers_masel_calculation "
                     "uses a stateless Solution and evaluates the prescribed reaction-enthalpy shift "
                     "directly; dusty_gas_calculation runs on a fresh workspace per repetition and copies "
                     "every returned array out of workspace buffers."),
    "native_scope_note": ("native blowers_masel_calculation prescribes the reaction-enthalpy shift directly "
                          "in activation_energy/rate_constant, physically equivalent to Cantera's H-thermo "
                          "rewrite; the native side therefore skips the per-step species thermo mutation"),
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

any_failure = False
for name, setup, run in (("blowers_masel", reset_bm, lambda: bm_calculation(bm_gas)),
                         ("dusty_gas", build_dusty, dusty_calculation)):
    setup()
    started = time.perf_counter()
    output = run()
    first_seconds = time.perf_counter()-started
    samples = []
    if not args.validate_only:
        setup()
        warmup = run()
        for key, rtol, atol, _ in COMPARISONS[name]:
            np.testing.assert_allclose(warmup[key], output[key], rtol=rtol, atol=atol,
                                       err_msg=f"{name} reset changed {key}")
        gc.collect()
        for _ in range(args.repetitions):
            setup()
            started = time.perf_counter()
            repeated_output = run()
            samples.append(time.perf_counter()-started)
            # Check outside the measured region: a fast, mutated workload is
            # not a valid repetition of the original calculation.
            for key, rtol, atol, _ in COMPARISONS[name]:
                np.testing.assert_allclose(repeated_output[key], output[key], rtol=rtol, atol=atol,
                                           err_msg=f"{name} repetition changed {key}")
    outputs = {}
    case_ok = True
    for key, rtol, atol, note in COMPARISONS[name]:
        got = np.asarray(output[key], dtype=float)
        want = np.asarray(native[f"{name}_{key}"], dtype=float)
        if got.shape != want.shape:
            raise SystemExit(f"{name}_{key}: shape {got.shape} != native {want.shape}")
        error = np.abs(got-want)
        ok = bool(np.allclose(got, want, rtol=rtol, atol=atol))
        case_ok &= ok
        outputs[key] = {"correctness_pass": ok, "rtol": rtol, "atol": atol, "tolerance_note": note,
                        "max_absolute_error": float(error.max()),
                        "max_relative_error": float((error/np.maximum(np.abs(want), max(atol, 1e-300))).max()),
                        "cantera_output": got.tolist(), "native_output": want.tolist()}
    native_samples = native[f"{name}_seconds"].tolist()
    ratio = statistics.median(samples)/statistics.median(native_samples) if samples and native_samples else None
    controlled = bool(controlled_base and len(samples) >= MIN_REPETITIONS
                      and len(native_samples) >= MIN_REPETITIONS)
    report["cases"][name] = {
        "correctness_pass": case_ok, "outputs": outputs, "timed_scope": TIMED_SCOPE[name],
        "cantera_first_seconds": first_seconds, "julia_first_seconds": float(native[f"{name}_first_seconds"][0]),
        "cantera_warm_seconds": samples, "julia_warm_seconds": native_samples,
        "warm_outputs_checked": len(samples), "native_warm_outputs_checked": len(native_samples),
        "speed_ratio": ratio, "performance_pass": bool(controlled and case_ok and ratio is not None and ratio >= .95),
        "qualification": "controlled" if controlled else "not_qualified",
    }
    any_failure |= not case_ok
    print(name, "correctness", case_ok, "speed ratio", ratio, "qualified", controlled, flush=True)
report["hardware"]["load_average_end"] = os.getloadavg()
args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
if any_failure:
    raise SystemExit("native and Cantera property example outputs differ")
