#!/usr/bin/env python3
import json
import os
import sys
import numpy as np
import cantera as ct

def compare_case(spec_path, native_path, cantera_path):
    spec = json.load(open(spec_path, encoding="utf-8"))
    arr = np.load(native_path, allow_pickle=False)
    ref = np.load(cantera_path, allow_pickle=False)
    gas = ct.Solution(spec["mechanism"] + ".yaml")
    names = [gas.element_name(m) for m in range(gas.n_elements)]
    comp = np.array([[gas.n_atoms(k, m) * gas.atomic_weights[m] / gas.molecular_weights[k]
                      for m in range(gas.n_elements)] for k in range(gas.n_species)])

    def unpack_names(d):
        if "species_names_utf8" not in d.files:
            return list(gas.species_names)
        return bytes(d["species_names_utf8"].tolist()).decode("utf-8").split("\n")

    def errors(d):
        order = [unpack_names(d).index(name) for name in gas.species_names]
        if "inlet_Y" in d.files:
            inlet_y = d["inlet_Y"][order]
        else:
            gas.TPX = spec["T"], spec["P_atm"] * ct.one_atm, spec["X"]
            inlet_y = gas.Y.copy()
        inlet = comp.T @ inlet_y
        outlet = comp.T @ d["Y"][order, -1]
        return {name: abs(float(a - b)) for name, a, b in zip(names, inlet, outlet)}

    native_error, ct_error = errors(arr), errors(ref)
    allowed = {name: max(1.02 * ct_error[name], 1e-5) for name in names}
    inlet_temperature = max(float(spec.get("T", 300.0)), float(spec.get("T_surface", 0.0) or 0.0))
    burning_rise = float(np.max(ref["T"])) - inlet_temperature
    reference_burning = burning_rise >= 500.0
    counts = {}
    ct_count_keys = {"jacobians": "stats_jacobians", "time_steps": "stats_time_steps",
                     "regrids": "stats_regrids", "failures": "failures"}
    for key, ct_key in ct_count_keys.items():
        if key in arr.files and ct_key in ref.files:
            counts[key] = {"arrhenius": int(np.asarray(arr[key]).reshape(-1)[0]),
                           "cantera": int(np.asarray(ref[ct_key]).reshape(-1)[0])}
    def flux_profiles(data, is_native):
        data_names = unpack_names(data)
        order = [data_names.index(name) for name in gas.species_names]
        grid = np.asarray(data["grid"], dtype=float)
        y = np.asarray(data["Y"], dtype=float)[order]
        t = np.asarray(data["T"], dtype=float)
        nsp, npts = y.shape
        diff_face = np.zeros((nsp, max(npts - 1, 0)))
        if is_native:
            diff_face[:, :] = np.asarray(data["species_diffusive_flux"], dtype=float)[order, :npts-1]
        else:
            for j in range(npts - 1):
                ymid = 0.5 * (y[:, j] + y[:, j + 1])
                gas.TPY = 0.5 * (t[j] + t[j + 1]), spec["P_atm"] * ct.one_atm, ymid
                rho = gas.density
                mean_mw = gas.mean_molecular_weight
                x0 = y[:, j] / gas.molecular_weights
                x0 /= x0.sum()
                x1 = y[:, j + 1] / gas.molecular_weights
                x1 /= x1.sum()
                raw = gas.mix_diff_coeffs * gas.molecular_weights * rho / mean_mw * (x0 - x1) / (grid[j + 1] - grid[j])
                correction = -raw.sum()
                diff_face[:, j] = raw + correction * y[:, j]

        conv_node = float(spec["mdot"]) * y
        diff_node = np.empty((nsp, npts))
        diff_node[:, 0], diff_node[:, -1] = diff_face[:, 0], diff_face[:, -1]
        if npts > 2:
            diff_node[:, 1:-1] = 0.5 * (diff_face[:, :-1] + diff_face[:, 1:])
        total_node = conv_node + diff_node

        conv_element = comp.T @ conv_node
        diff_element = comp.T @ diff_node
        total_element = conv_element + diff_element
        return grid, {"convective": conv_element, "diffusive": diff_element, "total": total_element}

    z_native, flux_native = flux_profiles(arr, True)
    z_ct, flux_ct = flux_profiles(ref, False)
    z_union = np.unique(np.concatenate((z_native, z_ct)))
    left_cut = max(z_native[1], z_ct[1])
    right_cut = min(z_native[-2], z_ct[-2])
    regions = {"boundary": (z_union <= left_cut) | (z_union >= right_cut),
               "interior": (z_union > left_cut) & (z_union < right_cut)}
    flux_comparison = {}
    for field in ("convective", "diffusive", "total"):
        rows = {}
        for k, name in enumerate(names):
            an = np.interp(z_union, z_native, flux_native[field][k])
            ct_values = np.interp(z_union, z_ct, flux_ct[field][k])
            delta = an - ct_values
            i = int(np.argmax(np.abs(delta)))
            record = {"max_abs_delta": float(abs(delta[i])), "max_z": float(z_union[i]),
                "max_by_region": {region: float(np.max(np.abs(delta[mask]))) if np.any(mask) else 0.0
                    for region, mask in regions.items()}}
            if os.environ.get("ROUND7_FLUX_FULL", "0") == "1":
                record.update({"z": z_union.tolist(), "arrhenius": an.tolist(),
                    "cantera": ct_values.tolist(), "delta": delta.tolist()})
            rows[name] = record
        flux_comparison[field] = rows

    result = {"case": spec["id"], "reference_burning_rise_K": burning_rise,
              "reference_burning": reference_burning, "native_element_error": native_error,
              "cantera_element_error": ct_error, "gate_allowance": allowed,
              "element_pass": {name: native_error[name] <= allowed[name] for name in names},
              "counts": counts, "flux_comparison": flux_comparison}
    print(json.dumps(result, sort_keys=True))


args = sys.argv[1:]
if len(args) < 3 or len(args) % 3:
    raise SystemExit("usage: fixed_t_six_compare.py spec.json native.npz cantera.npz [spec.json native.npz cantera.npz ...]")
for i in range(0, len(args), 3):
    compare_case(*args[i:i + 3])
