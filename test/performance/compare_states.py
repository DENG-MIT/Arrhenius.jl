import glob, sys
import numpy as np

base_dir, candidate_dir = sys.argv[1:3]
fields = ("grid", "T", "Y", "velocity", "inlet_Y", "SL")
counts = ("stats_jacobians", "stats_newton_iterations", "stats_time_steps",
          "stats_regrids", "stats_final_points", "stats_failures")
paths = sorted(glob.glob(base_dir + "/*.npz"))
failures = 0
max_sum_y_error = 0.0
print("case,field,max_abs,max_relative,shape_equal")
for base_path in paths:
    name = base_path.rsplit("/", 1)[-1]
    candidate_path = candidate_dir + "/" + name
    with np.load(base_path, allow_pickle=False) as base, np.load(candidate_path, allow_pickle=False) as cand:
        if "Y" in cand:
            max_sum_y_error = max(
                max_sum_y_error,
                float(np.max(np.abs(np.sum(cand["Y"], axis=0) - 1.0), initial=0.0)),
            )
        for field in fields:
            a, b = np.asarray(base[field]), np.asarray(cand[field])
            same_shape = a.shape == b.shape
            if same_shape:
                delta = np.abs(b.astype(np.float64) - a.astype(np.float64))
                max_abs = float(np.max(delta, initial=0.0))
                scale = float(np.max(np.abs(a), initial=0.0))
                max_rel = max_abs / max(scale, 1e-300)
            else:
                max_abs = max_rel = float("inf")
            print(f"{name[:-4]},{field},{max_abs:.17g},{max_rel:.17g},{same_shape}")
            if not same_shape or max_rel > 1e-10:
                failures += 1
        for field in counts:
            if field not in base or field not in cand:
                continue
            a, b = np.asarray(base[field]), np.asarray(cand[field])
            equal = a.shape == b.shape and np.array_equal(a, b)
            print(f"{name[:-4]},{field},{a.tolist()},{b.tolist()},same={equal}")
            if not equal:
                print(f"{name[:-4]},{field},COUNT_MISMATCH,{a.tolist()},{b.tolist()}")
                failures += 1
        if ("batch_means" in base and "batch_means" in cand
                and np.size(base["batch_means"]) and np.size(cand["batch_means"])):
            base_time = float(np.mean(base["batch_means"]))
            candidate_time = float(np.mean(cand["batch_means"]))
            print(f"{name[:-4]},timing_seconds,{base_time:.9g},{candidate_time:.9g},ratio={candidate_time/base_time:.9g}")
print(f"compared_cases={len(paths)} failures={failures}")
print(f"max_abs_sum_Y_error={max_sum_y_error:.17g}")
sys.exit(bool(failures))
