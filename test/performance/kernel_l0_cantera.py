import json, statistics, sys, time
import cantera as ct
import numpy as np

spec = json.load(open(sys.argv[1]))
gas = ct.Solution(spec["mechanism"] + ".yaml")
X = {k: float(v) for k, v in (item.split(":") for item in spec["X"].split(","))}
T = 1500.0
P = spec["P_atm"] * ct.one_atm
if len(sys.argv) > 2 and sys.argv[2].endswith(".npz"):
    with np.load(sys.argv[2], allow_pickle=False) as state:
        state_index = max(0, (state["T"].size - 1) // 2)
        T = float(state["T"][state_index])
        Y = np.asarray(state["Y"][:, state_index], dtype=np.float64)
    gas.TPY = T, P, Y
else:
    state_index = -1
    T = float(sys.argv[2]) if len(sys.argv) > 2 else T
    gas.TPX = T, P, X
Y_base = gas.Y.copy()

def measure(fn, n=3000, reps=7):
    fn()
    times = []
    for _ in range(reps):
        t0 = time.perf_counter_ns()
        for i in range(n):
            fn()
        times.append((time.perf_counter_ns() - t0) / n)
    return round(statistics.median(times)), [round(v) for v in times]

toggle = [0]
def set_state():
    toggle[0] ^= 1
    composition = Y_base.copy()
    o2 = gas.species_index("O2")
    n2 = gas.species_index("N2")
    delta = (1.0 if toggle[0] else -1.0) * 1e-10
    composition[o2] += delta
    composition[n2] -= delta
    gas.TPY = T, P, composition

setter, setter_samples = measure(set_state)
def measured(prop):
    def call():
        set_state()
        getattr(gas, prop)
    total, samples = measure(call)
    return max(0, total - setter), total, samples

print(f"Cantera={ct.__version__} mechanism={spec['mechanism']} nspecies={gas.n_species} nreactions={gas.n_reactions} T={T} P={P} state_index={state_index + 1}")
print(f"state_setter_ns={setter} samples_ns={setter_samples}")
for label, prop in (("cp", "partial_molar_cp"), ("h", "partial_molar_enthalpies"),
                    ("s", "partial_molar_entropies"), ("forward_rates", "forward_rate_constants"),
                    ("reverse_rates", "reverse_rate_constants"), ("equilibrium", "equilibrium_constants"),
                    ("net_production", "net_production_rates")):
    delta, total, samples = measured(prop)
    print(f"{label}_adjusted_ns={delta} total_ns={total} samples_ns={samples}")
