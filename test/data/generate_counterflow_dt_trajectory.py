"""Generate the Cantera residual-ratio time-step decision trajectory fixture.

Runs the h2o2 counterflow diffusion flame of the timestep-controls example
(15 MPa, fixed 20-point grid, finite-difference Jacobian, residual-ratio
growth) through its transient fallback with a bounded step budget, and
records for every accepted pseudo-time step the solution before and after
the step, the step size, and the growth factor Cantera applied.

The callback reports dt after growth and before the tmax clip
(SteadyStateSystem::timeStep), so the applied factor is callback_dt divided
by the size of the accepted attempt. The attempt size is reconstructed at
full precision from the previous callback dt and the intervening failed
attempts (a failure multiplies dt by m_tfactor=0.5, except the fourth
successive failure which resets values without halving), starting from the
dt printed for the first attempt at loglevel 2; each reconstruction is
cross-checked against the dt printed for that attempt (4 significant
digits).

States are stored as Cantera component values in the layout
['T', species..., 'spreadRate', 'Lambda', 'velocity']; the native test
rescales them like generate_counterflow_growth_norm.py.
    python test/data/generate_counterflow_dt_trajectory.py
"""
from pathlib import Path
import io
import sys
import contextlib

import numpy as np
import cantera as ct

root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'test'/'cantera_parity'))

pressure=15e6
budget=300          # accepted steps recorded; Cantera throws when exhausted
gas=ct.Solution(str(root/'mechanism'/'h2o2.yaml'))
flame=ct.CounterflowDiffusionFlame(gas,grid=np.linspace(0.,30e-3,20))
flame.flame.jacobian_mode="finite-difference"
flame.max_time_step_count=budget
flame.set_refine_criteria(ratio=2.,slope=.06,curve=.08,prune=.02)
flame.time_step_regrid=0
flame.time_step_growth_strategy="residual-ratio"
flame.time_step_growth_factor=1.5
flame.P=pressure
flame.fuel_inlet.X="H2:1.0"; flame.fuel_inlet.T=800.
flame.oxidizer_inlet.X="O2:1.0"; flame.oxidizer_inlet.T=711.
# The example reads both inlet phases as one gas object, so the oxidizer
# mdot becomes 3*mdot_fuel despite different inlet states.
flame.fuel_inlet.mdot=.3
flame.oxidizer_inlet.mdot=.3/flame.fuel_inlet.phase.density* \
    flame.oxidizer_inlet.phase.density*3.
names=['T']+gas.species_names+['spreadRate','Lambda','velocity']

def snapshot(flame):
    return np.array([flame.flame.values(name) for name in names])

states=[snapshot(flame)]; callback_dts=[]
def record(dt):
    callback_dts.append(float(dt))
    states.append(snapshot(flame))
flame.set_time_step_callback(lambda dt: (record(dt),0)[1])

log=io.StringIO()
try:
    with contextlib.redirect_stdout(log):
        flame.solve(loglevel=3,auto=False,refine_grid=False)
    print("solve returned before exhausting the step budget")
except ct.CanteraError:
    pass  # budget exhausted: the trajectory up to it is the fixture

# SteadyStateSystem::solve passes loglevel-1 to timeStep, so loglevel=3 makes
# timeStep print its level-2 rows: "Timestep (n) dt= ... log(ss)= ..." starts
# an attempt and "... succeeded"/"... failed" reports its outcome (the failure
# row concatenates its debug suffix). Jacobian evaluations print exactly two
# messages in MultiNewton::solve ("Maximum Jacobian age reached" and
# "Re-evaluating Jacobian"); they are counted into the currently open attempt
# only, so evaluations of the between-retries steady solves are excluded.
lines=log.getvalue().splitlines()
Path('/tmp/dt_trajectory.log').write_text(log.getvalue())
attempts=[]
for line in lines:
    line=line.strip()
    if line.startswith("Timestep (") and "log(ss)=" in line:
        attempts.append(dict(dt=float(line.split("dt=")[1].split()[0]),
            ss=float(line.split("log(ss)=")[1]),jac=0,ok=None))
    elif "Maximum Jacobian age reached" in line or "Re-evaluating Jacobian" in line:
        if attempts and attempts[-1]["ok"] is None:
            attempts[-1]["jac"]+=1
    elif line.startswith("Timestep (") and "failed" in line:
        attempts[-1]["ok"]=False
    elif line.startswith("Timestep (") and "succeeded" in line:
        attempts[-1]["ok"]=True

# Reconstruct the accepted attempt dt at full precision alongside the
# printed 4-significant-digit value. A failure multiplies dt by m_tfactor=0.5
# while successiveFailures<=2; the third successive failure resets values
# and the counter instead (SteadyStateSystem::timeStep). resetBadValues may
# modify the state, which would desynchronize the recorded before-states;
# an unchanged log(ss) across a streak proves it did not.
reconstructed=[]; printed=[]; halvings=[]; jac_evals=[]; mismatches=0
cursor=None; streak=0; cuts=0; streak_ss=None; pending=False
for a in attempts:
    if cursor is None:
        cursor=a["dt"]  # seed: the first attempt size is only known as printed
    if not pending:
        streak_ss=a["ss"]
    elif abs(a["ss"]-streak_ss)>5e-5:
        raise SystemExit("state changed across a failure streak "
            "(resetBadValues modified x); before-states would desynchronize")
    pending=True
    if a["ok"] is False:
        streak+=1
        if streak>2:
            streak=0  # resetBadValues: no halving, counter restarts
        else:
            cursor*=0.5; cuts+=1
    elif a["ok"] is True:
        streak=0; pending=False
        reconstructed.append(cursor); printed.append(a["dt"])
        halvings.append(cuts); cuts=0; jac_evals.append(a["jac"])
        if len(reconstructed)<=len(callback_dts):
            # the next attempt continues from the dt the callback reported
            cursor=callback_dts[len(reconstructed)-1]
assert len(reconstructed)==len(callback_dts),(len(reconstructed),len(callback_dts),
    [l for l in lines if l.strip()][:8])
for k,(rec,pri) in enumerate(zip(reconstructed,printed)):
    if abs(rec/pri-1)>2e-4:
        mismatches+=1
        print(f"dt reconstruction mismatch at step {k+1}: {rec!r} vs printed {pri!r}")
        break
if mismatches:
    reconstructed=printed  # fall back to the printed values throughout

out=dict(grid=flame.grid,pressure=np.array([pressure]),mdot_fuel=np.array([.3]),
    mdot_oxidizer=np.array([flame.oxidizer_inlet.mdot]),budget=np.array([budget]),
    species_names_utf8=np.frombuffer('\n'.join(gas.species_names).encode(),dtype=np.uint8),
    cantera_version_utf8=np.frombuffer(ct.__version__.encode(),dtype=np.uint8))
out['values']=np.array(states)               # (steps+1, components, points)
out['attempt_dt']=np.array(reconstructed)    # dt of each accepted attempt
out['callback_dt']=np.array(callback_dts)    # grown dt reported per step
out['halvings']=np.array(halvings,dtype=np.int64)  # cuts before each accepted step
out['jac_evals']=np.array(jac_evals,dtype=np.int64)  # Jacobian evals per attempt
out['factor']=np.array([c/a for c,a in zip(callback_dts,reconstructed)])
out['reconstruction_ok']=np.array([mismatches==0])  # attempt_dt full precision
# With full-precision attempt sizes every applied factor is at least 1 after
# the seeded first step; a violation means an unmodeled dt change. Steps that
# evaluated a Jacobian cannot grow (the caller gates on m_jac->nEvals()).
assert mismatches or out['factor'][1:].min()>=1-1e-12
for k,jac in enumerate(jac_evals):
    assert not jac or abs(out['factor'][k]-1)<=(2e-3 if k==0 else 1e-9),k
np.savez(root/'test'/'data'/'cantera_counterflow_dt_trajectory.npz',**out)
grown=sum(1 for f in out['factor'] if f>1+1e-12)
print(f"accepted={len(callback_dts)} grown={grown} "
      f"ungated_grown={sum(1 for f,j in zip(out['factor'],jac_evals) if j==0 and f>1+1e-12)} "
      f"jacobian_steps={int((out['jac_evals']>0).sum())} "
      f"max_factor={out['factor'].max():.6f} "
      f"final_dt={callback_dts[-1]:.3e} Tmax={max(flame.T):.1f} "
      f"mismatches={mismatches}")
