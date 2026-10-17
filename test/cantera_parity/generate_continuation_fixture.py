"""Generate a fixed-grid Cantera two-point-control H2/O2 continuation fixture.

Requires Cantera 4.0.0a2 and writes a compressed NPZ fixture to argv[1]. The control update follows
samples/python/onedim/diffusion_flame_continuation.py at commit
726522be4e2a13454d8415b7ef799d621f665cf3, with refinement disabled so every
stored solution uses the same fixed grid.
"""
import json, os, sys
from pathlib import Path
import cantera as ct
import numpy as np

out=sys.argv[1]
mechanism_name=os.environ.get("CONTINUATION_MECHANISM","h2o2.yaml")
mechanism=Path(mechanism_name)
if not mechanism.is_absolute():
    mechanism=Path(__file__).resolve().parents[2]/"mechanism"/mechanism
fuel=os.environ.get("CONTINUATION_FUEL","H2:1")
oxidizer=os.environ.get("CONTINUATION_OXIDIZER","O2:1")
T_fuel=float(os.environ.get("CONTINUATION_T_FUEL","300"))
T_oxidizer=float(os.environ.get("CONTINUATION_T_OXIDIZER","500"))
pressure=float(os.environ.get("CONTINUATION_PRESSURE_PA","100000"))
mdot_fuel=float(os.environ.get("CONTINUATION_MDOT_FUEL","0.5"))
mdot_oxidizer=float(os.environ.get("CONTINUATION_MDOT_OXIDIZER","3.0"))
npoints=int(sys.argv[2]) if len(sys.argv)>2 else 12
maxsteps=int(sys.argv[3]) if len(sys.argv)>3 else 1000
f=ct.CounterflowDiffusionFlame(ct.Solution(str(mechanism)),width=.018)
f.P=pressure
f.fuel_inlet.mdot=mdot_fuel; f.fuel_inlet.X=fuel; f.fuel_inlet.T=T_fuel
f.oxidizer_inlet.mdot=mdot_oxidizer; f.oxidizer_inlet.X=oxidizer; f.oxidizer_inlet.T=T_oxidizer
f.flame.grid=np.linspace(0,.018,npoints)
f.solve(loglevel=0,refine_grid=False,auto=False)
initial_state={"grid_m":np.asarray(f.grid),"T_K":np.asarray(f.T),"Y":np.asarray(f.Y).T,
    "velocity_m_s":np.asarray(f.velocity),"spreadRate":np.asarray(f.flame.values("spreadRate")),
    "Lambda":np.asarray(f.flame.values("Lambda")),"fuel_mdot":float(f.fuel_inlet.mdot),
    "oxidizer_mdot":float(f.oxidizer_inlet.mdot)}
f.two_point_control_enabled=True
f.flame.set_bounds(spread_rate=(-1e-5,1e20))
f.max_time_step_count=100
initial_max=max(f.T); amax=strain=f.strain_rate("max")
step=20.; rows=[]; errors=0; pending_retries=0
for i in range(maxsteps):
    spacing=.6 if strain>.98*amax else .95
    target=float(np.min(f.T)+spacing*(np.max(f.T)-np.min(f.T)))
    backup=f.to_array()
    f.set_left_control_point(target); f.set_right_control_point(target)
    f.left_control_point_temperature-=step; f.right_control_point_temperature-=step
    if f.left_control_point_temperature<f.fuel_inlet.T+100 or f.right_control_point_temperature<f.oxidizer_inlet.T+100:
        break
    try:
        f.solve(loglevel=0,refine_grid=False,auto=False)
        errors=0
    except ct.CanteraError:
        f.from_array(backup); step*=.7; errors+=1; pending_retries+=1
        if strain/amax<.10 or errors>=3: break
        continue
    current=f.strain_rate("max")
    left_index=int(np.argmin(np.abs(np.asarray(f.grid)-f.left_control_point_coordinate)))
    right_index=int(np.argmin(np.abs(np.asarray(f.grid)-f.right_control_point_coordinate)))
    rows.append({"index":i,"Tmax_K":float(max(f.T)),"strain_rate_per_s":float(current),
        "control_temperature_K":target,"temperature_decrement_K":float(step),
        "retries_before_success":pending_retries,
        "left_control_point_coordinate_m":float(f.left_control_point_coordinate),
        "right_control_point_coordinate_m":float(f.right_control_point_coordinate),
        "left_control_point_index":left_index,"right_control_point_index":right_index,
        "left_control_point_temperature_K":float(f.T[left_index]),
        "right_control_point_temperature_K":float(f.T[right_index]),
        "grid_m":np.asarray(f.grid).tolist(),"T_K":np.asarray(f.T).tolist(),
        "Y":np.asarray(f.Y).T.tolist(),"velocity_m_s":np.asarray(f.velocity).tolist(),
        "fuel_mdot":float(f.fuel_inlet.mdot),"oxidizer_mdot":float(f.oxidizer_inlet.mdot),
        "spreadRate":np.asarray(f.flame.values("spreadRate")).tolist(),
        "Lambda":np.asarray(f.flame.values("Lambda")).tolist(),
        "solver_stats":f.solver_stats})
    if abs(max(f.T)-initial_max)<16: step=min(step+3.,100.)
    elif abs(max(f.T)-initial_max)>20: step*=18./abs(max(f.T)-initial_max)
    strain=current; amax=max(current,amax); initial_max=max(f.T)
turn=int(np.argmax([r["strain_rate_per_s"] for r in rows])) if rows else -1
turnrow=rows[turn] if turn>=0 else {}
np.savez_compressed(Path(out).with_suffix(".initial.npz"),
    grid=initial_state["grid_m"],T=initial_state["T_K"],Y=initial_state["Y"],
    velocity=initial_state["velocity_m_s"],spreadRate=initial_state["spreadRate"],
    Lambda=initial_state["Lambda"])
np.savez_compressed(out,grid_points=np.asarray([npoints]),pressure_Pa=np.asarray([pressure]),
    fuel_mdot=np.asarray([mdot_fuel]),oxidizer_mdot=np.asarray([mdot_oxidizer]),
    mechanism_utf8=np.frombuffer(mechanism.name.encode(),np.uint8),
    fuel_utf8=np.frombuffer(fuel.encode(),np.uint8),oxidizer_utf8=np.frombuffer(oxidizer.encode(),np.uint8),
    T_fuel_K=np.asarray([T_fuel]),T_oxidizer_K=np.asarray([T_oxidizer]),
    grid_m=np.asarray(turnrow.get("grid_m",f.grid)),
    Tmax_curve_K=np.asarray([r["Tmax_K"] for r in rows]),
    strain_rate_curve_per_s=np.asarray([r["strain_rate_per_s"] for r in rows]),
    turning_index=np.asarray([turn]),turning_T_K=np.asarray(turnrow.get("T_K",f.T)),
    turning_Y=np.asarray(turnrow.get("Y",f.Y.T)),
    turning_velocity_m_s=np.asarray(turnrow.get("velocity_m_s",f.velocity)),
    turning_Tmax_K=np.asarray([turnrow.get("Tmax_K",max(f.T))]),
    turning_strain_rate_per_s=np.asarray([turnrow.get("strain_rate_per_s",f.strain_rate("max"))]),
    turning_fuel_mdot=np.asarray([turnrow.get("fuel_mdot",f.fuel_inlet.mdot)]),
    turning_oxidizer_mdot=np.asarray([turnrow.get("oxidizer_mdot",f.oxidizer_inlet.mdot)]),
    turning_spreadRate=np.asarray(turnrow.get("spreadRate",f.flame.values("spreadRate"))),
    turning_Lambda=np.asarray(turnrow.get("Lambda",f.flame.values("Lambda"))),
    grid=np.asarray(turnrow.get("grid_m",f.grid)),
    initial_grid_m=initial_state["grid_m"], initial_T_K=initial_state["T_K"],
    initial_Y=initial_state["Y"], initial_velocity_m_s=initial_state["velocity_m_s"],
    initial_spreadRate=initial_state["spreadRate"], initial_Lambda=initial_state["Lambda"],
    initial_fuel_mdot=np.asarray([initial_state["fuel_mdot"]]),
    initial_oxidizer_mdot=np.asarray([initial_state["oxidizer_mdot"]]),
    curve_control_temperature_K=np.asarray([r["control_temperature_K"] for r in rows]),
    curve_temperature_decrement_K=np.asarray([r["temperature_decrement_K"] for r in rows]),
    curve_fuel_mdot=np.asarray([r["fuel_mdot"] for r in rows]),
    curve_oxidizer_mdot=np.asarray([r["oxidizer_mdot"] for r in rows]),
    curve_left_control_z_m=np.asarray([r["left_control_point_coordinate_m"] for r in rows]),
    curve_right_control_z_m=np.asarray([r["right_control_point_coordinate_m"] for r in rows]),
    curve_left_control_index=np.asarray([r["left_control_point_index"] for r in rows]),
    curve_right_control_index=np.asarray([r["right_control_point_index"] for r in rows]),
    curve_left_control_T_K=np.asarray([r["left_control_point_temperature_K"] for r in rows]),
    curve_right_control_T_K=np.asarray([r["right_control_point_temperature_K"] for r in rows]),
    curve_retry_count=np.asarray([r["retries_before_success"] for r in rows]),
    curve_T_K=np.asarray([r["T_K"] for r in rows]),
    curve_Y=np.asarray([r["Y"] for r in rows]),
    curve_velocity_m_s=np.asarray([r["velocity_m_s"] for r in rows]),
    curve_spreadRate=np.asarray([r["spreadRate"] for r in rows]),
    curve_Lambda=np.asarray([r["Lambda"] for r in rows]),
    T=np.asarray(turnrow.get("T_K",f.T)),
    Y=np.asarray(turnrow.get("Y",f.Y.T)),
    velocity=np.asarray(turnrow.get("velocity_m_s",f.velocity)),
    spreadRate=np.asarray(turnrow.get("spreadRate",f.flame.values("spreadRate"))),
    Lambda=np.asarray(turnrow.get("Lambda",f.flame.values("Lambda"))))
print(json.dumps({"points":len(rows),"turning_index":turn,
    "turning_strain":float(turnrow.get("strain_rate_per_s",float("nan"))),
    "turning_Tmax":float(turnrow.get("Tmax_K",float("nan"))),
    "last_strain":rows[-1]["strain_rate_per_s"] if rows else None,"grid_points":len(f.grid)}))
