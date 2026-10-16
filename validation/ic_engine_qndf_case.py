"""Validate the full public QNDF engine against the pinned Cantera 4 source.

python validation/ic_engine_qndf_case.py prepare --source-mechanism CANTERA/data/nDodecane_Reitz.yaml --output engine-reference
julia --project=CALLER_ENV validation/ic_engine_qndf_case.jl engine-reference/dodecane_IG.yaml engine-reference/native.npz
python validation/ic_engine_qndf_case.py compare --native engine-reference/native.npz --native-mechanism engine-reference/dodecane_IG.yaml --source-example CANTERA/samples/python/reactors/ic_engine.py --source-mechanism CANTERA/data/nDodecane_Reitz.yaml --output engine-reference/comparison.json

Use Cantera commit 726522be4e2a13454d8415b7ef799d621f665cf3 and one numerical
thread. Preparation exports numerical data only. Comparison starts every CT
trajectory from the original source state; native states are never CT inputs.
This is a correctness driver, with no timing or performance qualification.
"""
from pathlib import Path
from itertools import islice
import argparse
import hashlib
import importlib.util
import inspect
import json
import tomllib
import numpy as np
import cantera as ct
import cantera._cantera as compiled
from ic_engine_source import (source_programs,new_source_state,source_snapshot,reference_segments,
    integral_terms,compare_native,verify_mechanisms,sha,SOURCE_COMMIT,SOURCE_SHA256,
    REFINED_RTOL,REFINED_ATOL)
from benchmark_environment import cantera_library_hashes,verify_numerical_threads

SOURCE_MECHANISM_SHA256="3d3b59ed91dec0d0bcbac2fa2ef2cba13fbd565bfff8f266ca847ed6aa92f7f1"


class NativeView(dict):
    @property
    def files(self):
        return self.keys()


def prepare(source,output):
    output.mkdir(parents=True,exist_ok=True)
    gas=ct.Solution(str(source),"nDodecane_IG")
    if gas.n_species!=100 or gas.n_reactions!=553:
        raise ValueError("the source nDodecane_IG phase requires 100 species and 553 reactions")
    mechanism=output/"dodecane_IG.yaml"
    gas.write_yaml(str(mechanism))
    path=Path(__file__).resolve().parents[1]/"mechanism"/"export_sidecar.py"
    spec=importlib.util.spec_from_file_location("engine_sidecar",path)
    exporter=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    exporter.export(mechanism,Path(str(mechanism)+".npz"))
    (output/"preparation.json").write_text(json.dumps(dict(cantera_version=ct.__version__,
        source_mechanism_sha256=sha(source),mechanism_sha256=sha(mechanism),
        sidecar_sha256=sha(str(mechanism)+".npz")),indent=2)+"\n")


def complete_histories(native,reference):
    names=("time","temperature","pressure","volume","mass","Y","entropy_mass",
        "mdot_in","mdot_out","mdot_fuel","work_rate","heat_release_rate")
    actual=NativeView({"output_"+name:np.concatenate([
        native[f"segment_{i:02d}_output_{name}"] for i in range(1,26)]) for name in names})
    actual["output_time"]=np.concatenate([native[f"segment_{i:02d}_time"] for i in range(1,26)])
    expected={name:np.concatenate([reference[f"segment_{i:02d}_{name}"]
        for i in range(1,26)]) for name in names}
    return actual,expected


def require_native(native,meta):
    if not meta["checks_pass"] or not meta["checks"]["checks_pass"]:
        raise ValueError("native calculation did not pass its physical gates")
    assert meta["inputs_unchanged"] and meta["inputs_before"]==meta["inputs_after"]
    for name in ("mechanism_sha256","sidecar_sha256"):
        assert meta[name]==meta["inputs_before"][name]
    for name in ("runtime_before","runtime_after"):
        assert meta[name]["checks_pass"] and all(value==1 for value in meta[name]["settings"].values())
    species=bytes(native["species_names_utf8"]).decode().splitlines()
    assert len(species)==len(set(species))==100
    initial=native["exact_source_initial_state"]
    assert initial.shape==(105,) and np.isfinite(initial).all()
    assert native["quadrature_terms"].shape==native["coarse_quadrature_terms"].shape==(4,)
    assert native["ledger_work"].shape==native["checks_pass"].shape==(1,)
    assert native["checks_pass"][0]
    assert native["states"].ndim==2 and native["states"].shape[0]==105
    s=meta["summary"];c=meta["checks"];records=meta["segments"]
    assert s["end_time_s"]==.16 and len(records)==25
    assert [r["integrated_coordinate_count"] for r in records]==[8]*3+[105]*22
    assert all(r["rms_norm_denominator"]==r["full_state_count"]==105 for r in records)
    assert c["mass_history_relative_drift"]<2e-7 and c["energy_history_relative_drift"]<2e-6
    assert max(c["quadrature_relative_change"])<1e-4
    assert c["efficiency_quadrature_change"]<1e-4 and c["CO_quadrature_change"]<1e-4
    assert c["work_quadrature_ledger_relative_error"]<1e-4 and c["prescribed_fuel_source_pass"]
    assert s["minimum_mass_fraction"]>=-1e-13 and s["volume_identity_error_m3"]<1e-10
    assert s["maximum_output_interval_s"]<=1/(360*50)*(1+1e-12)
    assert s["maximum_output_temperature_change_K"]<=20+1e-8
    for i,r in enumerate(records,1):
        assert r["element_flux_relative_error"]<2e-7 and r["element_flux_quadrature_change"]<2e-7
        assert r["maximum_reference_to_mass_ratio"]<=1
        state=native[f"segment_{i:02d}_state"]
        times=native[f"segment_{i:02d}_time"]
        assert times.ndim==1 and len(times)>=2 and np.all(np.diff(times)>0)
        assert state.shape==(105,len(times)) and np.isfinite(state).all()
        for name in ("time","temperature","pressure","volume","mass","entropy_mass",
                "mdot_in","mdot_out","mdot_fuel","work_rate","heat_release_rate"):
            assert native[f"segment_{i:02d}_output_{name}"].shape==times.shape
        assert native[f"segment_{i:02d}_output_Y"].shape==(len(times),100)
        if i>1:
            assert np.array_equal(state[:,0].copy().view(np.uint64),
                native[f"segment_{i-1:02d}_state"][:,-1].copy().view(np.uint64))


def pointwise_reference(programs,native,save_path):
    ns=new_source_state(programs);network=ns["sim"]
    assert ns["rtol"]==1e-12 and ns["atol"]==1e-16
    assert bytes(native["species_names_utf8"]).decode().splitlines()==ns["cyl"].phase.species_names
    initial=native["exact_source_initial_state"]
    np.testing.assert_allclose(initial[:100]/sum(initial[:100]),ns["cyl"].phase.Y,rtol=2e-14,atol=0.)
    assert initial[100]==ns["cyl"].T==300.
    saved={}
    try:
        for index,(start,stop) in enumerate(reference_segments(ns),1):
            key=f"segment_{index:02d}";times=native[key+"_time"]
            assert abs(times[0]-start)<=2*np.spacing(max(start,1e-300))
            assert abs(times[-1]-stop)<=2*np.spacing(stop) and np.all(np.diff(times)>0)
            rows=[];ys=[]
            for t in times:
                local=float(t-start)
                if local>network.time:
                    network.advance(local,apply_limit=False)
                rows.append(source_snapshot(ns,start));ys.append(ns["cyl"].phase.Y.copy())
            reference={name:np.array([row[name] for row in rows]) for name in rows[0]}
            reference["Y"]=np.asarray(ys)
            saved.update({key+"_"+name:value for name,value in reference.items()})
    finally:
        np.savez(save_path,**saved)
    return saved


def reference_integrals(programs,trace_path):
    """Accepted CT states plus endpoints from a second, source-only CT instance."""
    ns=new_source_state(programs)
    endpoint_ns=new_source_state(programs)
    network=ns["sim"]
    endpoint_network=endpoint_ns["sim"]
    total=np.zeros(4);coarse=np.zeros(4);saved={};counts=[];step_records=[]
    # This cap refines only the independent integral oracle, after the separate
    # source-trajectory comparison. It does not change native solver settings.
    quadrature_step_cap=1/(4*360*50)
    intervals=zip(islice(reference_segments(ns),25),islice(reference_segments(endpoint_ns),25))
    for index,((start,stop),endpoint_interval) in enumerate(intervals,1):
        assert endpoint_interval==(start,stop)
        # Both instances start from the exact source state, then from the same
        # independently propagated CT endpoint. No native state is an input.
        entry_state=ns["cyl"].phase.state.copy()
        endpoint_entry=endpoint_ns["cyl"].phase.state.copy()
        saved[f"segment_{index:02d}_step_entry_phase_state"]=entry_state
        saved[f"segment_{index:02d}_endpoint_entry_phase_state"]=endpoint_entry
        saved[f"segment_{index:02d}_entry_volumes"]=np.array([ns["cyl"].volume,endpoint_ns["cyl"].volume])
        np.savez(trace_path,**saved)
        assert np.array_equal(entry_state.view(np.uint64),endpoint_entry.view(np.uint64))
        assert ns["cyl"].volume==endpoint_ns["cyl"].volume
        rows=[source_snapshot(ns,start)]
        local_times=[network.time]
        duration=stop-start
        # This separate instance has not crossed the event: advance() moves
        # forward from local time zero and returns its exact dense endpoint.
        assert endpoint_network.time==0.
        endpoint_network.advance(duration,apply_limit=False)
        assert endpoint_network.time==duration
        endpoint_row=source_snapshot(endpoint_ns,start)
        endpoint_state=endpoint_ns["cyl"].phase.state.copy()
        saved[f"segment_{index:02d}_endpoint_phase_state"]=endpoint_state
        saved[f"segment_{index:02d}_endpoint_reactor_state"]=endpoint_ns["cyl"].get_state()
        network.max_time_step=quadrature_step_cap
        try:
            while network.time<duration:
                previous=network.time
                network.step()
                assert network.time>previous
                if network.time>=duration:
                    saved[f"segment_{index:02d}_excluded_crossing_time"]=np.array([network.time])
                    saved[f"segment_{index:02d}_excluded_crossing_reactor_state"]=ns["cyl"].get_state()
                    break
                rows.append(source_snapshot(ns,start))
                local_times.append(network.time)
            rows.append(endpoint_row)
            local_times.append(duration)
        finally:
            output={name:np.array([row[name] for row in rows]) for name in rows[0]}
            output["local_time"]=np.asarray(local_times)
            saved.update({f"segment_{index:02d}_integral_{name}":value for name,value in output.items()})
            np.savez(trace_path,**saved)
        assert local_times[-1]==duration and np.all(np.diff(local_times)>0)
        assert max(np.diff(local_times))<=quadrature_step_cap*(1+1e-10)
        quadrature_output={**output,"time":output["local_time"]}
        total+=integral_terms(quadrature_output)
        selected=np.unique(np.r_[np.arange(0,len(rows),2),len(rows)-1])
        coarse+=integral_terms(quadrature_output,selected)
        counts.append(len(rows))
        step_records.append(dict(max_step_cap_s=quadrature_step_cap,
            pre_event_accepted_observations=len(rows)-2,
            endpoint_origin="separate CT instance advanced forward from identical CT interval entry",
            entry_phase_and_volume_bitwise_equal=True,
            excluded_crossing_local_time_s=network.time,
            final_panel_span_s=local_times[-1]-local_times[-2],
            final_local_time_s=local_times[-1],interval_duration_s=duration,
            maximum_observed_step_s=float(np.max(np.diff(local_times))),
            solver_stats=network.solver_stats,endpoint_solver_stats=endpoint_network.solver_stats))
        # Cantera4 phase.state invokes ThermoPhase::restoreState directly;
        # volume is a separate exact scalar. reference_segments next resets
        # local time and initializes the next frozen source equation branch.
        ns["cyl"].phase.state=endpoint_state
        ns["cyl"].volume=endpoint_ns["cyl"].volume
        assert np.array_equal(ns["cyl"].phase.state.view(np.uint64),endpoint_state.view(np.uint64))
    return total,coarse,saved,counts,step_records


def compare(args,report):
    native=np.load(args.native,allow_pickle=False)
    meta=tomllib.loads(args.native.with_suffix(".toml").read_text())
    assert meta["native_sha256"]==sha(args.native)
    require_native(native,meta)
    provenance=verify_mechanisms(args.source_mechanism,args.native_mechanism,meta)
    report["mechanism_provenance"]={key:value for key,value in provenance.items() if key.endswith("sha256")}
    programs=source_programs(args.source_example,args.source_mechanism)
    reference=pointwise_reference(programs,native,args.output.with_suffix(".reference.npz"))
    actual,expected=complete_histories(native,reference)
    # Check the complete-history trajectory and rate limits before the separate
    # quadrature solve. Equal integral inputs isolate these existing gates.
    pointwise=compare_native(actual,meta["summary"],expected,meta["summary"]["integrals"])
    report["pointwise_comparison"]=pointwise
    if not pointwise["correctness_pass"]:
        raise ValueError("engine trajectory failed original complete-source gates")
    terms,coarse,_,counts,steps=reference_integrals(programs,args.output.with_suffix(".integrals.npz"))
    changes=np.abs(coarse/terms-1)
    efficiency_change=abs((coarse[1]/coarse[0])/(terms[1]/terms[0])-1)
    co_change=abs((coarse[2]/coarse[3])/(terms[2]/terms[3])-1)
    converged=bool(np.all(terms>0) and np.all(coarse>0) and max(changes)<1e-4
        and efficiency_change<1e-4 and co_change<1e-4)
    report["reference_convergence"]=dict(checks_pass=converged,relative_limit=1e-4,
        terms=terms.tolist(),coarse_terms=coarse.tolist(),term_changes=changes.tolist(),
        efficiency_change=float(efficiency_change),CO_change=float(co_change),points=counts,
        steps=steps,rtol=REFINED_RTOL,atol=REFINED_ATOL)
    if not converged:
        raise ValueError("independent reference integrals did not converge; comparison skipped")
    values=dict(heat_J=float(terms[0]),work_J=float(terms[1]),
        efficiency=float(terms[1]/terms[0]),CO_ppm=float(1e6*terms[2]/terms[3]))
    report["comparison"]=compare_native(actual,meta["summary"],expected,values)
    report["checks_pass"]=report["comparison"]["correctness_pass"]
    if not report["checks_pass"]:
        raise ValueError("engine failed original full-source comparison gates")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest="command",required=True)
    prep=commands.add_parser("prepare")
    prep.add_argument("--source-mechanism",type=Path,required=True)
    prep.add_argument("--output",type=Path,required=True)
    check=commands.add_parser("compare")
    for name in ("native","native-mechanism","source-example","source-mechanism","output"):
        check.add_argument("--"+name,type=Path,required=True)
    args=parser.parse_args()
    if not ct.__version__.startswith("4.0."):
        raise RuntimeError("Cantera 4.0 is required")
    if sha(args.source_mechanism)!=SOURCE_MECHANISM_SHA256:
        raise ValueError("nDodecane_Reitz.yaml differs from the pinned source mechanism")
    if args.command=="prepare":
        prepare(args.source_mechanism,args.output)
        return
    args.output.parent.mkdir(parents=True,exist_ok=True)
    report=dict(checks_pass=False,performance_qualified=False,cantera_version=ct.__version__,
        source_commit=SOURCE_COMMIT,source_sha256=SOURCE_SHA256,
        source_mechanism_sha256=sha(args.source_mechanism),native_sha256=sha(args.native),
        native_record_sha256=sha(args.native.with_suffix(".toml")),
        driver_sha256=sha(__file__),comparison_file_sha256=sha(inspect.getsourcefile(compare_native)),
        comparison_function_sha256=hashlib.sha256(inspect.getsource(compare_native).encode()).hexdigest(),
        cantera_library_sha256={**cantera_library_hashes(compiled.__file__),Path(compiled.__file__).name:sha(compiled.__file__)},
        thread_checks={"before":verify_numerical_threads(set_accelerate=True)})
    try:
        compare(args,report)
        report["thread_checks"]["after"]=verify_numerical_threads()
    except BaseException as error:
        report["checks_pass"]=False
        report["error"]=str(error)
        raise
    finally:
        args.output.write_text(json.dumps(report,indent=2)+"\n")
    print("Complete eight-revolution engine comparison passed.")


if __name__=="__main__":
    main()
