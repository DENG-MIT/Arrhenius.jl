"""Compare complete custom2 source/default/refined/native histories and CLI output."""
from pathlib import Path
import argparse,hashlib,json,tomllib,traceback
import numpy as np

sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
decode=lambda data,key:bytes(data[key]).decode().splitlines()
PRIMARY_LIMITS=dict(temperature=.05,pressure_relative=2e-5,volume=2e-6,Y=2e-6,velocity=2e-3,time=1e-14)
THERMO_KEYS=('density','internal_energy_mass','enthalpy_mass','entropy_mass','gibbs_mass','cp_mass','cv_mass',
    'mean_molecular_weight','gas_internal_energy','wall_kinetic_energy','mechanical_energy','pressure_work_rate')
THERMO_LIMIT=2e-5
PRIOR_PUBLIC_ENGINE_SHA='5301ebcb2f44beb2b773efbf1ea5c8d7860c1e5e0bd0c3e553d46fba98a9dec2'
PINNED_LAZY_LIBRARIES={
    'libmkl_avx2.so.2':'0696e8d6caa451052c1898ef5bbfdbccf6f3c2548cf049e55f9db213979f1e17',
    'libmkl_vml_avx2.so.2':'91c94d26b346837422a703c3b6183f88f6c9a96b545cf08158258a7926851180'}

def pinned_library_receipt(path):
    assert sha(path)==PRIOR_PUBLIC_ENGINE_SHA,'prior published engine receipt changed'
    recorded=json.loads(path.read_text())['public_reproduction']['loaded_numerical_library_sha256']
    assert all(recorded.get(name)==digest for name,digest in PINNED_LAZY_LIBRARIES.items())
    return dict(receipt_sha256=PRIOR_PUBLIC_ENGINE_SHA,pinned_lazy_additions=PINNED_LAZY_LIBRARIES)

def numerical_library_changes(before,after,allowed=None):
    pins={} if allowed is None else allowed
    removed={name:digest for name,digest in before.items() if name not in after}
    changed={name:[digest,after[name]] for name,digest in before.items() if name in after and after[name]!=digest}
    added={name:digest for name,digest in after.items() if name not in before}
    unlisted={name:digest for name,digest in added.items() if name not in pins}
    incorrectly_pinned={name:digest for name,digest in added.items()
        if name in pins and digest!=pins[name]}
    return dict(passed=not (removed or changed or unlisted or incorrectly_pinned),removed=removed,
        changed=changed,added=added,unlisted_additions=unlisted,incorrectly_pinned_additions=incorrectly_pinned)

def identity(data):
    times=data['time'];n=len(times)
    assert n==101 and times.shape==(n,) and times[0]==0 and abs(times[-1]-.5)<=1e-14
    assert np.all(np.diff(times)>0) and len(decode(data,'species_names_utf8'))==10
    assert data['Y'].shape==data['X'].shape==data['species_mass'].shape==(n,10)
    assert data['elements'].shape==(n,len(decode(data,'element_names_utf8')))
    assert data['molecular_weights'].shape==(10,) and np.all(data['molecular_weights']>0)
    assert np.isfinite(data['Y']).all() and np.isfinite(data['X']).all()
    for key in THERMO_KEYS+('temperature','pressure','volume','mass','velocity','environment_pressure'):
        assert data[key].shape==(n,) and np.isfinite(data[key]).all()
    assert np.all(data['mass']>0) and np.all(data['volume']>0) and np.all(data['temperature']>0)

def comparison(left,right):
    identity(left);identity(right)
    assert decode(left,'species_names_utf8')==decode(right,'species_names_utf8')
    np.testing.assert_array_equal(left['molecular_weights'],right['molecular_weights'])
    errors=dict(time=float(np.max(np.abs(left['time']-right['time']))),
        temperature=float(np.max(np.abs(left['temperature']-right['temperature']))),
        pressure_relative=float(np.max(np.abs(left['pressure']/right['pressure']-1))),
        volume=float(np.max(np.abs(left['volume']-right['volume']))),
        Y=float(np.max(np.abs(left['Y']-right['Y']))),velocity=float(np.max(np.abs(left['velocity']-right['velocity']))))
    primary={key:dict(value=value,limit=PRIMARY_LIMITS[key],passed=value<PRIMARY_LIMITS[key]) for key,value in errors.items()}
    thermo={key:dict(value=float(np.max(np.abs(left[key]-right[key]))),
        scale=max(float(np.max(np.abs(right[key]))),1.),relative_limit=THERMO_LIMIT) for key in THERMO_KEYS}
    for row in thermo.values():row['passed']=row['value']<THERMO_LIMIT*row['scale']
    element_names=decode(left,'element_names_utf8');other=decode(right,'element_names_utf8')
    assert set(element_names)==set(other)
    mapped=right['elements'][:,[other.index(e) for e in element_names]]
    extra=dict(mass_relative=float(np.max(np.abs(left['mass']/right['mass']-1))),
        X=float(np.max(np.abs(left['X']-right['X']))),
        species_mass_relative_to_total=float(np.max(np.abs(left['species_mass']-right['species_mass'])/right['mass'][:,None])),
        elements_kmol=float(np.max(np.abs(left['elements']-mapped))))
    extra_limits=dict(mass_relative=1e-10,X=2e-6,species_mass_relative_to_total=2e-6,elements_kmol=1e-11)
    derived={key:dict(value=value,limit=extra_limits[key],passed=value<extra_limits[key]) for key,value in extra.items()}
    return dict(passed=all(r['passed'] for group in (primary,thermo,derived) for r in group.values()),primary=primary,
        additional_thermodynamic_outputs=thermo,additional_state_outputs=derived)

def conservation(data,native=False):
    scale=max(abs(float(data['gas_internal_energy'][0])),1.)
    mass=float(np.max(np.abs(data['mass']/data['mass'][0]-1)))
    elements=float(np.max(np.abs(data['elements']-data['elements'][0])))
    mechanical=float(np.max(np.abs(data['mechanical_energy']-data['mechanical_energy'][0]))/scale)
    checks=dict(mass_drift=dict(value=mass,limit=1e-10,passed=mass<1e-10),
        element_drift_kmol=dict(value=elements,limit=1e-11,passed=elements<1e-11),
        mechanical_energy_relative_drift=dict(value=mechanical,limit=2e-7,passed=mechanical<2e-7))
    if native:
        assert data['state_full'].shape==(16,101) and np.isfinite(data['state_full']).all()
        rawY=data['state_full'][:10,:]/np.sum(data['state_full'][:10,:],axis=0)
        energy=float(np.max(np.abs(data['energy_balance']-data['energy_balance'][0]))/scale)
        ledger_mass=float(np.max(np.abs(data['mass_balance']/data['mass_balance'][0]-1)))
        checks.update(energy_balance_relative_drift=dict(value=energy,limit=2e-7,passed=energy<2e-7),
            mass_ledger_relative_drift=dict(value=ledger_mass,limit=1e-10,passed=ledger_mass<1e-10),
            raw_species_domain=dict(minimum_Y=float(np.min(rawY)),minimum=-1e-13,passed=bool(np.min(rawY)>=-1e-13)),
            normalized_mass_fractions=dict(error=float(np.max(np.abs(np.sum(data['Y'],axis=1)-1))),limit=1e-12,passed=bool(np.max(np.abs(np.sum(data['Y'],axis=1)-1))<1e-12)))
    return dict(passed=all(row['passed'] for row in checks.values()),checks=checks)

def main(args):
    n=np.load(args.native/'native.npz',allow_pickle=False)
    original=np.load(args.reference/'original.npz',allow_pickle=False)
    refined=np.load(args.reference/'refined.npz',allow_pickle=False)
    native=tomllib.loads((args.native/'native.toml').read_text())
    reference=json.loads((args.reference/'reference.json').read_text())
    binding=(dict(mode='pinned_lazy_additions',**pinned_library_receipt(args.lazy_library_receipt))
        if args.lazy_library_receipt is not None else dict(mode='strict',pinned_lazy_additions={}))
    libraries=numerical_library_changes(native['runtime_before']['numerical_library_sha256'],
        native['runtime_after']['numerical_library_sha256'],binding['pinned_lazy_additions'])
    checks=dict(native_completed=native['completed'],reference_completed=reference['passed'],
        native_inputs_frozen=native['inputs_unchanged'] and native['inputs_before']==native['inputs_after'],
        native_numerical_libraries_validated=libraries['passed'],
        reference_inputs_frozen=reference['inputs_before']==reference['inputs_after'],
        reference_libraries_frozen=reference['libraries_before']==reference['libraries_after'],
        native_archive_binding=sha(args.native/'native.npz')==native['native_sha256'],
        reference_archives_bound=all(sha(args.reference/(name+'.npz'))==reference['cases'][name]['archive_sha256'] for name in ('original','refined')),
        actual_cli_csv_bitwise_match=(args.native/'native.csv').read_bytes()==args.cli.read_bytes())
    checks['native_actual_threads']=all(all(v==1 for v in native[key]['settings'].values()) for key in ('runtime_before','runtime_after'))
    checks['reference_actual_threads']=reference['numerical_threads_before']==reference['numerical_threads_after'] and all(x['threads']==1 for x in reference['numerical_threads_after']['threadpools']) and reference['numerical_threads_after'].get('accelerate_threading_mode',1)==1
    for name in ('h2o2.yaml','h2o2.yaml.npz'):
        digest=sha(args.reference/'inputs'/name)
        matches=[v for p,v in native['inputs_before'].items() if Path(p).name==name]
        checks['input_binding_'+name]=matches==[digest]
    comparisons=dict(reference_refinement=comparison(original,refined),native_original=comparison(n,original),native_refined=comparison(n,refined))
    balances=dict(native=conservation(n,native=True),original=conservation(original),refined=conservation(refined))
    checks['all_trajectory_outputs']=all(row['passed'] for row in comparisons.values())
    checks['all_conservation_gates']=all(row['passed'] for row in balances.values())
    report=dict(passed=all(checks.values()),checks=checks,comparisons=comparisons,conservation=balances,
        numerical_library_validation=dict(**binding,**libraries),
        time_mapping='One-to-one 101-point grid comparison with the unchanged 1e-14 s limit; the exact source repeated-addition endpoint may exceed 0.5 s by roundoff.',
        additional_thermo_gate='Absolute maximum error divided by max(full-history reference absolute peak,1 SI unit) < 2e-5. Existing T/P/V/Y/velocity/time and conservation limits unchanged.',
        artifacts={str(path):sha(path) for path in [args.native/'native.npz',args.native/'native.toml',args.native/'native.csv',
            args.reference/'original.npz',args.reference/'refined.npz',args.reference/'reference.json',args.cli,Path(__file__)]},
        performance_qualified=False)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    if not report['passed']:raise ValueError('custom2 complete comparison failed; report saved')
    print('Complete custom2 source/native/CLI comparison passed.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('native','reference','cli','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--lazy-library-receipt',type=Path,
        help='Optional independently pinned public WSL engine report; default rejects every newly loaded library.')
    args=parser.parse_args()
    try:main(args)
    except Exception:
        if not args.output.exists():
            args.output.write_text(json.dumps(dict(passed=False,failure=traceback.format_exc()),indent=2)+'\n')
        raise
