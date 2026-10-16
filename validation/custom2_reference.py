"""Original-source and refined custom2 histories, with no native-state inputs."""
from pathlib import Path
import argparse,ast,copy,hashlib,importlib.util,json,os,traceback
import numpy as np
import cantera as ct
import cantera._cantera as compiled

SOURCE_COMMIT='726522be4e2a13454d8415b7ef799d621f665cf3'
SOURCE_SHA='181db876f6af58844be3e3430045f7787478d8b4d9401ac33e2ef02e72f28d41'
MECHANISM_SHA='0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4'
SCALARS=('time','temperature','pressure','volume','mass','density','velocity',
    'internal_energy_mass','enthalpy_mass','entropy_mass','gibbs_mass','cp_mass','cv_mass',
    'mean_molecular_weight','gas_internal_energy','wall_kinetic_energy','mechanical_energy',
    'environment_pressure','pressure_work_rate')
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def source_program(source):
    tree=ast.parse(source.read_text())
    loops=[(i,n) for i,n in enumerate(tree.body) if isinstance(n,ast.While)]
    assert len(loops)==1
    index,loop=loops[0]
    assert ast.unparse(loop.test)=='net.time < 0.5'
    assert len(loop.body)==2 and ast.unparse(loop.body[0])=='net.advance(net.time + 0.005)'
    setup=ast.Module(body=tree.body[:index],type_ignores=[])
    observed=copy.deepcopy(loop)
    observed.body.append(ast.Expr(value=ast.Call(func=ast.Name(id='_capture_custom2',ctx=ast.Load()),args=[],keywords=[])))
    assert ast.dump(observed.body[0])==ast.dump(loop.body[0])
    assert ast.dump(observed.body[1])==ast.dump(loop.body[1])
    execution=ast.Module(body=[observed],type_ignores=[]);ast.fix_missing_locations(execution)
    return compile(setup,str(source),'exec'),compile(execution,str(source),'exec'),hashlib.sha256(ast.dump(loop).encode()).hexdigest()

def snapshot(ns):
    reactor=ns['r'];phase=reactor.phase;mass=reactor.mass;volume=reactor.volume;velocity=reactor.v_wall
    pressure=ns['res'].phase.P;energy=mass*phase.int_energy_mass;kinetic=50*velocity**2
    atoms=np.array([[phase.n_atoms(k,e) for k in range(phase.n_species)] for e in range(phase.n_elements)])
    return dict(time=ns['net'].time,temperature=phase.T,pressure=phase.P,volume=volume,mass=mass,
        density=phase.density,velocity=velocity,internal_energy_mass=phase.int_energy_mass,
        enthalpy_mass=phase.enthalpy_mass,entropy_mass=phase.entropy_mass,gibbs_mass=phase.gibbs_mass,
        cp_mass=phase.cp_mass,cv_mass=phase.cv_mass,mean_molecular_weight=phase.mean_molecular_weight,
        gas_internal_energy=energy,wall_kinetic_energy=kinetic,mechanical_energy=energy+kinetic+pressure*volume,
        environment_pressure=pressure,pressure_work_rate=phase.P*velocity,
        Y=phase.Y.copy(),X=phase.X.copy(),species_mass=mass*phase.Y,
        elements=atoms@(mass*phase.Y/phase.molecular_weights))

def arrays(rows,phase):
    data={key:np.asarray([row[key] for row in rows]) for key in SCALARS+('Y','X','species_mass','elements')}
    data['species_names_utf8']=np.frombuffer('\n'.join(phase.species_names).encode(),dtype=np.uint8)
    data['element_names_utf8']=np.frombuffer('\n'.join(phase.element_names).encode(),dtype=np.uint8)
    data['molecular_weights']=phase.molecular_weights
    return data

def main(args):
    source=args.cantera_source/'samples/python/reactors/custom2.py'
    mechanism=args.cantera_source/'data/h2o2.yaml'
    assert sha(source)==SOURCE_SHA and sha(mechanism)==MECHANISM_SHA and ct.__version__=='4.0.0a2'
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    inputs=output/'inputs';inputs.mkdir();(inputs/'h2o2.yaml').write_bytes(mechanism.read_bytes())
    exporter=load(args.source_root/'mechanism/export_sidecar.py','custom2_sidecar')
    exporter.export(inputs/'h2o2.yaml',inputs/'h2o2.yaml.npz')
    environment=load(args.source_root/'validation/benchmark_environment.py','custom2_environment')
    inventory=lambda:{str(p):sha(p) for p in [source,mechanism,Path(__file__),args.source_root/'mechanism/export_sidecar.py',args.source_root/'validation/benchmark_environment.py',*sorted(inputs.iterdir())]}
    libraries=lambda:{**environment.cantera_library_hashes(compiled.__file__),Path(compiled.__file__).name:sha(compiled.__file__)}
    before=inventory();threads=environment.verify_numerical_threads(set_accelerate=True);loaded=libraries()
    setup,loop,loop_sha=source_program(source)
    report=dict(passed=False,source_commit=SOURCE_COMMIT,source_sha256=SOURCE_SHA,source_loop_AST_sha256=loop_sha,
        original_loop_statements_unchanged=True,observation='One initial snapshot and one snapshot after each original SolutionArray append; no extra advances.',
        cantera_version=ct.__version__,host=environment.host_metadata(),inputs_before=before,
        numerical_threads_before=threads,libraries_before=loaded,cases={})
    prior=Path.cwd()
    try:
        os.chdir(inputs)
        for name,refined in (('original',False),('refined',True)):
            ns={'__name__':'__main__'};rows=[]
            exec(setup,ns)
            if refined:ns['net'].rtol=1e-11;ns['net'].atol=1e-20
            settings=dict(rtol=ns['net'].rtol,atol=ns['net'].atol)
            def capture():rows.append(snapshot(ns))
            ns['_capture_custom2']=capture;capture()
            try:exec(loop,ns)
            finally:
                data=arrays(rows,ns['r'].phase)
                # Preserve the example's actual plotted/stored arrays separately:
                # SolutionArray's TPY setter may renormalize at roundoff level.
                data.update(source_time=np.array(ns['states'].t),source_temperature=np.array(ns['states'].T),
                    source_volume=np.array(ns['states'].V),source_Y=np.array(ns['states'].Y))
                np.savez(output/(name+'.npz'),**data)
                report['cases'][name]=dict(**settings,points=len(rows),end_time=rows[-1]['time'],
                    archive_sha256=sha(output/(name+'.npz')))
            assert len(rows)==101 and rows[0]['time']==0 and .5<=rows[-1]['time']<.5+1e-14
            np.testing.assert_array_equal(data['time'],ns['states'].t)
            np.testing.assert_array_equal(data['temperature'],ns['states'].T)
            np.testing.assert_array_equal(data['volume'],data['source_volume'])
            array_Y_error=float(np.max(np.abs(data['Y']-data['source_Y'])))
            # This only checks duplicate serialization of the same state. The
            # independent trajectory Y gate remains the original 2e-6.
            report['cases'][name]['source_solution_array_Y_roundoff']=dict(max_absolute=array_Y_error,limit=1e-14)
            assert array_Y_error<=1e-14
            report['cases'][name]['source_solution_array_matches']=True
        report.update(inputs_after=inventory(),numerical_threads_after=environment.verify_numerical_threads(),libraries_after=libraries())
        assert report['inputs_after']==before and report['libraries_after']==loaded and report['numerical_threads_after']==threads
        report['passed']=True
    except Exception:
        report['failure']=traceback.format_exc();raise
    finally:
        os.chdir(prior)
        (output/'reference.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Original and refined pinned custom2 source histories saved.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--cantera-source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.source_root=args.source_root.resolve();args.cantera_source=args.cantera_source.resolve();main(args)
