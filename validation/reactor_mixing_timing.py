"""Time the complete pinned mixer on prepared models, checking every output."""
import argparse,ast,hashlib,importlib.util,json,statistics,time
from pathlib import Path
import cantera as ct
import cantera._cantera as compiled
import numpy as np

SOURCE_SHA256='87ccbfa629cb6cd494afae82d206c60b215261e16982675264db017be4314e11'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def snapshot(ns):
    p=ns['mixer'].phase
    return dict(Y=p.Y.copy(),X=p.X.copy(),temperature=np.array([p.T]),pressure=np.array([p.P]),
        mass=np.array([ns['mixer'].mass]),density=np.array([p.density]),enthalpy=np.array([p.enthalpy_mass]),
        internal_energy=np.array([p.int_energy_mass]),entropy=np.array([p.entropy_mass]),gibbs=np.array([p.gibbs_mass]),
        cp=np.array([p.cp_mass]),cv=np.array([p.cv_mass]),mean_molecular_weight=np.array([p.mean_molecular_weight]),
        flows=np.array([ns[name].mass_flow_rate for name in ('mfc1','mfc2','outlet')]))

def run(args):
    args.output.mkdir(parents=True,exist_ok=False)
    source=args.cantera_source/'samples/python/reactors/mix1.py'
    if sha(source)!=SOURCE_SHA256 or ct.__version__!='4.0.0a2':
        raise ValueError('pinned Cantera4 mixer source/runtime required')
    ref=json.loads((args.reference/'reference.json').read_text())
    assert ref['passed'] and ref['provenance']['source']==sha(source)
    inputs=args.reference/'inputs'
    assert all(sha(inputs/name)==digest for name,digest in ref['provenance']['inputs'].items())
    environment_path=args.source_root/'validation/benchmark_environment.py'
    spec=importlib.util.spec_from_file_location('mixer_environment',environment_path)
    environment=importlib.util.module_from_spec(spec);spec.loader.exec_module(environment)
    threads_before=environment.verify_numerical_threads(set_accelerate=True)
    paths=[source,Path(__file__),environment_path,args.reference/'reference.json',args.reference/'original.npz']
    paths.extend(inputs/name for name in ref['provenance']['inputs'])
    inventory=lambda:{str(p):sha(p) for p in paths}
    libraries=lambda:{**environment.cantera_library_hashes(compiled.__file__),Path(compiled.__file__).name:sha(compiled.__file__)}
    before=inventory();library_before=libraries()
    assert library_before==ref['cantera_loaded_sha256']
    tree=ast.parse(source.read_text())
    solve=next(i for i,node in enumerate(tree.body) if isinstance(node,ast.Expr)
        and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Attribute)
        and node.value.func.attr=='solve_steady')
    tree.body=tree.body[:solve+1]
    class Prepared(ast.NodeTransformer):
        replacements=0
        def visit_Call(self,node):
            if isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name) and node.func.value.id=='ct' and node.func.attr=='Solution':
                self.replacements+=1
                assert len(node.args)==1 and node.args[0].value in ('air.yaml','gri30.yaml') and not node.keywords
                return ast.copy_location(ast.Subscript(value=ast.Name(id='prepared',ctx=ast.Load()),slice=node.args[0],ctx=ast.Load()),node)
            return self.generic_visit(node)
    transform=Prepared();tree=transform.visit(tree);ast.fix_missing_locations(tree)
    assert transform.replacements==2
    program=compile(tree,str(source),'exec')
    prepared={name:ct.Solution(str(inputs/name)) for name in ('air.yaml','gri30.yaml')}
    assert prepared['air.yaml'].n_species==8 and prepared['gri30.yaml'].n_species==53
    reference=np.load(args.reference/'original.npz',allow_pickle=False)
    def calculation():
        ns={'prepared':prepared};exec(program,ns);return ns
    times=[];saved={};report=dict(passed=False,source_before=before,libraries_before=library_before,
        threads_before=threads_before,cantera_version=ct.__version__,source_sha256=sha(source),
        host=environment.host_metadata(),scope='all original network construction and solve statements; two mechanism-file loads replaced by prepared models; printing/diagram excluded',
        repetitions=9)
    try:
        for i in range(10):
            start=time.perf_counter();ns=calculation();elapsed=time.perf_counter()-start
            if i:times.append(elapsed)
            values=snapshot(ns)
            saved.update({f'run_{i:02d}_{k}':v for k,v in values.items()})
            np.savez(args.output/'outputs.npz',**saved)
            for key,value in values.items():np.testing.assert_allclose(value,reference[key],rtol=2e-12,atol=1e-12)
        report.update(source_after=inventory(),libraries_after=libraries(),threads_after=environment.verify_numerical_threads())
        assert report['source_after']==before and report['libraries_after']==library_before
        report['passed']=True
    finally:
        report.update(warm_seconds=times,median_seconds=statistics.median(times) if times else None)
        (args.output/'timing.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'Complete Cantera mixer median: {report["median_seconds"]} s; all output replays pass')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1])
    for name in ('cantera-source','reference','output'):parser.add_argument('--'+name,type=Path,required=True)
    run(parser.parse_args())
