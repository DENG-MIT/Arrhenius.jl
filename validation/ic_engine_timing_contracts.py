"""Small engine timing contracts: no Cantera reactor construction or solve."""
import ast
import copy
import inspect
import tempfile
import unittest
from pathlib import Path
import numpy as np
import ic_engine_timing as timing
import ic_engine_source as source
import ic_engine_qndf_case as reference

class EngineTimingContracts(unittest.TestCase):
    def test_unique_reference_preload_preserves_strict_inventory(self):
        timing.prepare_reference_imports()
        before=timing.dependency_snapshot()
        for name in ("numpy.ma","numpy.ma.core","numpy.ma.extras"):
            self.assertIn(name,before["loaded_module_files"])
        np.unique(np.array([0.,1.,0.,2.]))
        self.assertTrue(timing.require_inventory(timing.dependency_snapshot(),before,"unique preparation"))
        changed=copy.deepcopy(before)
        changed["loaded_module_files"]["unexpected_late_module"]={"path":"unknown","sha256":"unknown"}
        with self.assertRaises(ValueError):timing.require_inventory(changed,before,"late module")

    def test_no_circular_reference_import(self):
        tree=ast.parse(inspect.getsource(reference))
        imports=[n.module for n in tree.body if isinstance(n,ast.ImportFrom)]
        self.assertIn("ic_engine_source",imports)
        self.assertNotIn("ic_engine_timing",imports)
        self.assertIs(timing.reference.compare,reference.compare)

    def test_fresh_source_namespace(self):
        setup=compile("history=[]\nidentity=object()","<fixture>","exec")
        first=source.new_source_state((setup,));first["history"].append("previous solve")
        second=source.new_source_state((setup,))
        self.assertEqual(second["history"],[])
        self.assertIsNot(first["identity"],second["identity"])
        self.assertIsNot(first["history"],second["history"])

    def test_exact_inventory_keys_values_and_dependencies(self):
        expected={"source":{"entry.jl":"abc"},"dependencies":{"pkg":"one"}}
        self.assertTrue(timing.require_inventory(copy.deepcopy(expected),expected,"same"))
        for actual in ({}, {**expected,"added.jl":"extra"},
                       {**expected,"source":{}},{**expected,"dependencies":{"pkg":"two"}}):
            with self.subTest(actual=actual),self.assertRaises(ValueError):
                timing.require_inventory(actual,expected,"fixture")

    def test_save_source_before_replay_rejection(self):
        result={"time":np.linspace(1e-6,.16,2881),"Y":np.ones((2881,1)),
                "integrals":{"heat":1.},"solver_stats":{"steps":42},"rtol":1e-12,"atol":1e-16}
        self.assertTrue(timing.require_source_result(result))
        changed=copy.deepcopy(result);changed["solver_stats"]["steps"]+=1
        with tempfile.TemporaryDirectory() as directory:
            saved=timing.save_source(changed,Path(directory),"warm-1")
            with self.assertRaises(ValueError):timing.require_source_result(changed,result)
            self.assertTrue(Path(saved["array_path"]).is_file())
            self.assertEqual(timing.sha(saved["array_path"]),saved["array_sha256"])
            with self.assertRaises(ValueError):timing.save_source(result,Path(directory),"warm-1")

    def test_source_replay_all_outputs_and_counters(self):
        first={"a":np.array([0.,1.]),"integrals":{"Q":2.},"solver_stats":{"steps":42},"rtol":1e-12,"atol":1e-16}
        self.assertTrue(source.output_replay(first,copy.deepcopy(first)))
        for key,value in (("a",np.array([0.,2.])),("integrals",{"Q":3.}),("solver_stats",{"steps":43}),
                          ("rtol",1e-11),("atol",1e-15)):
            changed={**first,key:value}
            with self.subTest(key=key):self.assertFalse(source.output_replay(changed,first))
        for changed in ({**first,"added":1}, {**first,"a":np.array([-0.,1.])},
                        {**first,"a":np.array([0.,1.],dtype=np.float32)}):
            self.assertFalse(timing.exact_source_equal(changed,first))

    def test_native_dependency_and_pinned_library_current_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/"src").mkdir();file=root/"src/a.jl";file.write_text("original")
            library=root/"libfixture.so";library.write_bytes(b"original library")
            artifact=root/"input.yaml";artifact.write_text("mechanism")
            snapshot=dict(environment={"threads":"1"},settings={"blas":1},mapped=[str(library)],
                          mapped_sha256={str(library):timing.sha(library)})
            record=dict(input_before=dict(file_paths={"mechanism":str(artifact)},mechanism_sha256=timing.sha(artifact),
                dependencies={"fixture":dict(path=str(root),files=timing.native_dependency_files(root))}),
                library_pins={str(library):timing.sha(library)},runtime_before=snapshot,runtime_after=snapshot,
                runs=[dict(runtime_before=snapshot,runtime_after=snapshot)])
            self.assertTrue(timing.require_native_files(record))
            file.write_text("changed")
            with self.assertRaises(ValueError):timing.require_native_files(record)
            file.write_text("original");extra=root/"src/added.jl";extra.write_text("new membership")
            with self.assertRaises(ValueError):timing.require_native_files(record)
            extra.unlink();library.write_bytes(b"different library")
            with self.assertRaises(ValueError):timing.require_native_files(record)

    def test_saved_source_bytes_rechecked(self):
        with tempfile.TemporaryDirectory() as directory:
            saved=timing.save_source({"time":np.array([.16]),"solver_stats":{"steps":1}},Path(directory),"first")
            self.assertTrue(timing.require_saved_source(saved))
            Path(saved["array_path"]).with_suffix(".json").write_text("changed")
            with self.assertRaises(ValueError):timing.require_saved_source(saved)

    def test_native_identity_and_mode_rejections_precede_io(self):
        template={"complete":True,"mode":"smoke","warm_repetitions":1,
                  "public_callable":"NativeEngineQNDF.solve_ic_engine_qndf",
                  "wall_clock_exclusions":["seconds_including_first_specialization"],
                  "runs":[],"warm_seconds":[]}
        for mutation in ({"complete":False},{"mode":"controlled"},{"public_callable":"native_engine_sdirk"},
                         {"wall_clock_exclusions":["all_time_fields"]},{}):
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                timing.require_native_timing({**template,**mutation},"smoke",1,Path("missing"),Path("missing.npz"))

if __name__=="__main__":unittest.main()
