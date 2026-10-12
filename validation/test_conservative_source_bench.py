"""Fail-closed provenance checks without running a numerical benchmark."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch
import numpy as np

with patch.dict(sys.modules,{"cantera":types.ModuleType("cantera")}):
    spec=importlib.util.spec_from_file_location("flame_bench_guard",Path(__file__).with_name("conservative_source_bench.py"))
    bench=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)


class ProvenanceChecks(unittest.TestCase):
    def test_ion_dispatch_reuses_loaded_phase_without_transport_reset(self):
        class Phase:
            @property
            def transport_model(self):return "ionized-gas"
        gas=Phase();snapshots={"frozen":{"T":np.array([600.])},"field":{"E":np.array([0.])}}
        calls=[]
        def calculation(actual,case,**kwargs):
            calls.append((actual,case,kwargs));return [.2,.1],[7,7],snapshots
        helper=types.ModuleType("ion_source_sequence");helper.run_ion_source_sequence=calculation
        with patch.dict(sys.modules,{"ion_source_sequence":helper}):
            values=bench.cantera_sequence(gas,"ion-burner",None,True)
            self.assertIs(values[2],snapshots)
            self.assertEqual(bench.cantera_sequence(gas,"ion-free",None,False)[2],{})
        self.assertEqual([c[1] for c in calls],["burner","free"])
        self.assertTrue(all(c[0] is gas and c[2]["clock_ns"] is bench.benchmark_elapsed_ns for c in calls))

    def test_ion_dispatch_rejects_wrong_transport_before_calculation(self):
        helper=types.ModuleType("ion_source_sequence")
        def unexpected(*args,**kwargs):raise AssertionError("calculation must not start")
        helper.run_ion_source_sequence=unexpected
        with patch.dict(sys.modules,{"ion_source_sequence":helper}):
            with self.assertRaisesRegex(RuntimeError,"transport model mismatch"):
                bench.cantera_sequence(types.SimpleNamespace(transport_model="mixture-averaged"),"ion-burner",None,True)

    def test_elapsed_clock_uses_raw_on_linux_and_propagates_failure(self):
        with patch.object(bench.platform,"system",return_value="Linux"), patch.object(bench.time,"clock_gettime_ns",return_value=123456789) as raw:
            self.assertEqual(bench.benchmark_elapsed_ns(),123456789)
            raw.assert_called_once_with(bench.time.CLOCK_MONOTONIC_RAW)
        with patch.object(bench.platform,"system",return_value="Linux"), patch.object(bench.time,"clock_gettime_ns",side_effect=OSError("clock unavailable")):
            with self.assertRaises((OSError,RuntimeError)):
                bench.benchmark_elapsed_ns()

    def test_elapsed_clock_preserves_nonlinux_fallback(self):
        with patch.object(bench.platform,"system",return_value="Darwin"), patch.object(bench.time,"perf_counter_ns",return_value=987654321) as mono, patch.dict(sys.modules,{"cantera":types.ModuleType("cantera")}):
            alternate=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(alternate)
            self.assertEqual(alternate.benchmark_elapsed_ns(),987654321)
            self.assertEqual(alternate.ELAPSED_CLOCK,"monotonic")
            mono.assert_called_once()

    def test_missing_or_mismatched_native_clock_rejected(self):
        self.assertEqual(bench.checked_elapsed_clock({"elapsed_clock_utf8":np.frombuffer(bench.ELAPSED_CLOCK.encode(),dtype=np.uint8)}),bench.ELAPSED_CLOCK)
        for data in [{},{"elapsed_clock_utf8":np.frombuffer(b"different",dtype=np.uint8)},{"elapsed_clock_utf8":np.array([255],dtype=np.uint8)}]:
            with self.assertRaisesRegex(RuntimeError,"elapsed clock"):
                bench.checked_elapsed_clock(data)

    def test_process_clock_observations_reject_invalid_or_missing_rows(self):
        samples=np.tile([.2,.18,0,1],(6,1))
        self.assertEqual(bench.checked_clock_observations(samples,5),samples.tolist())
        for invalid in [samples[1:],samples[:,:3],np.full((6,4),np.nan),
                np.tile([0,.18,0,1],(6,1)),np.tile([.2,-.1,0,1],(6,1)),
                np.tile([.2,.18,-1,1],(6,1)),np.tile([.2,.18,.5,1],(6,1))]:
            with self.assertRaisesRegex(RuntimeError,"invalid process-clock"):
                bench.checked_clock_observations(invalid,5)


    def test_incomplete_or_invalid_timing_samples_rejected(self):
        samples=np.ones((10,4))
        np.testing.assert_array_equal(bench.checked_times(samples,9,4),np.full(10,4.))
        for invalid in [samples[1:],samples[:,:3],np.full((10,4),np.nan),np.zeros((10,4))]:
            with self.assertRaisesRegex(RuntimeError,"invalid first/warm"):
                bench.checked_times(invalid,9,4)

    def test_changed_and_removed_input_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"source.jl"
            path.write_text("original")
            expected={str(path):bench.digest(path)}
            bench.require_unchanged(expected)
            path.write_text("changed")
            with self.assertRaisesRegex(RuntimeError,"bytes changed"):bench.require_unchanged(expected)
            path.unlink()
            with self.assertRaisesRegex(RuntimeError,"bytes changed"):bench.require_unchanged(expected)

    def test_native_file_addition_detected(self):
        with tempfile.TemporaryDirectory() as folder:
            project=Path(folder);(project/"src").mkdir()
            (project/"src/A.jl").write_text("x=1")
            expected=bench.current_source_hashes(project)
            (project/"src/B.jl").write_text("x=2")
            self.assertNotEqual(expected,bench.current_source_hashes(project))

    def test_actual_loaded_libraries_must_match_record(self):
        with tempfile.TemporaryDirectory() as folder:
            library=Path(folder)/"libcantera_shared.so";library.write_bytes(b"pristine shared")
            module=Path(folder)/"_cantera.so";module.write_bytes(b"pristine extension")
            compiled=types.ModuleType("cantera._cantera");compiled.__file__=str(module)
            cantera=types.ModuleType("cantera");cantera._cantera=compiled
            record={"library_hashes":{"shared":bench.digest(library),"extension":bench.digest(module)}}
            with patch.dict(sys.modules,{"cantera":cantera,"cantera._cantera":compiled}),patch.object(bench,"loaded_library_paths",return_value=[library]):
                self.assertEqual(len(bench.verified_cantera_libraries(record)),2)
                library.write_bytes(b"different build")
                with self.assertRaisesRegex(RuntimeError,"build-record mismatch"):bench.verified_cantera_libraries(record)

    def test_missing_shared_library_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            module=Path(folder)/"_cantera.so";module.write_bytes(b"extension only")
            compiled=types.ModuleType("cantera._cantera");compiled.__file__=str(module)
            cantera=types.ModuleType("cantera");cantera._cantera=compiled
            with patch.dict(sys.modules,{"cantera":cantera,"cantera._cantera":compiled}),patch.object(bench,"loaded_library_paths",return_value=[]):
                with self.assertRaisesRegex(RuntimeError,"build-record mismatch"):
                    bench.verified_cantera_libraries({"library_hashes":{"extension":bench.digest(module)}})

    def test_exact_prescribed_profile_values_required(self):
        with tempfile.TemporaryDirectory() as folder:
            source=Path(folder)/"flame_fixed_T.py"
            source.write_text("zloc=np.array([0.,.005,.01])\ntvalues=np.array([373.7,1000.,1589.])\n")
            mechanism=Path(folder)/"gri30.yaml";mechanism.write_text("test mechanism")
            sha=bench.digest(source)
            provenance={"source_commit":bench.PINNED_COMMIT,"examples":{"flame_fixed_T":{"sha256":sha}}}
            profile={"positions":np.array([0.,.005,.01]),"temperatures":np.array([373.7,1000.,1589.])}
            with patch.dict(bench.SOURCE_EXAMPLES,{"fixed":("flame_fixed_T",sha)}),patch.dict(bench.MECHANISM_HASHES,{"gri30":bench.digest(mechanism)}):
                bench.verify_source_inputs("fixed",source,mechanism,provenance,profile)
                profile["temperatures"][1]+=1e-9
                with self.assertRaisesRegex(ValueError,"exact source values"):
                    bench.verify_source_inputs("fixed",source,mechanism,provenance,profile)
                provenance["source_commit"]="unverified"
                with self.assertRaisesRegex(ValueError,"source/provenance"):
                    bench.verify_source_inputs("fixed",source,mechanism,provenance,profile)


if __name__=="__main__":unittest.main()
