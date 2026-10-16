"""Run with the Cantera 4 benchmark Python environment; no ODE calculations."""
from pathlib import Path
import copy
import json
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import real_gas_timing as timing


def as_toml(values):
    return "\n".join(f"{json.dumps(k)} = {json.dumps(v)}" for k,v in values.items())


class ShockTubeProvenanceTests(unittest.TestCase):
    def test_pristine_phase_reset_preserves_complete_inlets(self):
        gases={phase:timing.ct.Solution("nDodecane_Reitz.yaml","nDodecane_"+phase)
               for phase in ("RK","IG")}
        states=timing.capture_phase_states(gases)
        original={phase:state.copy() for phase,state in states.items()}
        inlets={}
        for phase,gas in gases.items():
            timing.make_network(gas,1000)
            inlets[phase]=(gas.state.copy(),gas.P,gas.X.copy())
            self.assertTrue(np.array_equal(states[phase],original[phase]))
        for repetition in range(2):
            # A completed prior sweep leaves a different composition/density.
            for gas in gases.values():
                gas.TPX=2200.,40*timing.ct.one_atm,{"co2":12.,"h2o":13.,"n2":69.56}
            timing.restore_phase_states(gases,states)
            for phase,gas in gases.items():
                with self.subTest(phase=phase,repetition=repetition):
                    self.assertTrue(np.array_equal(gas.state,original[phase]))
                    timing.make_network(gas,1000)
                    self.assertTrue(np.array_equal(gas.state,inlets[phase][0]))
                    self.assertEqual(gas.P,inlets[phase][1])
                    self.assertTrue(np.array_equal(gas.X,inlets[phase][2]))
                    self.assertTrue(np.array_equal(states[phase],original[phase]))
        with self.assertRaisesRegex(ValueError,"membership"):
            timing.restore_phase_states(gases,{"RK":states["RK"]})

    def test_complete_native_inventory(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            paths=[root/"src/first.jl",root/"validation/real_gas_timing.jl"]
            paths += [root/"example/reactors"/name for name in timing.EXAMPLE_HELPERS]
            for path in paths:
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(path.name)
            hashes=timing.current_source_hashes(root)
            meta={"source_hashes_toml":as_toml(hashes),"source_hashes_after_toml":as_toml(hashes)}
            self.assertEqual(timing.verify_native_source_inventory(root,meta),hashes)
            for path in paths:
                original=path.read_bytes()
                path.write_bytes(original+b" changed")
                with self.subTest(changed=path.name),self.assertRaisesRegex(ValueError,"inventory or bytes"):
                    timing.verify_native_source_inventory(root,meta)
                path.write_bytes(original)
            added=root/"src/new-source.dat"
            added.write_text("new inventory member")
            with self.assertRaisesRegex(ValueError,"inventory or bytes"):
                timing.verify_native_source_inventory(root,meta)
            added.unlink()
            removed=paths[0];content=removed.read_bytes();removed.unlink()
            with self.assertRaisesRegex(ValueError,"inventory or bytes"):
                timing.verify_native_source_inventory(root,meta)
            removed.write_bytes(content)
            for key in meta:
                with self.subTest(missing=key),self.assertRaisesRegex(ValueError,"inventory or bytes"):
                    timing.verify_native_source_inventory(root,dict(meta,**{key:""}))
            self.assertEqual(timing.verify_native_source_inventory(root,meta),hashes)

    def test_input_bytes(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"mechanism.yaml";path.write_text("initial")
            hashes={str(path):timing.sha(path)}
            timing.require_unchanged(hashes)
            path.write_text("changed")
            with self.assertRaisesRegex(ValueError,"bytes changed"):
                timing.require_unchanged(hashes)
            path.unlink()
            with self.assertRaisesRegex(ValueError,"bytes changed"):
                timing.require_unchanged(hashes)

    def test_mapped_extension_and_shared_library(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            extension=root/"_cantera.test.so";extension.write_text("extension")
            library=root/"libcantera_shared.so";library.write_text("library")
            extra=root/"libcantera_extra.so";extra.write_text("extra")
            record={"library_hashes":{p.name:timing.sha(p) for p in (extension,library)}}
            with patch.object(timing.compiled,"__file__",str(extension)):
                with patch.object(timing,"loaded_library_paths",return_value=[extension,library]):
                    baseline=timing.mapped_cantera_hashes(record)
                    self.assertEqual(set(baseline),{str(extension),str(library)})
                    for path in (extension,library):
                        original=path.read_bytes();path.write_bytes(original+b" changed")
                        with self.subTest(changed=path.name),self.assertRaisesRegex(ValueError,"build record"):
                            timing.mapped_cantera_hashes(record)
                        path.write_bytes(original)
                    self.assertEqual(timing.mapped_cantera_hashes(record),baseline)
                for mapped in ([library],[extension],[]):
                    with patch.object(timing,"loaded_library_paths",return_value=mapped):
                        with self.assertRaisesRegex(ValueError,"actually mapped"):
                            timing.mapped_cantera_hashes(record)
                with patch.object(timing,"loaded_library_paths",return_value=[extension,library,extra]):
                    with self.assertRaisesRegex(ValueError,"build record"):
                        timing.mapped_cantera_hashes(record)

    def test_native_mechanism_before_after_membership(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for name in timing.NATIVE_MECHANISMS:
                (root/name).write_text(name)
            hashes={name:timing.sha(root/name) for name in timing.NATIVE_MECHANISMS}
            for changed in ({},dict(hashes,unexpected="0"),dict(hashes,**{"dodecane_RK.yaml":"wrong"})):
                for key in ("mechanism_hashes_toml","mechanism_hashes_after_toml"):
                    meta={"mechanism_hashes_toml":as_toml(hashes),"mechanism_hashes_after_toml":as_toml(hashes)}
                    meta[key]=as_toml(changed)
                    with self.subTest(field=key,hashes=changed),self.assertRaisesRegex(ValueError,"differ"):
                        timing.verify_native_mechanisms(root,meta,{})

    def test_replay_fails_closed(self):
        one={"RK_1000":dict(time=np.array([1.]),state=np.array([[2.]]),final_state=np.array([3.]),
                            delay=1.,steps=20,final_time=1.1)}
        self.assertTrue(timing.require_replay(one,copy.deepcopy(one),"test"))
        for field in ("time","state","final_state","delay","steps","final_time"):
            changed=copy.deepcopy(one);changed["RK_1000"][field]+=1
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,"differs"):
                timing.require_replay(changed,one,"test")
        for value in (-1.,0.,np.nan,np.inf):
            with self.subTest(duration=value),self.assertRaisesRegex(ValueError,"positive"):
                timing.native_replay_checks({"warm_seconds":np.array([value]),"warm_matches_first":np.array([1])})
        with self.assertRaisesRegex(ValueError,"replay"):
            timing.native_replay_checks({"warm_seconds":np.array([1.]),"warm_matches_first":np.array([0])})


if __name__=="__main__":
    unittest.main()
