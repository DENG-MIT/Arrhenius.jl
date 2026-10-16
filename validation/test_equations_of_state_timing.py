"""Local benchmark safeguards; requires the Cantera/CoolProp reference environment."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
import equations_of_state_timing as timing


class EquationsOfStateBenchmarkTests(unittest.TestCase):
    def test_replay_covers_every_model_property_and_pressure(self):
        output = tuple(tuple(np.arange(1000, dtype=float) for _ in range(5)) for _ in range(3))
        first = timing.stacked(output)
        timing.check_bitwise(first, output)
        for model in range(3):
            for prop in range(5):
                changed = [[row.copy() for row in rows] for rows in output]
                changed[model][prop][-1] = np.nextafter(changed[model][prop][-1], np.inf)
                with self.subTest(model=model, property=prop), self.assertRaises(RuntimeError):
                    timing.check_bitwise(first, changed)
        with self.assertRaises(ValueError):
            timing.check_bitwise(first, output[:-1])

    def test_native_inventory_rejects_changed_missing_and_extra_sources(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root/"src/Thermo/model.jl"
            source.parent.mkdir(parents=True)
            source.write_text("original")
            meta = {"core_source_paths": "Thermo/model.jl", "core_source_sha256": timing.sha(source)}
            self.assertEqual(timing.verify_core_sources(meta, root), 1)
            source.write_text("changed")
            with self.assertRaisesRegex(RuntimeError, "source changed"):
                timing.verify_core_sources(meta, root)
            source.unlink()
            with self.assertRaisesRegex(RuntimeError, "source tree differs"):
                timing.verify_core_sources(meta, root)
            source.write_text("original")
            (root/"src/extra.jl").write_text("extra")
            with self.assertRaisesRegex(RuntimeError, "source tree differs"):
                timing.verify_core_sources(meta, root)

    def test_phase_reset_restores_temperature_density_and_composition(self):
        gas = timing.ct.Solution("h2o2.yaml")
        gas.TPX = 300., timing.ct.one_atm, "H2:2,O2:1,N2:3.76"
        pristine = timing.capture_state(gas)
        saved = pristine.copy()
        for temperature in (1800., 2300.):
            gas.TPX = temperature, 4*timing.ct.one_atm, "H2O:2,N2:3.76"
            timing.restore_state(gas, pristine)
            np.testing.assert_array_equal(gas.state, saved)
            np.testing.assert_array_equal(pristine, saved)

    def test_property_gate_rejects_partial_or_nonfinite_results(self):
        expected = np.ones((5,1000))
        self.assertTrue(timing.property_check(expected, expected)["pass"])
        self.assertFalse(timing.property_check(expected[:,:-1], expected)["pass"])
        for value in (np.nan, np.inf, 1.01):
            actual = expected.copy()
            actual[-1,-1] = value
            self.assertFalse(timing.property_check(actual, expected)["pass"])


if __name__ == "__main__":
    unittest.main()
