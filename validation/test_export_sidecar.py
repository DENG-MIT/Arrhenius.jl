"""Import provenance checks; run with the Cantera preprocessing environment."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ruamel.yaml import YAML

_spec = importlib.util.spec_from_file_location("export_sidecar", Path(__file__).resolve().parents[1] / "mechanism/export_sidecar.py")
exporter = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(exporter)


class ImportProvenanceTests(unittest.TestCase):
    def test_selected_phase_dependencies_and_search_order(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            local = root / "model"
            local.mkdir()
            data = root / "data"
            data.mkdir()
            mechanism = local / "selected.yaml"
            (local / "shared.yaml").write_text("local bytes\n")
            (data / "shared.yaml").write_text("fallback bytes\n")
            (root / "parent.yaml").write_text("parent bytes\n")
            document = {"phases": [
                {"species": [{"shared.yaml/species": ["B", "A"]}],
                 "reactions": [{"shared.yaml/reactions": "declared-species"},
                               {"../parent.yaml/named-reactions": "all"}]},
                {"species": [{"missing.yaml/species": "all"}]}]}
            with mechanism.open("w") as stream:
                YAML().dump(document, stream)
            gas = SimpleNamespace(species_names=["B", "A"], n_reactions=7)
            with patch.object(exporter.ct, "get_data_directories", return_value=[str(data)]):
                selected = json.loads(exporter.phase_selection(mechanism, gas))
                self.assertEqual(selected["species_names"], ["B", "A"])
                self.assertEqual(selected["n_reactions"], 7)
                self.assertEqual(set(selected["dependencies"]), {"shared.yaml", "../parent.yaml"})
                self.assertEqual(selected["dependencies"]["shared.yaml"],
                                 hashlib.sha256((local / "shared.yaml").read_bytes()).hexdigest())
                (local / "shared.yaml").unlink()
                selected = json.loads(exporter.phase_selection(mechanism, gas))
                self.assertEqual(selected["dependencies"]["shared.yaml"],
                                 hashlib.sha256((data / "shared.yaml").read_bytes()).hexdigest())
                self.assertEqual(exporter.resolve_import(mechanism, str(root / "parent.yaml")), root / "parent.yaml")
                (data / "shared.yaml").unlink()
                with self.assertRaises(ValueError):
                    exporter.phase_selection(mechanism, gas)

    def test_local_phase_has_no_import_dependencies(self):
        with tempfile.TemporaryDirectory() as folder:
            mechanism = Path(folder) / "local.yaml"
            mechanism.write_text("phases: [{species: [A, B], reactions: all}]\n")
            gas = SimpleNamespace(species_names=["A", "B"], n_reactions=0)
            selected = json.loads(exporter.phase_selection(mechanism, gas))
            self.assertEqual(selected, {"species_names": ["A", "B"], "n_reactions": 0, "dependencies": {}})


if __name__ == "__main__":
    unittest.main()
