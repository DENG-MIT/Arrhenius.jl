"""Synthetic library-policy checks; no chemistry or solver calls."""
from pathlib import Path
import tempfile,unittest
from custom2_compare import numerical_library_changes,pinned_library_receipt,PINNED_LAZY_LIBRARIES

class LibraryContract(unittest.TestCase):
    def test_strict_unchanged(self):
        self.assertTrue(numerical_library_changes({'old.so':'a'},{'old.so':'a'})['passed'])
    def test_strict_rejects_known_but_not_enabled_addition(self):
        self.assertFalse(numerical_library_changes({'old.so':'a'},{'old.so':'a',**PINNED_LAZY_LIBRARIES})['passed'])
    def test_pinned_unchanged(self):
        self.assertTrue(numerical_library_changes({'old.so':'a'},{'old.so':'a'},PINNED_LAZY_LIBRARIES)['passed'])
    def test_one_allowed_addition(self):
        name,digest=next(iter(PINNED_LAZY_LIBRARIES.items()))
        self.assertTrue(numerical_library_changes({'old.so':'a'},{'old.so':'a',name:digest},PINNED_LAZY_LIBRARIES)['passed'])
    def test_both_allowed_additions(self):
        self.assertTrue(numerical_library_changes({'old.so':'a'},{'old.so':'a',**PINNED_LAZY_LIBRARIES},PINNED_LAZY_LIBRARIES)['passed'])
    def test_removed(self):
        self.assertFalse(numerical_library_changes({'old.so':'a'},dict(PINNED_LAZY_LIBRARIES),PINNED_LAZY_LIBRARIES)['passed'])
    def test_changed(self):
        self.assertFalse(numerical_library_changes({'old.so':'a'},{'old.so':'b',**PINNED_LAZY_LIBRARIES},PINNED_LAZY_LIBRARIES)['passed'])
    def test_unlisted_addition(self):
        self.assertFalse(numerical_library_changes({'old.so':'a'},{'old.so':'a','unlisted.so':'b'},PINNED_LAZY_LIBRARIES)['passed'])
    def test_wrong_pinned_bytes(self):
        name=next(iter(PINNED_LAZY_LIBRARIES))
        self.assertFalse(numerical_library_changes({'old.so':'a'},{'old.so':'a',name:'wrong'},PINNED_LAZY_LIBRARIES)['passed'])
    def test_prior_receipt_binding(self):
        path=Path(__file__).with_name('results')/'cantera4_wsl_ic_engine.json'
        self.assertEqual(pinned_library_receipt(path)['pinned_lazy_additions'],PINNED_LAZY_LIBRARIES)
    def test_altered_receipt_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'changed.json';path.write_text('{}')
            with self.assertRaises(AssertionError):pinned_library_receipt(path)

if __name__=='__main__':unittest.main(verbosity=2)
