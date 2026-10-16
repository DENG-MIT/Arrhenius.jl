"""Small no-chemistry checks for qualification rejection paths."""
from copy import deepcopy
import controlled_reference as driver

library = {"libcantera.so": "bound-library"}
row = dict(pid=1, initializer_completed=1.0, source_path_matches=True, source_sha256=driver.SOURCE_SHA,
           mechanism_sha256=driver.MECHANISM_SHA, transport_model="multicomponent",
           cantera_loaded_sha256=library,
           numerical_threads={"threadpools": [{"threads": 1}], "accelerate_threading_mode": 1})
rows = [dict(row, pid=i) for i in range(4)]
driver.validate_workers(rows, library)
def rejects(candidate):
    try:
        driver.validate_workers(candidate, library)
    except RuntimeError:
        return
    raise AssertionError("invalid worker evidence accepted")
rejects(rows[:3])
rejects([row] * 4)
for key, invalid in (("initializer_completed", None), ("source_path_matches", False), ("source_sha256", "changed"),
                     ("mechanism_sha256", "changed"), ("transport_model", "mixture-averaged"),
                     ("cantera_loaded_sha256", {})):
    changed = deepcopy(rows)
    changed[0][key] = invalid
    rejects(changed)
for invalid in ({"threadpools": []}, {"threadpools": [{"threads": 2}]},
                {"threadpools": [{"threads": 1}], "accelerate_threading_mode": 0}):
    changed = deepcopy(rows)
    changed[0]["numerical_threads"] = invalid
    rejects(changed)
print("PASS: valid evidence accepted; 11 invalid worker evidence variants rejected")
