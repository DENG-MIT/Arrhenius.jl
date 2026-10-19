#!/usr/bin/env python3
"""Create a fixed-grid LoKI-B++ probe setup from its bundled Argon template."""
import argparse
import re
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("template", type=Path)
p.add_argument("output", type=Path)
p.add_argument("cross_sections", type=Path)
p.add_argument("result_dir", type=Path)
p.add_argument("--cells", type=int, default=2000)
p.add_argument("--emax", type=float, default=100.0)
p.add_argument("--field-td", type=float, default=100.0)
p.add_argument("--eedf-rel-error", type=float, default=1e-5)
a = p.parse_args()
s = a.template.read_text()

def replace_line(pattern, replacement):
    global s
    s, n = re.subn(pattern, replacement, s, count=1, flags=re.M)
    if n != 1:
        raise ValueError(f"expected exactly one match for {pattern!r}, got {n}")

replace_line(r"^(\s*reducedField:\s*).*$", rf"\g<1>{a.field_td:g}")
replace_line(r"^(\s*ionizationOperatorType:\s*).*$", r"\g<1>equalSharing")
replace_line(r"^(\s*- Ar_LXCat\.txt\s*)$", f"    - {a.cross_sections}")
replace_line(r"^(\s*maxEnergy:\s*).*$", rf"\g<1>{a.emax:g}")
replace_line(r"^(\s*cellNumber:\s*).*$", rf"\g<1>{a.cells}")
replace_line(r"^(\s*minEedfDecay:\s*).*$", r"\g<1>0")
replace_line(r"^(\s*maxEedfDecay:\s*).*$", r"\g<1>100")
replace_line(r"^(\s*maxEedfRelError:\s*).*$", rf"\g<1>{a.eedf_rel_error:g}")
gui_start = s.index("gui:")
output_start = s.index("output:", gui_start)
gui = s[gui_start:output_start]
gui, n = re.subn(r"^(\s*isOn:\s*)true(.*)$", r"\g<1>false\2", gui, count=1, flags=re.M)
if n != 1:
    raise ValueError("could not disable LoKI GUI")
s = s[:gui_start] + gui + s[output_start:]
output_start = s.index("output:", gui_start)
output = s[output_start:]
output, n = re.subn(r"^(\s*isOn:\s*)false(.*)$", r"\g<1>true\2", output, count=1, flags=re.M)
if n != 1:
    raise ValueError("could not enable LoKI output")
s = s[:output_start] + output
replace_line(r"^(\s*folder:\s*).*$", rf"\g<1>{a.result_dir}")
a.output.parent.mkdir(parents=True, exist_ok=True)
a.output.write_text(s)
print(f"wrote={a.output} cells={a.cells} emax_eV={a.emax:g} field_Td={a.field_td:g}")
print(f"cross_sections={a.cross_sections} output_dir={a.result_dir}")



