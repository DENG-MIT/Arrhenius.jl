#!/usr/bin/env python3
"""Run repository xresid for a snapshot pair."""
import sys
from pathlib import Path

project = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project / "test" / "cantera_parity"))
import xresid
xresid.PROJECT = str(project)
xresid.main()
