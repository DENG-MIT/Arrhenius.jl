#!/usr/bin/env python3
"""Run the shared xresid harness against this worktree's native project.

Set XRESID_PATH to the checkout containing test/cantera_parity/xresid.py.
ARRHENIUS_PROJECT points the harness's Julia subprocess at this worktree.
"""
import importlib, os, sys
from pathlib import Path

xroot_value = os.environ.get("XRESID_PATH")
if not xroot_value:
    raise SystemExit("set XRESID_PATH to the checkout containing test/cantera_parity/xresid.py")
xroot = Path(xroot_value).expanduser()
sys.path.insert(0, str(xroot / "test" / "cantera_parity"))
import xresid
xresid.PROJECT = os.environ.get("ARRHENIUS_PROJECT", str(Path(__file__).resolve().parents[3]))
xresid.main()
