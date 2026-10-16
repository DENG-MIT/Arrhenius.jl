#!/usr/bin/env python3
"""Run the repository's bundled xresid harness for this project checkout."""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import xresid

xresid.PROJECT = os.environ.get("ARRHENIUS_PROJECT", str(Path(__file__).resolve().parents[3]))
xresid.main()
