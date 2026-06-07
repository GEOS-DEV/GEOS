#!/usr/bin/env python3
"""Folder-owned post-process entry point.

The current implementation delegates to the common generic post-processor.
Keep folder-specific quantitative metrics in this file as the analytical
comparison for this test is tightened.
"""
from __future__ import annotations

import runpy
from pathlib import Path

root = Path(__file__).resolve()
for parent in [root.parent] + list(root.parents):
    candidate = parent / "mpm_vv_postprocess.py"
    if candidate.is_file():
        runpy.run_path(str(candidate), run_name="__main__")
        raise SystemExit
raise SystemExit("could not find mpm_vv_postprocess.py")
