#!/usr/bin/env python3
"""Compatibility wrapper for the renamed manual metadata extractor."""
from pathlib import Path
import runpy
import sys

script = Path(__file__).resolve().with_name("extract_mpm_metadata.py")
sys.argv[0] = str(script)
runpy.run_path(str(script), run_name="__main__")
