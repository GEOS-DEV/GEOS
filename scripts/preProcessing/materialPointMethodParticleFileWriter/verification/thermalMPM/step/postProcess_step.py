#!/usr/bin/env python3
"""Post-process step: analytical/conservation comparison and plots."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from thermal_post import run_case
if __name__=='__main__':raise SystemExit(run_case('step'))
