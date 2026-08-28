#!/usr/bin/env python3
"""Generate variants for the surface-informed polymer continuum model validation of kelF parameterization versus temperature and loading conditions."""
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = Path('/usr/workspace/crook5/GEOS_DEV/mpm') #SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from pfw_replace import *

method = ["PIC"] #, "PIC", "XPIC", "FMPM"]
precompute = ["Precomputed"] # ["OnTheFly"]

comboObj = Combination([
    ParameterSet({'method': method}, method),
    ParameterSet({"precompute": precompute}, precompute)])
combos, names = comboObj.generateCombinations()

VARIANTS = []
for i in range(comboObj.numCombinations):
    method = combos["method"][i]
    precompute = combos["precompute"][i]
    variant = {
        "name": f"{method}_{precompute}",
        "label": f"{method}_{precompute}",
        "case_name": f"brazilianDisk_shapeFunc_{method}_{precompute}",
        "env": {
            "DISK_CASE_NAME": f"brazilianDisk_shapeFunc_{method}_{precompute}",
            "DISK_VARIANT_LABEL": f"{method}_{precompute}",
            "DISK_METHOD": method,
            "DISK_PRECOMPUTE": precompute, 
        }
    }
    VARIANTS.append(variant)

print(VARIANTS)
