#!/usr/bin/env python3
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_weak_trace_post import postprocess


VISIT_PLOTS = [
    {
        "name": "sigma_xx",
        "variable": "particleStressXX",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle stressXX field for the elliptical-inclusion diagnostic.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "sigma_yy",
        "variable": "particleStressYY",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle stressYY field.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "sigma_xy",
        "variable": "particleStressXY",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle stressXY field.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "trace_active",
        "variable": "gridWeakInterfaceTraceActive",
        "states": "final",
        "view": "xy",
        "range_mode": "unit",
        "strict_variable": True,
        "caption": "final weak-interface trace active nodes, with explicit surface-normal vectors overlaid when present.",
        "vector_variables": ["gridExplicitSurfaceNormal"],
        "vector_scale": 0.20,
        "vector_count": 400,
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "normal_jump",
        "variable": "gridInterfaceNormalDisplacementJump",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final nodal normal displacement jump on interface-capable nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
]


if __name__ == "__main__":
    raise SystemExit(
        postprocess(
            source_dir=SOURCE_DIR,
            default_case_id="weakTraceEllipticalInclusion",
            title="Weak Trace Elliptical Inclusion",
            result_prefix="weakTraceEllipticalInclusion",
            visit_plots=VISIT_PLOTS,
        )
    )
