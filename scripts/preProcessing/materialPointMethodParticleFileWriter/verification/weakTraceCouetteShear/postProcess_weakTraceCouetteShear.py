#!/usr/bin/env python3
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_weak_trace_post import postprocess


VISIT_PLOTS = [
    {
        "name": "sigma_xy",
        "variable": "particleStressXY",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle stressXY field for the no-slip Couette shear diagnostic.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "velocity_x",
        "variable": "particleVelocityX",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle velocityX field.",
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
        "vector_count": 350,
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "tangential_jump",
        "variable": "gridInterfaceTangentialDisplacementJump",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final nodal tangential displacement jump on interface-capable nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "tangential_force",
        "variable": "gridInterfaceTangentialForce",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final nodal tangential interface force on interface-capable nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
]


if __name__ == "__main__":
    raise SystemExit(
        postprocess(
            source_dir=SOURCE_DIR,
            default_case_id="weakTraceCouetteShear",
            title="Weak Trace Couette Shear",
            result_prefix="weakTraceCouetteShear",
            visit_plots=VISIT_PLOTS,
        )
    )
