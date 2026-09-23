#!/usr/bin/env python3
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_weak_trace_post import postprocess


VISIT_PLOTS = [
    {
        "name": "material_type_xy",
        "variable": "particleMaterialType",
        "states": "initial,final",
        "view": "xy",
        "range_mode": "explicit",
        "color_min": 0.0,
        "color_max": 1.0,
        "strict_variable": True,
        "caption": "particle material type for the x-normal 3D split, viewed along z.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "material_type_xz",
        "variable": "particleMaterialType",
        "states": "initial,final",
        "view": "xz",
        "range_mode": "explicit",
        "color_min": 0.0,
        "color_max": 1.0,
        "strict_variable": True,
        "caption": "particle material type for the x-normal 3D split, viewed along y.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "surface_flag",
        "variable": "particleSurfaceFlag",
        "states": "initial,final",
        "view": "xy",
        "range_mode": "explicit",
        "color_min": 0.0,
        "color_max": 5.0,
        "strict_variable": True,
        "caption": "particle WeakDiscontinuity flags on both sides of the split plane.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "stress_xx",
        "variable": "particleStressXX",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle stressXX field.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "stress_xz",
        "variable": "particleStressXZ",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle stressXZ field.",
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
        "vector_variables": ["gridExplicitSurfaceNormal"],
        "vector_scale": 0.20,
        "vector_count": 450,
        "caption": "final 3D weak-interface trace active nodes, with explicit normal vectors.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "trace_skip_reason",
        "variable": "gridWeakInterfaceTraceSkipReason",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final trace skip/fallback reason codes.",
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
        "caption": "final normal displacement jump on interface-capable nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
]


if __name__ == "__main__":
    raise SystemExit(
        postprocess(
            source_dir=SOURCE_DIR,
            default_case_id="weakTracePlanarInterface3D",
            title="3D Weak Trace Planar Interface",
            result_prefix="weakTracePlanarInterface3D",
            visit_plots=VISIT_PLOTS,
        )
    )
