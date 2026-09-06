#!/usr/bin/env python3
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_weak_trace_post import postprocess


VISIT_PLOTS = [
    {
        "name": "sigma_yz",
        "variable": "particleStressYZ",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle stressYZ field.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "sigma_xz",
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
        "caption": "final weak-interface trace active nodes, with explicit surface-normal vectors overlaid when present.",
        "vector_variables": ["gridExplicitSurfaceNormal"],
        "vector_scale": 0.20,
        "vector_count": 400,
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
            default_case_id="weakTraceCircularInclusionAntiPlane",
            title="Weak Trace Circular Inclusion, Anti-Plane Shear",
            result_prefix="weakTraceCircularInclusionAntiPlane",
            circular_inclusion={
                "mode": "anti_plane",
                "inclusion_radius": 0.22,
                "outer_radius": 0.60,
                "matrix_E": 1.0,
                "inclusion_E": 4.0,
                "poisson": 0.25,
                "plate_velocity_z": 0.001,
                "gap": 1.2,
            },
            visit_plots=VISIT_PLOTS,
        )
    )
