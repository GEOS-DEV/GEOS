#!/usr/bin/env python3
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_weak_trace_post import postprocess


VISIT_PLOTS = [
    {
        "name": "material_type",
        "variable": "particleMaterialType",
        "states": "initial,final",
        "view": "xy",
        "range_mode": "explicit",
        "color_min": 0.0,
        "color_max": 1.0,
        "strict_variable": True,
        "caption": "particle material type: platen material 0 and fluid material 1.",
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
        "caption": "particle surface flags marking the two platen-fluid interfaces.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "particle_rank",
        "variable": "particleRank",
        "states": "initial,final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "particle owning rank for the periodic 3-by-3 decomposition diagnostic.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
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
        "name": "deformation_gradient_xy",
        "variable": "particleDeformationGradientXY",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final particle deformation-gradient Fxy component for shear-strain diagnostics.",
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


COUETTE_REFERENCE = {
    "fluid_ymin": -0.20,
    "fluid_ymax": 0.20,
    "top_plate_velocity": 0.002,
    "fluid_viscosity": 0.1,
}


if __name__ == "__main__":
    raise SystemExit(
        postprocess(
            source_dir=SOURCE_DIR,
            default_case_id="weakTraceCouetteShear",
            title="Weak Trace Couette Shear",
            result_prefix="weakTraceCouetteShear",
            couette_shear=COUETTE_REFERENCE,
            visit_plots=VISIT_PLOTS,
        )
    )
