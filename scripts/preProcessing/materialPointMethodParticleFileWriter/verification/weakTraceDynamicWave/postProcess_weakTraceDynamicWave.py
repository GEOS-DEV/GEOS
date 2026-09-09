#!/usr/bin/env python3
import os
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_weak_trace_post import postprocess


def float_env(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, ""))
    except Exception:
        return default


def float_list_env(name: str) -> list[float]:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return []
    values: list[float] = []
    for item in raw.replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        values.append(float(item))
    return values


LEFT_LENGTH = float_env("WEAK_TRACE_LEFT_LENGTH", 1.0)
RIGHT_LENGTH = float_env("WEAK_TRACE_RIGHT_LENGTH", 1.5)

DYNAMIC_WAVE_ANALYTICAL = {
    "source": "left_boundary_gaussian_velocity",
    "source_x": -(LEFT_LENGTH + RIGHT_LENGTH),
    "interface_x": -RIGHT_LENGTH,
    "left_density": float_env("WEAK_TRACE_LEFT_DENSITY", 1.0),
    "right_density": float_env("WEAK_TRACE_RIGHT_DENSITY", 1.0),
    "left_E": float_env("WEAK_TRACE_LEFT_E", 1.0),
    "right_E": float_env("WEAK_TRACE_RIGHT_E", 4.0),
    "poisson": float_env("WEAK_TRACE_POISSON", 0.25),
    "pulse_time_center": float_env("WEAK_TRACE_PULSE_TIME_CENTER", 0.20),
    "pulse_time_width": float_env("WEAK_TRACE_PULSE_TIME_WIDTH", 0.05),
    "pulse_velocity": float_env("WEAK_TRACE_PULSE_VELOCITY", 0.02),
    "stop_time": float_env("WEAK_TRACE_STOP_TIME", 1.6),
}
PROFILE_TIMES = float_list_env("WEAK_TRACE_PROFILE_TIMES")
if PROFILE_TIMES:
    DYNAMIC_WAVE_ANALYTICAL["profile_times"] = PROFILE_TIMES


VISIT_PLOTS = [
    {
        "name": "sigma_xx",
        "variable": "particleStressXX",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "particle stressXX field during and after wave-interface interaction.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "velocity_x",
        "variable": "particleVelocityX",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "particle velocityX field during and after wave-interface interaction.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "trace_active",
        "variable": "gridWeakInterfaceTraceActive",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "unit",
        "strict_variable": True,
        "caption": "weak-interface trace active nodes, with explicit surface-normal vectors overlaid when present.",
        "vector_variables": ["gridExplicitSurfaceNormal"],
        "vector_scale": 0.20,
        "vector_count": 350,
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "normal_jump",
        "variable": "gridInterfaceNormalDisplacementJump",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "nodal normal displacement jump on interface-capable nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "primal_mass_ratio",
        "variable": "gridWeakPrimalMassRatio",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "primal enrichment mass ratio. Low values identify small enriched modes that may limit stability or create noisy updates.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "primal_active",
        "variable": "gridWeakPrimalActive",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "unit",
        "strict_variable": True,
        "caption": "active primal enrichment nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "primal_geometry_active",
        "variable": "gridWeakPrimalGeometryActive",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "unit",
        "strict_variable": True,
        "vector_variables": ["gridWeakPrimalGeometryNormal"],
        "vector_scale": 0.08,
        "vector_count": 350,
        "caption": "cells with reconstructed primal weak-interface geometry, overlaid with reconstructed geometry normals.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "primal_geometry_distance",
        "variable": "gridWeakPrimalGeometryDistance",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "distance from grid support to the reconstructed primal weak-interface geometry.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "primal_geometry_fit_residual",
        "variable": "gridWeakPrimalGeometryFitResidual",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "least-squares residual for local reconstructed primal weak-interface geometry fits.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
    {
        "name": "primal_geometry_segment_count",
        "variable": "gridWeakPrimalGeometrySegmentCount",
        "states": "middle,final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "number of reconstructed primal weak-interface geometry segments touching each grid support.",
        "width": "0.47\\linewidth",
        "max_frames": 2,
    },
]


if __name__ == "__main__":
    raise SystemExit(
        postprocess(
            source_dir=SOURCE_DIR,
            default_case_id="weakTraceDynamicWave",
            title="Weak Trace Dynamic Wave Transmission",
            result_prefix="weakTraceDynamicWave",
            dynamic_wave=DYNAMIC_WAVE_ANALYTICAL,
            visit_plots=VISIT_PLOTS,
        )
    )
