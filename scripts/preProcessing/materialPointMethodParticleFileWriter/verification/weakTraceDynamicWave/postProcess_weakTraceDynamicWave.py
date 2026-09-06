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
