#!/usr/bin/env python3
from pathlib import Path
import sys
import math

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_weak_trace_post import postprocess


def float_env(name: str, default: float) -> float:
    import os
    try:
        return float(os.environ.get(name, ""))
    except Exception:
        return default


INCLUSION_RADIUS = float_env("WEAK_TRACE_INCLUSION_RADIUS", 0.22)
OUTER_RADIUS = float_env("WEAK_TRACE_OUTER_RADIUS", 0.72)
MATRIX_E = float_env("WEAK_TRACE_MATRIX_E", 1.0)
INCLUSION_E = float_env("WEAK_TRACE_INCLUSION_E", 4.0)
POISSON = float_env("WEAK_TRACE_POISSON", 0.25)
EXTERNAL_PRESSURE = float_env("WEAK_TRACE_EXTERNAL_PRESSURE", 0.004)
STOP_TIME = float_env("WEAK_TRACE_STOP_TIME", 8.0)
REFINE = int(float_env("WEAK_TRACE_REFINE", 1))
CPP = int(float_env("WEAK_TRACE_CPP", 64))
DOMAIN_HALF_WIDTH = float_env("WEAK_TRACE_DOMAIN_HALF_WIDTH", 1.15 * OUTER_RADIUS)
CELL_SIZE = 2.0 * DOMAIN_HALF_WIDTH / float(max(1, REFINE * CPP - 2))
SURFACE_TRACER_OFFSET = float_env("WEAK_TRACE_SURFACE_TRACER_OFFSET", 0.5 * CELL_SIZE)


def radial_stress_range() -> tuple[float, float]:
    shear_m = MATRIX_E / (2.0 * (1.0 + POISSON))
    h_m = MATRIX_E * POISSON / ((1.0 + POISSON) * (1.0 - 2.0 * POISSON)) + shear_m
    h_i = INCLUSION_E * POISSON / ((1.0 + POISSON) * (1.0 - 2.0 * POISSON)) + INCLUSION_E / (2.0 * (1.0 + POISSON))
    a = INCLUSION_RADIUS
    b = OUTER_RADIUS
    denominator = 1.0 / (h_m * b * b) + 1.0 / (shear_m * a * a) - (1.0 / (b * b) - 1.0 / (a * a)) / h_i
    d_m = EXTERNAL_PRESSURE * (1.0 / h_m - 1.0 / h_i) / denominator
    c_m = -EXTERNAL_PRESSURE + d_m / (b * b)
    c_i = c_m - d_m / (a * a)
    samples = [c_i]
    for i in range(200):
        r = a + (b - a) * i / 199.0
        samples.extend([c_m - d_m / (r * r), c_m + d_m / (r * r)])
    return min(samples), max(samples)


STRESS_MIN, STRESS_MAX = radial_stress_range()


VISIT_PLOTS = [
    {
        "name": "sigma_xx",
        "variable": "particleStressXX",
        "states": "final",
        "view": "xy",
        "range_mode": "explicit",
        "color_min": STRESS_MIN,
        "color_max": STRESS_MAX,
        "strict_variable": True,
        "caption": "final particle stressXX field.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "sigma_yy",
        "variable": "particleStressYY",
        "states": "final",
        "view": "xy",
        "range_mode": "explicit",
        "color_min": STRESS_MIN,
        "color_max": STRESS_MAX,
        "strict_variable": True,
        "caption": "final particle stressYY field.",
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
        "caption": "final weak-interface trace active anchor nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "trace_skip_reason",
        "variable": "gridWeakInterfaceTraceSkipReason",
        "states": "final",
        "view": "xy",
        "range_mode": "explicit",
        "color_min": 0.0,
        "color_max": 10.0,
        "strict_variable": True,
        "caption": "final weak-interface trace skip reason.  Code 10 marks active projection anchors.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "trace_anchor_weight",
        "variable": "gridWeakInterfaceTraceAnchorWeight",
        "states": "final",
        "view": "xy",
        "range_mode": "explicit",
        "color_min": 0.0,
        "color_max": 1.0,
        "strict_variable": True,
        "caption": "final weak-interface trace interpolation weight at the anchor node.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "trace_support_weight",
        "variable": "gridWeakInterfaceTraceSupportWeight",
        "states": "final",
        "view": "xy",
        "range_mode": "explicit",
        "color_min": 0.0,
        "color_max": 1.0,
        "strict_variable": True,
        "caption": "final active support weight for the weak-interface trace projection.",
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
    {
        "name": "normal_velocity_jump",
        "variable": "gridInterfaceNormalVelocityJump",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final post-projection normal velocity jump on interface-capable nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "tangential_velocity_jump",
        "variable": "gridInterfaceTangentialVelocityJump",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final post-projection tangential velocity jump magnitude on interface-capable nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
    {
        "name": "normal_force",
        "variable": "gridInterfaceNormalForce",
        "states": "final",
        "view": "xy",
        "range_mode": "auto",
        "strict_variable": True,
        "caption": "final nodal normal interface force on interface-capable nodes.",
        "width": "0.47\\linewidth",
        "max_frames": 1,
    },
]


if __name__ == "__main__":
    raise SystemExit(
        postprocess(
            source_dir=SOURCE_DIR,
            default_case_id="weakTraceCircularInclusionRadial",
            title="Weak Trace Circular Inclusion, Radial Loading",
            result_prefix="weakTraceCircularInclusionRadial",
            circular_inclusion={
                "mode": "radial",
                "inclusion_radius": INCLUSION_RADIUS,
                "outer_radius": OUTER_RADIUS,
                "matrix_E": MATRIX_E,
                "inclusion_E": INCLUSION_E,
                "poisson": POISSON,
                "start_pressure": 0.0,
                "external_pressure": EXTERNAL_PRESSURE,
                "load_start_time": 0.0,
                "load_end_time": STOP_TIME,
                "stop_time": STOP_TIME,
                "load_interp": "cosine",
                "axis_tolerance": 0.35 * CELL_SIZE,
                "ring_tolerance": 1.5 * SURFACE_TRACER_OFFSET,
            },
            visit_plots=VISIT_PLOTS,
        )
    )
