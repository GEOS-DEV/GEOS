#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

VARIANTS = [{"name": "identity", "label": "identity kinematics", "case_name": "materialSwapEvent_identity"}]
DEFAULT_NAME = "materialSwap"
SWAP_TIME = 0.50
STOP_TIME = 1.00
WRITE_INTERVAL = STOP_TIME / 100.0
SKIP_AROUND_SWAP = 1.5 * WRITE_INTERVAL


def as_float(value: object, default: float | None = None) -> float | None:
    try:
        return float(value)
    except Exception:
        return default


def find_column(row: dict, contains: list[str], excludes: list[str] | None = None) -> str | None:
    excludes = excludes or []
    for key in row:
        low = key.lower()
        if all(token.lower() in low for token in contains) and not any(token.lower() in low for token in excludes):
            return key
    return None


def velocity_norm(row: dict) -> float | None:
    aliases = [
        ("velocityx", "velocityy", "velocityz"),
        ("particlevelocityx", "particlevelocityy", "particlevelocityz"),
        ("vx", "vy", "vz"),
    ]
    low_to_key = {str(k).lower(): k for k in row}
    for triplet in aliases:
        if all(name in low_to_key for name in triplet):
            vals = [as_float(row[low_to_key[name]]) for name in triplet]
            if all(v is not None for v in vals):
                return math.sqrt(sum(float(v) * float(v) for v in vals))
    cols = []
    for comp in ["x", "y", "z"]:
        for key in row:
            low = key.lower()
            if "velocity" in low and (low.endswith(comp) or low.endswith("_" + comp) or low.endswith("/" + comp)):
                cols.append(key)
                break
    if len(cols) == 3:
        vals = [as_float(row[key]) for key in cols]
        if all(v is not None for v in vals):
            return math.sqrt(sum(float(v) * float(v) for v in vals))
    return None


def expected_material_type(t: float) -> int:
    return 0 if t < SWAP_TIME else 1


def process_tracer(run_dir: Path, subcase: dict) -> tuple[list[dict], dict, list[tuple[str, list[float], list[float]]]]:
    candidates = tracer_files(run_dir, "materialSwap_center")
    if not candidates:
        return [], {"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": 0, "error": "materialSwap_center tracer file not found"}, []
    data = read_tracer(candidates[0])
    if not data:
        return [], {"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": 0, "error": "empty tracer file"}, []
    material_col = find_column(data[0], ["material", "type"])
    if material_col is None:
        return [], {"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": len(data), "error": "particleMaterialType tracer column not found"}, []

    x0, y0, z0 = data[0]["x"], data[0]["y"], data[0]["z"]
    rows = []
    mismatch_count = 0
    compared = 0
    max_displacement = 0.0
    max_velocity = 0.0
    times, measured, expected, displacements = [], [], [], []
    for row in data:
        t = float(row["t"])
        mval = as_float(row.get(material_col))
        if mval is None:
            continue
        mtype = int(round(mval))
        exp_type = expected_material_type(t)
        disp = math.sqrt((float(row["x"]) - x0) ** 2 + (float(row["y"]) - y0) ** 2 + (float(row["z"]) - z0) ** 2)
        vel = velocity_norm(row)
        max_displacement = max(max_displacement, disp)
        if vel is not None:
            max_velocity = max(max_velocity, vel)
        is_transition_sample = abs(t - SWAP_TIME) <= SKIP_AROUND_SWAP
        mismatch = int((not is_transition_sample) and mtype != exp_type)
        if not is_transition_sample:
            compared += 1
            mismatch_count += mismatch
        rows.append({
            "variant": subcase.get("name"),
            "label": subcase.get("label"),
            "time": t,
            "measured_material_type": mtype,
            "expected_material_type": exp_type,
            "material_type_mismatch": mismatch,
            "displacement_error": disp,
            "velocity_norm": vel if vel is not None else "",
        })
        times.append(t)
        measured.append(float(mtype))
        expected.append(float(exp_type))
        displacements.append(disp)
    summary = {
        "variant": subcase.get("name"),
        "label": subcase.get("label"),
        "num_samples": len(rows),
        "num_compared_samples": compared,
        "material_type_mismatch_count": mismatch_count,
        "max_displacement_error": max_displacement,
        "max_velocity_norm": max_velocity,
        "transition_skip_half_width": SKIP_AROUND_SWAP,
    }
    series = [
        (str(subcase.get("label")), times, measured),
        ("expected", times, expected),
        (str(subcase.get("label")) + " displacement", times, displacements),
    ]
    return rows, summary, series


def main() -> int:
    args = parse_common_args("Post-process MaterialSwap event verification")
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)

    metric_rows = []
    summaries = []
    visit_tex = []
    material_series = []
    displacement_series = []
    visit_frames = []
    for subcase in manifest.get("subcases", []):
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "error": "run directory not found"})
            continue
        rows, summary, series = process_tracer(run_dir, subcase)
        metric_rows.extend(rows)
        summaries.append(summary)
        if len(series) >= 2:
            material_series.extend(series[:2])
        if len(series) >= 3:
            displacement_series.append(series[2])
        frames = render_visit_frames(args, run_dir, output_dir, str(subcase.get("case_name")), "particleMaterialType", states="initial,middle,final", view="auto", colortable="viridis", range_mode="explicit", color_min=0.0, color_max=1.0)
        visit_frames.extend(frames)
        visit_tex.extend(visit_frames_tex(frames, output_dir, "Particle material-type frames show the pre-event, intermediate, and post-event material assignment."))

    plot_metric(output_dir, "material_swap_type_history.png", "MaterialSwap tracer material type", "time", "material type", material_series)
    plot_metric(output_dir, "material_swap_displacement_error.png", "MaterialSwap stationary-particle error", "time", "|x - x0|", displacement_series)
    write_rows(output_dir / "material_swap_metrics.csv", metric_rows)
    write_json(output_dir / "material_swap_summary.json", {"summaries": summaries, "swap_time": SWAP_TIME, "visit_frames": visit_frames})

    tex = [r"\paragraph{Quantitative result.} The center tracer should remain fixed while its material type changes from 0 to 1 at the MaterialSwap event time. Samples within one output interval of the event are excluded from the mismatch count because event ordering and output ordering can coincide at exactly the same time."]
    tex.append(r"{\scriptsize\begin{tabular}{lrrrr}\toprule Variant & samples & compared & mismatches & max $|x-x_0|$ \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("num_samples", 0), 0) + " & " + compact_float(s.get("num_compared_samples", 0), 0) + " & " + compact_float(s.get("material_type_mismatch_count", 0), 0) + " & " + compact_float(s.get("max_displacement_error")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["material_swap_type_history.png", "material_swap_displacement_error.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    if visit_frames:
        tex.append(r"\paragraph{VisIt frames.} Particle material type is rendered at the initial, intermediate, and final states to show the region transfer.")
    tex.extend(visit_tex)
    (output_dir / "material_swap_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
