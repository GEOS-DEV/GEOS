#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

VARIANTS = [{"name": "singleMaterial2d", "label": "Single-material 2D disk", "case_name": "partitionCrossing2D_singleMaterial"}]
DEFAULT_NAME = "diskThruPartitions"
VELOCITY = (0.20, 0.16, 0.0)


def tracer_velocity(row: dict, axis: str) -> float | None:
    aliases = {
        "x": ["velocityx", "velx", "vx", "v_x", "velocity_x"],
        "y": ["velocityy", "vely", "vy", "v_y", "velocity_y"],
        "z": ["velocityz", "velz", "vz", "v_z", "velocity_z"],
    }[axis]
    for key, value in row.items():
        low = key.lower()
        if low in aliases or any(low.endswith(suffix) for suffix in aliases):
            try:
                return float(value)
            except Exception:
                pass
    return None


def main() -> int:
    args = parse_common_args("Post-process partition-crossing free-flight verification")
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)

    rows_out: list[dict] = []
    summaries: list[dict] = []
    visit_tex = []
    position_series: list[tuple[str, list[float], list[float]]] = []
    velocity_series: list[tuple[str, list[float], list[float]]] = []

    for subcase in manifest.get("subcases", []):
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": subcase.get("name"), "error": "run directory not found"})
            continue

        by_time: dict[float, list[tuple[float, float, float | None, float | None]]] = {}
        for tracer_file in tracer_files(run_dir, "partitionCrossingTracer"):
            try:
                rows = read_tracer(tracer_file)
            except Exception as exc:
                rows_out.append({"variant": subcase.get("name"), "tracer": tracer_file.name, "error": str(exc)})
                continue
            if not rows:
                continue
            x0 = rows[0]["x"]
            y0 = rows[0]["y"]
            for row in rows:
                t = float(row["t"])
                ex = x0 + VELOCITY[0] * t
                ey = y0 + VELOCITY[1] * t
                dx = float(row["x"]) - ex
                dy = float(row["y"]) - ey
                vx = tracer_velocity(row, "x")
                vy = tracer_velocity(row, "y")
                by_time.setdefault(t, []).append((dx, dy, vx, vy))
                rows_out.append({
                    "variant": subcase.get("name"),
                    "label": subcase.get("label"),
                    "tracer": tracer_file.name,
                    "time": t,
                    "measured_x": row["x"],
                    "measured_y": row["y"],
                    "expected_x": ex,
                    "expected_y": ey,
                    "position_error_x": dx,
                    "position_error_y": dy,
                    "position_error_norm": math.hypot(dx, dy),
                    "measured_vx": vx if vx is not None else "",
                    "measured_vy": vy if vy is not None else "",
                })

        times = sorted(by_time)
        rms_pos: list[float] = []
        rms_vel: list[float] = []
        for t in times:
            samples = by_time[t]
            rms_pos.append(math.sqrt(sum(dx * dx + dy * dy for dx, dy, _vx, _vy in samples) / max(len(samples), 1)))
            vel_terms = [(vx - VELOCITY[0]) ** 2 + (vy - VELOCITY[1]) ** 2 for _dx, _dy, vx, vy in samples if vx is not None and vy is not None]
            rms_vel.append(math.sqrt(sum(vel_terms) / len(vel_terms)) if vel_terms else float("nan"))

        if times:
            finite_vel = [v for v in rms_vel if math.isfinite(v)]
            summary = {
                "variant": subcase.get("name"),
                "label": subcase.get("label"),
                "num_time_samples": len(times),
                "num_tracer_records": sum(len(v) for v in by_time.values()),
                "rms_position_error_final": rms_pos[-1],
                "max_rms_position_error": max(rms_pos),
            }
            if finite_vel:
                summary.update({"max_rms_velocity_error": max(finite_vel), "rms_velocity_error_final": finite_vel[-1]})
            summaries.append(summary)
            position_series.append((str(subcase.get("label")), times, rms_pos))
            velocity_series.append((str(subcase.get("label")), times, rms_vel))
        else:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "num_time_samples": 0, "error": "no tracer data found"})

        frames = render_visit_frames(args, run_dir, output_dir, str(subcase.get("case_name")), "Velocity", states="initial,middle,final", view="auto", colortable="bluehot", range_mode="auto")
        visit_tex.extend(visit_frames_tex(frames, output_dir, "Particle velocity frames for the partition-crossing subcase show the initial, intermediate, and final states."))

    plot_metric(output_dir, "partition_crossing_position_error.png", "Partition crossing position error", "time", "RMS position error", position_series)
    plot_metric(output_dir, "partition_crossing_velocity_error.png", "Partition crossing velocity error", "time", "RMS velocity error", velocity_series)
    write_rows(output_dir / "partition_crossing_metrics.csv", rows_out)
    write_json(output_dir / "partition_crossing_summary.json", {"summaries": summaries, "velocity": VELOCITY})

    tex = [r"\paragraph{Quantitative result.} Tracer positions are compared with the free-flight solution $\mathbf{x}(t)=\mathbf{x}_0+\mathbf{v}_0t$ while particles migrate across the 3 by 3 processor layout."]
    tex.append(r"{\scriptsize\begin{tabular}{lrrr}\toprule Variant & records & final RMS position error & max RMS position error \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("num_tracer_records", 0), 0) + " & " + compact_float(s.get("rms_position_error_final")) + " & " + compact_float(s.get("max_rms_position_error")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["partition_crossing_position_error.png", "partition_crossing_velocity_error.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    tex.extend(visit_tex)
    (output_dir / "partition_crossing_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
