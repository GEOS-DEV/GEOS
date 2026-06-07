#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

VARIANTS = [{"name": "disk2d", "label": "2D periodic disk", "case_name": "periodicAdvection_disk2d"}]
DEFAULT_NAME = "periodicBoundaries"
CENTER0 = (-0.18, -0.16, 0.0)
VELOCITY = (0.42, 0.31, 0.0)
DOMAIN = (-0.5, 0.5, -0.5, 0.5)


def wrap_coordinate(x: float, lo: float, hi: float) -> float:
    length = hi - lo
    return lo + ((x - lo) % length)


def minimum_image_error(measured: float, expected: float, lo: float, hi: float) -> float:
    length = hi - lo
    return ((measured - expected + 0.5 * length) % length) - 0.5 * length


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
    args = parse_common_args("Post-process periodic advection verification")
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)

    rows_out: list[dict] = []
    summaries: list[dict] = []
    rms_series: list[tuple[str, list[float], list[float]]] = []
    velocity_series: list[tuple[str, list[float], list[float]]] = []

    for subcase in manifest.get("subcases", []):
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": subcase.get("name"), "error": "run directory not found"})
            continue

        by_time: dict[float, list[tuple[float, float, float | None, float | None]]] = {}
        for tracer_file in tracer_files(run_dir, "periodicAdvectionTracer"):
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
                ex = wrap_coordinate(x0 + VELOCITY[0] * t, DOMAIN[0], DOMAIN[1])
                ey = wrap_coordinate(y0 + VELOCITY[1] * t, DOMAIN[2], DOMAIN[3])
                dx = minimum_image_error(float(row["x"]), ex, DOMAIN[0], DOMAIN[1])
                dy = minimum_image_error(float(row["y"]), ey, DOMAIN[2], DOMAIN[3])
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
                    "expected_x_modulo": ex,
                    "expected_y_modulo": ey,
                    "periodic_error_x": dx,
                    "periodic_error_y": dy,
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
            if vel_terms:
                rms_vel.append(math.sqrt(sum(vel_terms) / len(vel_terms)))
            else:
                rms_vel.append(float("nan"))
        if times:
            finite_vel = [v for v in rms_vel if math.isfinite(v)]
            summary = {
                "variant": subcase.get("name"),
                "label": subcase.get("label"),
                "num_time_samples": len(times),
                "num_tracer_records": sum(len(v) for v in by_time.values()),
                "rms_periodic_position_error_final": rms_pos[-1],
                "max_rms_periodic_position_error": max(rms_pos),
            }
            if finite_vel:
                summary.update({"max_rms_velocity_error": max(finite_vel), "rms_velocity_error_final": finite_vel[-1]})
            summaries.append(summary)
            rms_series.append((str(subcase.get("label")), times, rms_pos))
            velocity_series.append((str(subcase.get("label")), times[:len(rms_vel)], rms_vel))
        else:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "num_time_samples": 0, "error": "no tracer data found"})

        render_visit_frames(args, run_dir, output_dir, str(subcase.get("case_name")), "Velocity", states="initial,middle,final", view="auto", colortable="bluehot", range_mode="auto")

    plot_metric(output_dir, "periodic_position_error.png", "Periodic advection position error", "time", "RMS minimum-image position error", rms_series)
    plot_metric(output_dir, "periodic_velocity_error.png", "Periodic advection velocity error", "time", "RMS velocity error", velocity_series)
    write_rows(output_dir / "periodic_advection_metrics.csv", rows_out)
    write_json(output_dir / "periodic_advection_summary.json", {"summaries": summaries, "velocity": VELOCITY, "domain": DOMAIN})

    tex = [r"\paragraph{Quantitative result.} Tracer positions are compared with $x(t)=x_0+v_xt$ and $y(t)=y_0+v_yt$ using a minimum-image distance on the periodic unit square."]
    tex.append(r"{\scriptsize\begin{tabular}{lrrr}\toprule Variant & records & final RMS position error & max RMS position error \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("num_tracer_records", 0), 0) + " & " + compact_float(s.get("rms_periodic_position_error_final")) + " & " + compact_float(s.get("max_rms_periodic_position_error")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["periodic_position_error.png", "periodic_velocity_error.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    (output_dir / "periodic_advection_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
