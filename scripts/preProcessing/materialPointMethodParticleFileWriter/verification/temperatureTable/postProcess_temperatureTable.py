#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

VARIANTS = [{"name": "linear", "label": "linear temperature table", "case_name": "temperatureProfileEvent_linear"}]
DEFAULT_NAME = "temperatureTable"
TEMPERATURE_TABLE = [(0.0, 300.0), (0.4, 260.0), (1.0, 380.0)]
KINK_TIMES = [0.4]
STOP_TIME = 1.0
WRITE_INTERVAL = STOP_TIME / 100.0
SKIP_AROUND_KINK = 1.5 * WRITE_INTERVAL


def as_float(value: object, default: float | None = None) -> float | None:
    try:
        return float(value)
    except Exception:
        return default


def interp_temperature(t: float) -> float:
    if t <= TEMPERATURE_TABLE[0][0]:
        return TEMPERATURE_TABLE[0][1]
    for (t0, v0), (t1, v1) in zip(TEMPERATURE_TABLE[:-1], TEMPERATURE_TABLE[1:]):
        if t <= t1:
            a = (t - t0) / max(t1 - t0, 1.0e-30)
            return v0 + a * (v1 - v0)
    return TEMPERATURE_TABLE[-1][1]


def interp_rate(t: float) -> float:
    if t <= TEMPERATURE_TABLE[0][0]:
        t0, v0 = TEMPERATURE_TABLE[0]
        t1, v1 = TEMPERATURE_TABLE[1]
        return (v1 - v0) / (t1 - t0)
    for (t0, v0), (t1, v1) in zip(TEMPERATURE_TABLE[:-1], TEMPERATURE_TABLE[1:]):
        if t <= t1:
            return (v1 - v0) / (t1 - t0)
    t0, v0 = TEMPERATURE_TABLE[-2]
    t1, v1 = TEMPERATURE_TABLE[-1]
    return (v1 - v0) / (t1 - t0)


def find_column(row: dict, contains: list[str], excludes: list[str] | None = None) -> str | None:
    excludes = excludes or []
    for key in row:
        low = key.lower()
        if all(token.lower() in low for token in contains) and not any(token.lower() in low for token in excludes):
            return key
    return None


def read_temperature_rows(run_dir: Path, subcase: dict) -> tuple[list[dict], str]:
    rows = []
    for path in tracer_files(run_dir, "temperatureProfile_center"):
        data = read_tracer(path)
        if not data:
            continue
        temp_col = find_column(data[0], ["temperature"], ["rate"])
        rate_col = find_column(data[0], ["temperature", "rate"])
        if temp_col is None:
            continue
        for row in data:
            t = float(row["t"])
            rows.append({
                "variant": subcase.get("name"),
                "label": subcase.get("label"),
                "source": "tracer",
                "time": t,
                "measured_temperature": as_float(row.get(temp_col)),
                "measured_temperature_rate": as_float(row.get(rate_col)) if rate_col else "",
            })
        if rows:
            return rows, "tracer"

    box_path = find_first(run_dir, ["boxAverageHistory.csv"])
    data = read_csv_numeric(box_path) if box_path else []
    if data:
        temp_col = find_column(data[0], ["temperature"], ["rate"])
        rate_col = find_column(data[0], ["temperature", "rate"])
        if temp_col is not None:
            for row in data:
                t = float(row.get("time", row.get("t", 0.0)))
                rows.append({
                    "variant": subcase.get("name"),
                    "label": subcase.get("label"),
                    "source": "boxAverage",
                    "time": t,
                    "measured_temperature": as_float(row.get(temp_col)),
                    "measured_temperature_rate": as_float(row.get(rate_col)) if rate_col else "",
                })
            if rows:
                return rows, "boxAverage"
    return [], "none"


def enrich_and_summarize(rows: list[dict], subcase: dict, source: str) -> tuple[list[dict], dict, list[tuple[str, list[float], list[float]]], list[tuple[str, list[float], list[float]]]]:
    if not rows:
        return [], {"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": 0, "error": "temperature tracer or boxAverage data not found"}, [], []
    temp_errors = []
    rate_errors = []
    times = []
    measured_temp = []
    expected_temp = []
    measured_rate = []
    expected_rate = []
    rate_times = []
    out = []
    for row in rows:
        t = float(row["time"])
        temp = as_float(row.get("measured_temperature"))
        if temp is None:
            continue
        exp_t = interp_temperature(t)
        err_t = temp - exp_t
        rate = as_float(row.get("measured_temperature_rate"))
        exp_r = interp_rate(t)
        near_kink = any(abs(t - kink) <= SKIP_AROUND_KINK for kink in KINK_TIMES)
        err_r = ""
        if rate is not None and not near_kink:
            err_r = rate - exp_r
            rate_errors.append(err_r)
        temp_errors.append(err_t)
        times.append(t)
        measured_temp.append(temp)
        expected_temp.append(exp_t)
        if rate is not None:
            measured_rate.append(rate)
            expected_rate.append(exp_r)
            rate_times.append(t)
        out.append({
            **row,
            "expected_temperature": exp_t,
            "temperature_error": err_t,
            "expected_temperature_rate": exp_r,
            "temperature_rate_error": err_r,
            "near_table_kink": int(near_kink),
        })
    if not out:
        return [], {"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": 0, "error": "temperature column had no numeric samples"}, [], []
    summary = {
        "variant": subcase.get("name"),
        "label": subcase.get("label"),
        "source": source,
        "num_samples": len(out),
        "rms_temperature_error": math.sqrt(sum(e * e for e in temp_errors) / len(temp_errors)),
        "max_abs_temperature_error": max(abs(e) for e in temp_errors),
        "rms_temperature_rate_error": math.sqrt(sum(e * e for e in rate_errors) / len(rate_errors)) if rate_errors else "",
        "max_abs_temperature_rate_error": max(abs(e) for e in rate_errors) if rate_errors else "",
        "rate_samples_excluding_kinks": len(rate_errors),
    }
    temp_series = [(str(subcase.get("label")), times, measured_temp), ("expected", times, expected_temp)]
    rate_series = []
    if measured_rate:
        rate_series = [(str(subcase.get("label")), rate_times, measured_rate), ("expected", rate_times, expected_rate)]
    return out, summary, temp_series, rate_series


def main() -> int:
    args = parse_common_args("Post-process TemperatureProfile table verification")
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)

    metric_rows = []
    summaries = []
    visit_tex = []
    temp_series = []
    rate_series = []
    visit_frames = []
    for subcase in manifest.get("subcases", []):
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "error": "run directory not found"})
            continue
        rows, source = read_temperature_rows(run_dir, subcase)
        out, summary, this_temp_series, this_rate_series = enrich_and_summarize(rows, subcase, source)
        metric_rows.extend(out)
        summaries.append(summary)
        temp_series.extend(this_temp_series)
        rate_series.extend(this_rate_series)
        frames = render_visit_frames(args, run_dir, output_dir, str(subcase.get("case_name")), "particleTemperature", states="initial,middle,final", view="auto", colortable="hot_desaturated", range_mode="explicit", color_min=250.0, color_max=390.0)
        visit_frames.extend(frames)
        visit_tex.extend(visit_frames_tex(frames, output_dir, "Particle temperature frames show the initial, intermediate, and final temperature-table states."))

    plot_metric(output_dir, "temperature_profile_history.png", "TemperatureProfile table interpolation", "time", "temperature", temp_series)
    plot_metric(output_dir, "temperature_rate_history.png", "TemperatureProfile rate", "time", "temperature rate", rate_series)
    plot_metric(output_dir, "temperature_profile_error.png", "Temperature interpolation error", "time", "T - T_expected", [("error", [r["time"] for r in metric_rows], [r["temperature_error"] for r in metric_rows])])
    write_rows(output_dir / "temperature_profile_metrics.csv", metric_rows)
    write_json(output_dir / "temperature_profile_summary.json", {"summaries": summaries, "temperature_table": TEMPERATURE_TABLE, "visit_frames": visit_frames})

    tex = [r"\paragraph{Quantitative result.} The expected temperature is the piecewise-linear interpolation of the solver-level temperature table. Temperature-rate errors are not scored within one output interval of the table kink."]
    tex.append(r"{\scriptsize\begin{tabular}{lrrrr}\toprule Variant & samples & RMS $T$ error & max $T$ error & RMS $\dot T$ error \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("num_samples", 0), 0) + " & " + compact_float(s.get("rms_temperature_error")) + " & " + compact_float(s.get("max_abs_temperature_error")) + " & " + compact_float(s.get("rms_temperature_rate_error")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["temperature_profile_history.png", "temperature_rate_history.png", "temperature_profile_error.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.32\linewidth]{\CaseOutputDir/" + name + "}")
    if visit_frames:
        tex.append(r"\paragraph{VisIt frames.} Particle temperature is rendered at the initial, intermediate, and final states with a fixed 250--390 K color range.")
    tex.extend(visit_tex)
    (output_dir / "temperature_profile_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
