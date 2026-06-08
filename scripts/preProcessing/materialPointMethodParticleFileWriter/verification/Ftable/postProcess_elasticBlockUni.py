#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

VARIANTS = [{"name": "smoothstep", "label": "Smoothstep F-table", "case_name": "elasticFTableBoundarySwitch_smoothstep"}]
DEFAULT_NAME = "Ftable"
E = 1.0
NU = 0.25
M = E * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU))


def stretch_x(t: float) -> float:
    table = [(0.0, 1.000), (0.5, 1.010), (1.0, 1.020), (1.25, 1.000)]
    if t <= table[0][0]:
        return table[0][1]
    for (t0, v0), (t1, v1) in zip(table[:-1], table[1:]):
        if t <= t1:
            a = (t - t0) / max(t1 - t0, 1.0e-30)
            return v0 + a * (v1 - v0)
    return table[-1][1]


def expected_stress(t: float) -> float:
    eps = stretch_x(t) - 1.0
    if t <= 0.5:
        return M * eps
    if t <= 1.0:
        return M * 0.010 + E * (eps - 0.010)
    return M * 0.010 + E * (0.020 - 0.010) + E * (eps - 0.020)


def pick_stress_x(row: dict) -> float | None:
    keys = list(row.keys())
    preferred = []
    for k in keys:
        low = k.lower()
        if "stress" in low and any(tok in low for tok in ["xx", "11", "_0", "x"]):
            preferred.append(k)
    for k in preferred + keys:
        if "stress" not in k.lower():
            continue
        try:
            return float(row[k])
        except Exception:
            pass
    return None


def main() -> int:
    args = parse_common_args("Post-process elastic F-table verification")
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)
    rows_out = []
    summaries = []
    visit_tex = []
    plot_series = []
    expected_series = None
    for subcase in manifest.get("subcases", []):
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": subcase.get("name"), "error": "run directory not found"})
            continue
        box_path = find_first(run_dir, ["boxAverageHistory.csv"])
        data = read_csv_numeric(box_path) if box_path else []
        times, measured, expected, errors = [], [], [], []
        for row in data:
            t = float(row.get("time", row.get("t", 0.0)))
            sx = pick_stress_x(row)
            if sx is None:
                continue
            ex = expected_stress(t)
            times.append(t); measured.append(sx); expected.append(ex); errors.append(sx - ex)
            rows_out.append({"variant": subcase.get("name"), "label": subcase.get("label"), "time": t, "measured_sigma_x": sx, "expected_sigma_x": ex, "error_sigma_x": sx - ex})
        if times:
            rms = math.sqrt(sum(e*e for e in errors) / len(errors))
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": len(times), "rms_sigma_x_error": rms, "max_abs_sigma_x_error": max(abs(e) for e in errors)})
            plot_series.append((str(subcase.get("label")), times, measured))
            expected_series = ("expected", times, expected)
        else:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": 0, "error": "boxAverageHistory stress column not found"})
        frames = render_visit_frames(args, run_dir, output_dir, str(subcase.get("case_name")), "Stress", states="initial,middle,final", view="auto", colortable="hot_desaturated", range_mode="auto")
        visit_tex.extend(visit_frames_tex(frames, output_dir, "Particle stress frames for a representative F-table subcase show the initial, intermediate, and final states."))
    if expected_series:
        plot_series.insert(0, expected_series)
    plot_metric(output_dir, "elastic_f_table_stress_error.png", "Elastic F-table boundary switch", "time", "sigma_x", plot_series)
    plot_metric(output_dir, "elastic_f_table_error.png", "Axial stress error", "time", "sigma_x - sigma_x^expected", [("error", [r["time"] for r in rows_out], [r["error_sigma_x"] for r in rows_out])])
    write_rows(output_dir / "elastic_f_table_metrics.csv", rows_out)
    write_json(output_dir / "elastic_f_table_summary.json", {"summaries": summaries, "constrained_modulus": M, "young_modulus": E})
    tex = [r"\paragraph{Quantitative result.} The expected axial stress uses $M=E(1-\nu)/[(1+\nu)(1-2\nu)]$ before lateral release and $E$ after release."]
    tex.append(r"{\scriptsize\begin{tabular}{lrrr}\toprule Variant & samples & RMS stress error & max stress error \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("num_samples", 0), 0) + " & " + compact_float(s.get("rms_sigma_x_error")) + " & " + compact_float(s.get("max_abs_sigma_x_error")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["elastic_f_table_stress_error.png", "elastic_f_table_error.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    tex.extend(visit_tex)
    (output_dir / "elastic_f_table_results.tex").write_text("\n".join(tex) + "\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
