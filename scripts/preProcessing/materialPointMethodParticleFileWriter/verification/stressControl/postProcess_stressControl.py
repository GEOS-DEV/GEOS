#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

VARIANTS = [
    {"name": "pOnly", "label": "P controller", "case_name": "elasticStressControl_pOnly"},
    {"name": "pd", "label": "PD controller", "case_name": "elasticStressControl_pd"},
]
DEFAULT_NAME = "stressControl"
STRESS_TABLE = [(0.00, 0.0, 0.0, 0.0), (0.25, -0.01, -0.01, -0.01), (0.75, -0.01, -0.04, -0.01), (1.00, 0.0, 0.0, 0.0)]


def target(t: float, component: int) -> float:
    if t <= STRESS_TABLE[0][0]:
        return STRESS_TABLE[0][1 + component]
    for a, b in zip(STRESS_TABLE[:-1], STRESS_TABLE[1:]):
        if t <= b[0]:
            s = (t - a[0]) / max(b[0] - a[0], 1.0e-30)
            return a[1 + component] + s * (b[1 + component] - a[1 + component])
    return STRESS_TABLE[-1][1 + component]


def pick_component(row: dict, component: int) -> float | None:
    labels = [["xx", "11", "_0", "x"], ["yy", "22", "_1", "y"], ["zz", "33", "_2", "z"]][component]
    for k in row:
        low = k.lower()
        if "stress" in low and any(tok in low for tok in labels):
            try:
                return float(row[k])
            except Exception:
                pass
    return None


def main() -> int:
    args = parse_common_args("Post-process stress-control verification")
    source_dir = Path(args.source_dir); output_dir = Path(args.output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)
    rows_out, summaries = [], []
    visit_tex = []
    series_y = []
    target_y_series = None
    for subcase in manifest.get("subcases", []):
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": subcase.get("name"), "error": "run directory not found"}); continue
        box_path = find_first(run_dir, ["boxAverageHistory.csv"])
        data = read_csv_numeric(box_path) if box_path else []
        err2 = []; maxerr = 0.0; times = []; measured_y = []; target_y = []
        for row in data:
            t = float(row.get("time", row.get("t", 0.0)))
            comps = [pick_component(row, i) for i in range(3)]
            if any(v is None for v in comps):
                continue
            tgts = [target(t, i) for i in range(3)]
            errs = [float(c) - e for c, e in zip(comps, tgts)]
            err2.extend([e*e for e in errs]); maxerr = max(maxerr, max(abs(e) for e in errs))
            times.append(t); measured_y.append(float(comps[1])); target_y.append(tgts[1])
            rows_out.append({"variant": subcase.get("name"), "label": subcase.get("label"), "time": t, "sigma_x": comps[0], "sigma_y": comps[1], "sigma_z": comps[2], "target_x": tgts[0], "target_y": tgts[1], "target_z": tgts[2], "error_norm": math.sqrt(sum(e*e for e in errs))})
        if err2:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": len(times), "rms_stress_tracking_error": math.sqrt(sum(err2)/len(err2)), "max_abs_stress_tracking_error": maxerr})
            series_y.append((str(subcase.get("label")), times, measured_y)); target_y_series = ("target", times, target_y)
        else:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "num_samples": 0, "error": "stress columns not found"})
        frames = render_visit_frames(args, run_dir, output_dir, str(subcase.get("case_name")), "Stress", states="initial,middle,final", view="auto", range_mode="auto")
        visit_tex.extend(visit_frames_tex(frames, output_dir, "Particle stress frames for a representative stress-control subcase show the initial, intermediate, and final states."))
    if target_y_series:
        series_y.insert(0, target_y_series)
    plot_metric(output_dir, "stress_control_sigma_y.png", "Stress-control tracking", "time", "sigma_y", series_y)
    plot_metric(output_dir, "stress_control_error_norm.png", "Stress-control error norm", "time", "||sigma - target||", [("error", [r["time"] for r in rows_out], [r["error_norm"] for r in rows_out])])
    write_rows(output_dir / "stress_control_metrics.csv", rows_out)
    write_json(output_dir / "stress_control_summary.json", {"summaries": summaries, "stress_table": STRESS_TABLE})
    tex = [r"\paragraph{Quantitative result.} The stress target is the prescribed piecewise table; the metric is the RMS component-wise tracking error."]
    tex.append(r"{\scriptsize\begin{tabular}{lrrr}\toprule Variant & samples & RMS tracking error & max tracking error \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("num_samples", 0), 0) + " & " + compact_float(s.get("rms_stress_tracking_error")) + " & " + compact_float(s.get("max_abs_stress_tracking_error")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["stress_control_sigma_y.png", "stress_control_error_norm.png"]:
        if (output_dir / name).is_file(): tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    tex.extend(visit_tex)
    (output_dir / "stress_control_results.tex").write_text("\n".join(tex) + "\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
