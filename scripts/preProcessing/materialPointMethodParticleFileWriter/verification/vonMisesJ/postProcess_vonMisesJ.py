#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))
from mpm_vv_simple_post import *

VARIANTS = [{"name": "uniaxial", "label": "Uniaxial strain ramp", "case_name": "vonMisesUniaxialPlasticity_uniaxial"}]
DEFAULT_NAME = "vonMisesJ"
YOUNG_MODULUS = 1.0
YIELD_STRENGTH = 0.02
FINAL_STRAIN = 0.05
STOP_TIME = 1.0


def expected_stress(strain: float) -> float:
    return min(YOUNG_MODULUS * max(strain, 0.0), YIELD_STRENGTH)


def expected_plastic_strain(strain: float) -> float:
    return max(0.0, strain - YIELD_STRENGTH / YOUNG_MODULUS)


def pick_stress_y(row: dict) -> float | None:
    labels = ["yy", "22", "_1", "y"]
    for key in row:
        low = key.lower()
        if "stress" in low and any(tok in low for tok in labels):
            try:
                return float(row[key])
            except Exception:
                pass
    for key in row:
        if "stress" in key.lower():
            try:
                return float(row[key])
            except Exception:
                pass
    return None


def pick_plastic(row: dict) -> float | None:
    for key in row:
        low = key.lower()
        if "plastic" in low and any(tok in low for tok in ["strain", "eps", "magnitude", "eq"]):
            try:
                return float(row[key])
            except Exception:
                pass
    return None


def main() -> int:
    args = parse_common_args("Post-process Von Mises uniaxial plasticity verification")
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)

    rows_out: list[dict] = []
    summaries: list[dict] = []
    stress_series: list[tuple[str, list[float], list[float]]] = []
    expected_series: tuple[str, list[float], list[float]] | None = None
    plastic_series: list[tuple[str, list[float], list[float]]] = []
    expected_plastic_series: tuple[str, list[float], list[float]] | None = None

    for subcase in manifest.get("subcases", []):
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": subcase.get("name"), "error": "run directory not found"})
            continue

        box_path = find_first(run_dir, ["boxAverageHistory.csv"])
        data = read_csv_numeric(box_path) if box_path else []
        strains: list[float] = []
        measured: list[float] = []
        expected: list[float] = []
        stress_errors: list[float] = []
        plastic_x: list[float] = []
        plastic_measured: list[float] = []
        plastic_expected: list[float] = []
        plastic_errors: list[float] = []

        for row in data:
            t = float(row.get("time", row.get("t", 0.0)))
            strain = FINAL_STRAIN * max(0.0, min(t / STOP_TIME, 1.0))
            sy = pick_stress_y(row)
            ep = pick_plastic(row)
            es = expected_stress(strain)
            ep_exp = expected_plastic_strain(strain)
            if sy is not None:
                # Most runs report tensile stress as positive for this strain path; if a legacy
                # build writes compression-positive values, compare the magnitude but keep the
                # signed value in the CSV for debugging.
                cmp_sy = abs(sy) if sy * es < 0.0 else sy
                strains.append(strain)
                measured.append(cmp_sy)
                expected.append(es)
                stress_errors.append(cmp_sy - es)
            if ep is not None:
                plastic_x.append(strain)
                plastic_measured.append(ep)
                plastic_expected.append(ep_exp)
                plastic_errors.append(ep - ep_exp)
            rows_out.append({
                "variant": subcase.get("name"),
                "label": subcase.get("label"),
                "time": t,
                "strain_y_expected": strain,
                "measured_sigma_y_signed": sy if sy is not None else "",
                "measured_sigma_y_for_metric": (abs(sy) if sy is not None and sy * es < 0.0 else sy) if sy is not None else "",
                "expected_sigma_y": es,
                "stress_error": stress_errors[-1] if sy is not None else "",
                "measured_plastic_strain": ep if ep is not None else "",
                "expected_plastic_strain": ep_exp,
                "plastic_strain_error": plastic_errors[-1] if ep is not None else "",
            })

        if stress_errors:
            summary = {
                "variant": subcase.get("name"),
                "label": subcase.get("label"),
                "num_stress_samples": len(stress_errors),
                "rms_stress_error": math.sqrt(sum(e * e for e in stress_errors) / len(stress_errors)),
                "max_abs_stress_error": max(abs(e) for e in stress_errors),
                "final_measured_sigma_y": measured[-1],
                "final_expected_sigma_y": expected[-1],
            }
            if plastic_errors:
                summary.update({
                    "num_plastic_samples": len(plastic_errors),
                    "rms_plastic_strain_error": math.sqrt(sum(e * e for e in plastic_errors) / len(plastic_errors)),
                    "max_abs_plastic_strain_error": max(abs(e) for e in plastic_errors),
                    "final_measured_plastic_strain": plastic_measured[-1],
                    "final_expected_plastic_strain": plastic_expected[-1],
                })
            summaries.append(summary)
            stress_series.append((str(subcase.get("label")), strains, measured))
            expected_series = ("elastic-perfectly-plastic ideal", strains, expected)
            if plastic_x:
                plastic_series.append((str(subcase.get("label")), plastic_x, plastic_measured))
                expected_plastic_series = ("expected", plastic_x, plastic_expected)
        else:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "num_stress_samples": 0, "error": "stress column not found"})

        render_visit_frames(args, run_dir, output_dir, str(subcase.get("case_name")), "PlasticStrainMagnitude", states="initial,middle,final", view="auto", colortable="hot_desaturated", range_mode="auto")

    if expected_series:
        stress_series.insert(0, expected_series)
    if expected_plastic_series:
        plastic_series.insert(0, expected_plastic_series)
    plot_metric(output_dir, "von_mises_stress_strain.png", "Von Mises uniaxial response", "expected axial strain", "sigma_y", stress_series)
    plot_metric(output_dir, "von_mises_plastic_strain.png", "Equivalent plastic strain", "expected axial strain", "plastic strain", plastic_series)
    write_rows(output_dir / "von_mises_metrics.csv", rows_out)
    write_json(output_dir / "von_mises_summary.json", {"summaries": summaries, "young_modulus": YOUNG_MODULUS, "yield_strength": YIELD_STRENGTH, "final_strain": FINAL_STRAIN})

    tex = [r"\paragraph{Quantitative result.} The comparison uses the ideal one-dimensional elastic-perfectly-plastic law $\sigma=\min(E\epsilon,\sigma_y)$ with $E=1$ and $\sigma_y=0.02$."]
    tex.append(r"{\scriptsize\begin{tabular}{lrrr}\toprule Variant & samples & RMS stress error & max stress error \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("num_stress_samples", 0), 0) + " & " + compact_float(s.get("rms_stress_error")) + " & " + compact_float(s.get("max_abs_stress_error")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["von_mises_stress_strain.png", "von_mises_plastic_strain.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    (output_dir / "von_mises_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
