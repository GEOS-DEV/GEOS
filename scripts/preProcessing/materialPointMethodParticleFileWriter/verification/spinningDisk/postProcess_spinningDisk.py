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
    {"name": "pic", "label": "PIC", "case_name": "spinningDiskConservation_pic"},
    {"name": "flip", "label": "FLIP", "case_name": "spinningDiskConservation_flip"},
    {"name": "fmpm2", "label": "FMPM2", "case_name": "spinningDiskConservation_fmpm2"},
]
DEFAULT_NAME = "spinningDisk"


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


def safe_relative(value: float, reference: float) -> float:
    if abs(reference) < 1.0e-30:
        return value - reference
    return value / reference - 1.0


def main() -> int:
    args = parse_common_args("Post-process spinning-disk conservation verification")
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(source_dir, output_dir, args.case_id, DEFAULT_NAME, VARIANTS)

    rows_out: list[dict] = []
    summaries: list[dict] = []
    angular_series: list[tuple[str, list[float], list[float]]] = []
    energy_series: list[tuple[str, list[float], list[float]]] = []

    for subcase in manifest.get("subcases", []):
        run_dir = run_dir_from_subcase(source_dir, subcase)
        if run_dir is None:
            summaries.append({"variant": subcase.get("name"), "error": "run directory not found"})
            continue

        by_time: dict[float, list[tuple[float, float, float, float]]] = {}
        for tracer_file in tracer_files(run_dir, "spinningDiskTracer"):
            try:
                rows = read_tracer(tracer_file)
            except Exception as exc:
                rows_out.append({"variant": subcase.get("name"), "tracer": tracer_file.name, "error": str(exc)})
                continue
            for row in rows:
                vx = tracer_velocity(row, "x")
                vy = tracer_velocity(row, "y")
                if vx is None or vy is None:
                    continue
                by_time.setdefault(float(row["t"]), []).append((float(row["x"]), float(row["y"]), vx, vy))

        times = sorted(by_time)
        angular: list[float] = []
        energy: list[float] = []
        for t in times:
            samples = by_time[t]
            lz = sum(x * vy - y * vx for x, y, vx, vy in samples)
            ke = sum(0.5 * (vx * vx + vy * vy) for _x, _y, vx, vy in samples)
            angular.append(lz)
            energy.append(ke)

        if times:
            l0 = angular[0]
            e0 = energy[0]
            l_rel = [safe_relative(v, l0) for v in angular]
            e_rel = [safe_relative(v, e0) for v in energy]
            for t, lz, ke, dl, de in zip(times, angular, energy, l_rel, e_rel):
                rows_out.append({
                    "variant": subcase.get("name"),
                    "label": subcase.get("label"),
                    "time": t,
                    "tracer_count": len(by_time[t]),
                    "angular_momentum_proxy": lz,
                    "kinetic_energy_proxy": ke,
                    "relative_angular_momentum_drift": dl,
                    "relative_kinetic_energy_drift": de,
                })
            summaries.append({
                "variant": subcase.get("name"),
                "label": subcase.get("label"),
                "num_time_samples": len(times),
                "tracer_count_initial": len(by_time[times[0]]),
                "max_abs_relative_angular_momentum_drift": max(abs(v) for v in l_rel),
                "final_relative_angular_momentum_drift": l_rel[-1],
                "max_abs_relative_kinetic_energy_drift": max(abs(v) for v in e_rel),
                "final_relative_kinetic_energy_drift": e_rel[-1],
            })
            angular_series.append((str(subcase.get("label")), times, l_rel))
            energy_series.append((str(subcase.get("label")), times, e_rel))
        else:
            summaries.append({"variant": subcase.get("name"), "label": subcase.get("label"), "num_time_samples": 0, "error": "no tracer velocity data found"})

        render_visit_frames(args, run_dir, output_dir, str(subcase.get("case_name")), "Velocity", states="initial,middle,final", view="auto", colortable="bluehot", range_mode="auto")

    plot_metric(output_dir, "spinning_disk_angular_momentum_drift.png", "Spinning disk angular-momentum drift", "time", "relative drift", angular_series)
    plot_metric(output_dir, "spinning_disk_kinetic_energy_drift.png", "Spinning disk kinetic-energy drift", "time", "relative drift", energy_series)
    write_rows(output_dir / "spinning_disk_metrics.csv", rows_out)
    write_json(output_dir / "spinning_disk_summary.json", {"summaries": summaries})

    tex = [r"\paragraph{Quantitative result.} Tracer velocities define proxy angular momentum $L_z=\sum_i(x_iv_{y,i}-y_iv_{x,i})$ and kinetic energy $K=\sum_i|\mathbf{v}_i|^2/2$.  The ideal force-free solution has constant $L_z$ and $K$; PIC is expected to be dissipative, while FLIP and FMPM2 should reduce drift."]
    tex.append(r"{\scriptsize\begin{tabular}{lrrrr}\toprule Variant & samples & max $|\Delta L_z/L_z|$ & final $\Delta L_z/L_z$ & max $|\Delta K/K|$ \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("num_time_samples", 0), 0) + " & " + compact_float(s.get("max_abs_relative_angular_momentum_drift")) + " & " + compact_float(s.get("final_relative_angular_momentum_drift")) + " & " + compact_float(s.get("max_abs_relative_kinetic_energy_drift")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["spinning_disk_angular_momentum_drift.png", "spinning_disk_kinetic_energy_drift.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.48\linewidth]{\CaseOutputDir/" + name + "}")
    (output_dir / "spinning_disk_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
