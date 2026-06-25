#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

SOURCE_DIR = Path(__file__).resolve().parent
PFW_ROOT = SOURCE_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))

from mpm_vv_simple_post import compact_float, latex_escape, read_csv_numeric, write_json, write_rows, plot_metric

DEFAULT_CASE_ID = "weakTraceRotatedBar"
VARIANT_LABELS = {
    "singleField": "single-field baseline",
    "falseElasticCZ": "false-elastic CZ reference",
    "traceContactGroups": "trace projection, contact groups",
    "dfgSurfaceFlags": "DFG surface flags placeholder",
    "stressResidualPartition": "stress-residual partition placeholder",
}
LAYER_ORDER = ["compliant_near", "stiff_near", "compliant_far", "stiff_far"]
TRACERS_PER_LAYER = 3
N = (1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0), 0.0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Post-process weakTraceRotatedBar verification")
    p.add_argument("--source-dir", default=str(SOURCE_DIR))
    p.add_argument("--output-dir", default=str(SOURCE_DIR / "output" / DEFAULT_CASE_ID))
    p.add_argument("--case-id", default=DEFAULT_CASE_ID)
    p.add_argument("--output-prefix", default="weakTraceRotatedBar")
    p.add_argument("--run-dir", default="")
    p.add_argument("--python", dest="python_cmd", default=sys.executable)
    p.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", ""))
    p.add_argument("--visit-timeout", type=float, default=float(os.environ.get("MPM_VV_VISIT_TIMEOUT", "90")))
    p.add_argument("--no-visit", action="store_true")
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def find_manifest(source_dir: Path, output_dir: Path, case_id: str) -> dict:
    for path in [output_dir / "weakTraceRotatedBar_jobs.json", output_dir / f"{case_id}_jobs.json"]:
        if path.is_file():
            return json.loads(path.read_text())
    subcases = []
    for name, label in VARIANT_LABELS.items():
        jobs_path = source_dir / "output" / f"{case_id}__{name}" / f"weakTraceRotatedBar_{name}_jobs.json"
        jobs = {}
        if jobs_path.is_file():
            try:
                jobs = json.loads(jobs_path.read_text())
            except Exception:
                jobs = {}
        subcases.append({"name": name, "label": label, "case_id": f"{case_id}__{name}", "case_name": f"weakTraceRotatedBar_{name}", "jobs": jobs, "placeholder": name not in {"singleField", "falseElasticCZ", "traceContactGroups"}})
    return {"case_id": case_id, "subcases": subcases}


def run_dir_from_subcase(source_dir: Path, subcase: dict) -> Path | None:
    run_dir = str(subcase.get("jobs", {}).get("run_dir", "")).strip()
    if run_dir and Path(run_dir).is_dir():
        return Path(run_dir)
    for p in [source_dir / "output" / subcase.get("case_id", ""), source_dir / "output" / subcase.get("case_id", "") / subcase.get("case_name", "")]:
        if p.is_dir():
            return p
    return None


def tracer_files(run_dir: Path) -> list[Path]:
    return sorted(set(run_dir.glob("weakTraceRotatedBarTracer_*.csv")) | set(run_dir.glob("tracerPoint_*.csv")))


def parse_tracer(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            out = {}
            for k, v in row.items():
                if v in (None, ""):
                    continue
                try:
                    out[k] = float(v)
                except Exception:
                    out[k] = v
            rows.append(out)
    return rows


def stress_from_row(row: dict) -> tuple[float, float, float, float, float, float] | None:
    keys = ["stressXX", "stressYY", "stressZZ", "stressYZ", "stressXZ", "stressXY"]
    if not all(k in row for k in keys):
        return None
    return tuple(float(row[k]) for k in keys)


def von_mises(sig) -> float:
    sxx, syy, szz, syz, sxz, sxy = sig
    return math.sqrt(max(0.0, 0.5 * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2) + 3.0 * (sxy*sxy + sxz*sxz + syz*syz)))


def normal_stress(sig) -> float:
    sxx, syy, szz, syz, sxz, sxy = sig
    nx, ny, nz = N
    return nx*nx*sxx + ny*ny*syy + nz*nz*szz + 2.0*(nx*ny*sxy + nx*nz*sxz + ny*nz*syz)


def tracer_layer(index: int) -> str:
    return LAYER_ORDER[min(len(LAYER_ORDER) - 1, index // TRACERS_PER_LAYER)]


def process_variant(source_dir: Path, subcase: dict) -> tuple[list[dict], list[dict], dict]:
    if subcase.get("placeholder"):
        return [], [], {"variant": subcase.get("name"), "label": subcase.get("label"), "status": "placeholder"}
    run_dir = run_dir_from_subcase(source_dir, subcase)
    if run_dir is None:
        return [], [], {"variant": subcase.get("name"), "label": subcase.get("label"), "status": "missing run directory"}
    files = tracer_files(run_dir)
    if not files:
        return [], [], {"variant": subcase.get("name"), "label": subcase.get("label"), "status": "missing tracer files", "run_dir": str(run_dir)}

    rows = []
    for idx, path in enumerate(files):
        layer = tracer_layer(idx)
        for row in parse_tracer(path):
            sig = stress_from_row(row)
            if sig is None:
                continue
            t = float(row.get("t", row.get("time", 0.0)))
            rows.append({
                "variant": subcase.get("name"),
                "label": subcase.get("label"),
                "tracer": path.name,
                "tracerIndex": idx,
                "layer": layer,
                "time": t,
                "sigmaNN": normal_stress(sig),
                "vonMisesStress": von_mises(sig),
                "plasticStrainMagnitude": float(row.get("plasticStrainMagnitude", 0.0) or 0.0),
                "materialType": row.get("materialType", ""),
                "density": row.get("density", ""),
            })

    by = defaultdict(list)
    for r in rows:
        key = (r["variant"], r["label"], r["layer"], r["time"])
        by[key].append(r)
    layer_rows = []
    for (variant, label, layer, t), vals in sorted(by.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][3])):
        n = len(vals)
        layer_rows.append({
            "variant": variant,
            "label": label,
            "layer": layer,
            "time": t,
            "numTracers": n,
            "meanSigmaNN": sum(v["sigmaNN"] for v in vals) / n,
            "meanVonMisesStress": sum(v["vonMisesStress"] for v in vals) / n,
            "meanPlasticStrainMagnitude": sum(v["plasticStrainMagnitude"] for v in vals) / n,
        })
    final_by_layer = {r["layer"]: r for r in layer_rows if r["time"] == max((x["time"] for x in layer_rows), default=0.0)}
    summary = {
        "variant": subcase.get("name"),
        "label": subcase.get("label"),
        "status": "processed",
        "run_dir": str(run_dir),
        "numTracerRows": len(rows),
    }
    for layer in LAYER_ORDER:
        if layer in final_by_layer:
            summary[f"final_{layer}_vm"] = final_by_layer[layer]["meanVonMisesStress"]
            summary[f"final_{layer}_sigmaNN"] = final_by_layer[layer]["meanSigmaNN"]
            summary[f"final_{layer}_plastic"] = final_by_layer[layer]["meanPlasticStrainMagnitude"]
    return rows, layer_rows, summary


def plot_layers(output_dir: Path, layer_rows: list[dict], column: str, filename: str, ylabel: str):
    series = []
    for variant in sorted({r["variant"] for r in layer_rows}):
        for layer in LAYER_ORDER:
            data = sorted([r for r in layer_rows if r["variant"] == variant and r["layer"] == layer], key=lambda r: r["time"])
            if data:
                series.append((f"{variant} / {layer}", [r["time"] for r in data], [r[column] for r in data]))
    plot_metric(output_dir, filename, filename.replace("_", " ").replace(".png", ""), "time", ylabel, series)


def render_visit(args, subcases: list[dict], output_dir: Path):
    if args.no_visit:
        return []
    visit = args.visit_cmd or os.environ.get("VISIT_COMMAND", "") or shutil.which("visit")
    if not visit:
        return []
    frames = []
    for subcase in subcases:
        if subcase.get("placeholder"):
            continue
        run_dir = run_dir_from_subcase(Path(args.source_dir), subcase)
        if run_dir is None or not (run_dir / "siloFiles").is_dir():
            continue
        out = output_dir / f"weak_trace_rotated_bar_{subcase.get('name')}_vm_final.png"
        cmd = [visit, "-nowin", "-cli", "-s", str(run_dir / "visitRender_weakTraceRotatedBar.py"), "--run-dir", str(run_dir), "--output", str(out), "--state", "final"]
        try:
            subprocess.run(cmd, cwd=run_dir, timeout=args.visit_timeout, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if out.is_file():
                frames.append(out)
        except Exception:
            pass
    return frames


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = find_manifest(source_dir, output_dir, args.case_id)

    all_rows = []
    all_layers = []
    summaries = []
    for subcase in manifest.get("subcases", []):
        rows, layers, summary = process_variant(source_dir, subcase)
        all_rows.extend(rows)
        all_layers.extend(layers)
        summaries.append(summary)

    write_rows(output_dir / "weak_trace_rotated_bar_tracer_metrics.csv", all_rows)
    write_rows(output_dir / "weak_trace_rotated_bar_layer_metrics.csv", all_layers)
    final_layers = []
    for variant in sorted({r["variant"] for r in all_layers}):
        max_t = max((r["time"] for r in all_layers if r["variant"] == variant), default=None)
        if max_t is not None:
            final_layers.extend([r for r in all_layers if r["variant"] == variant and r["time"] == max_t])
    write_rows(output_dir / "weak_trace_rotated_bar_final_layer_metrics.csv", final_layers)
    write_json(output_dir / "weak_trace_rotated_bar_summary.json", {"summaries": summaries})

    plot_layers(output_dir, all_layers, "meanSigmaNN", "weak_trace_rotated_bar_layer_normal_stress.png", "layer mean normal stress")
    plot_layers(output_dir, all_layers, "meanVonMisesStress", "weak_trace_rotated_bar_layer_vm_stress.png", "layer mean von Mises stress")
    plot_layers(output_dir, all_layers, "meanPlasticStrainMagnitude", "weak_trace_rotated_bar_layer_plastic.png", "layer mean plastic strain")
    frames = render_visit(args, manifest.get("subcases", []), output_dir)

    tex = []
    tex.append(r"\paragraph{Local stress diagnostics.} The rotated-bar verification compares local tracer stresses on both sides of the inclined interface.  The primary quantities are the material-normal stress $\sigma_{nn}=\mathbf{n}\cdot\boldsymbol{\sigma}\mathbf{n}$ and von Mises stress reconstructed from tracer stress components.  The stiff phase is deliberately weak so that a mixed-cell stress spike can cause spurious plasticity even though the correct series response is elastic.")
    tex.append(r"{\scriptsize\begin{tabular}{lrrrr}\toprule Variant & stiff near VM & stiff far VM & stiff near plastic & status \\\midrule")
    for s in summaries:
        tex.append(latex_escape(s.get("label", s.get("variant", ""))) + " & " + compact_float(s.get("final_stiff_near_vm")) + " & " + compact_float(s.get("final_stiff_far_vm")) + " & " + compact_float(s.get("final_stiff_near_plastic")) + " & " + latex_escape(s.get("status", "")) + r" \\")
    tex.append(r"\bottomrule\end{tabular}}")
    for name in ["weak_trace_rotated_bar_layer_normal_stress.png", "weak_trace_rotated_bar_layer_vm_stress.png", "weak_trace_rotated_bar_layer_plastic.png"]:
        if (output_dir / name).is_file():
            tex.append(r"\includegraphics[width=0.49\linewidth]{\CaseOutputDir/" + name + "}")
    for frame in frames:
        try:
            rel = frame.relative_to(output_dir)
        except ValueError:
            rel = frame.name
        tex.append(r"\includegraphics[width=0.49\linewidth]{\CaseOutputDir/" + str(rel).replace(os.sep, "/") + "}")
    (output_dir / "weak_trace_rotated_bar_results.tex").write_text("\n".join(tex) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
