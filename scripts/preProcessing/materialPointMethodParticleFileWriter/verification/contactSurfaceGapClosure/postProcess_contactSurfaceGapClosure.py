#!/usr/bin/env python3
"""Aggregate post-processing for the contact surface/gap-closure verification test.

The script reduces each subcase reaction-history CSV to a common table, compares
measured y-reaction stress with the ideal zero-then-uniaxial-strain response,
optionally renders VisIt stress and surface-field diagnostics, and writes the LaTeX fragment
included by test.tex.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

FIELD_MODES = [
    ("multiField", "Multi-field"),
    ("dfg", "DFG"),
]
SURFACE_MODES = [
    ("implicit", "implicit surfaces"),
    ("explicit", "explicit surfaces"),
]
GAP_CORRECTIONS = ["Simple", "Implicit", "Softened"]
VARIANTS = []
for field_name, field_label in FIELD_MODES:
    for surface_name, surface_label in SURFACE_MODES:
        for gap_name in GAP_CORRECTIONS:
            name = f"{field_name}_{surface_name}_{gap_name}"
            VARIANTS.append(
                {
                    "name": name,
                    "label": f"{field_label} / {surface_label} / {gap_name} gap",
                    "case_name": f"contactGapClosure_{name}",
                    "field_mode": field_name,
                    "surface_mode": surface_name,
                    "gap_correction": gap_name,
                }
            )
DEFAULT_VARIANT_NAMES = [
    "multiField_explicit_Softened",
    "dfg_implicit_Implicit",
    "dfg_explicit_Softened",
]


DEFAULT_PARAMETERS = {
    "domain_width": 1.0,
    "domain_height": 1.0,
    "initial_gap": 0.10,
    "extra_compression": 0.12,
    "notch_radius": 0.18,
    "refine": 1,
    "cpp": 32,
    "ppc": 2,
    "cells_per_partition": 32,
    "young_modulus": 10.0,
    "poisson_ratio": 0.25,
    "density": 1.0,
}
REACTION_CSV = "reactionHistory.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Post-process contact surface/gap-closure variants")
    parser.add_argument("--suite", default="verification")
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--case-id", default="contactSurfaceGapClosure")
    parser.add_argument("--python", dest="python_cmd", default=sys.executable)
    parser.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", ""))
    parser.add_argument("--no-visit", action="store_true", help="Do not launch VisIt; still include any existing valid PNGs in the LaTeX fragment")
    parser.add_argument("--rerender-visit", action="store_true", help="Re-render VisIt frames even when matching PNGs already exist")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def compact_float(value: float | None) -> str:
    if value is None or not math.isfinite(float(value)):
        return "--"
    return f"{float(value):.6g}"


def latex_escape(value: object) -> str:
    text = str(value if value is not None else "")
    repl = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(repl.get(ch, ch) for ch in text)


def variant_case_id(folder_case_id: str, variant: dict[str, str]) -> str:
    return f"{folder_case_id}__{variant['name']}"


def load_manifest(source_dir: Path, output_dir: Path, case_id: str) -> dict:
    manifest_path = output_dir / "contactSurfaceGapClosure_jobs.json"
    if manifest_path.is_file():
        try:
            return json.loads(manifest_path.read_text())
        except Exception:
            pass
    subcases = []
    fallback_variants = [v for v in VARIANTS if v["name"] in set(DEFAULT_VARIANT_NAMES)]
    for variant in fallback_variants:
        cid = variant_case_id(case_id, variant)
        jobs_path = source_dir / "output" / cid / f"{variant['case_name']}_jobs.json"
        jobs = {}
        if jobs_path.is_file():
            try:
                jobs = json.loads(jobs_path.read_text())
            except Exception:
                jobs = {}
        subcases.append({**variant, "case_id": cid, "jobs": jobs})
    return {"case_id": case_id, "parameters": dict(DEFAULT_PARAMETERS), "subcases": subcases}


def run_dir_from_subcase(source_dir: Path, subcase: dict) -> Path | None:
    jobs = subcase.get("jobs", {}) if isinstance(subcase.get("jobs", {}), dict) else {}
    run_dir = jobs.get("run_dir")
    if run_dir:
        path = Path(run_dir)
        if path.is_dir():
            return path
    case_id = subcase.get("case_id", "")
    case_name = subcase.get("case_name", "")
    candidates = [
        source_dir / "output" / case_id / case_name,
        source_dir / "output" / case_id,
    ]
    for path in candidates:
        if path.is_dir():
            return path
    return None


def find_reaction_csv(run_dir: Path) -> Path | None:
    direct = run_dir / REACTION_CSV
    if direct.is_file():
        return direct
    matches = sorted(run_dir.glob(f"**/{REACTION_CSV}"))
    return matches[0] if matches else None


def parse_float(text: object, default: float = 0.0) -> float:
    try:
        return float(str(text).strip())
    except Exception:
        return default


def read_reaction_csv(path: Path, params: dict) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return rows
        # reactionHistory.csv historically uses spaces after commas. Normalize keys.
        field_map = {name: name.strip() for name in reader.fieldnames}
        for raw in reader:
            row = {field_map[key]: parse_float(value, 0.0) for key, value in raw.items()}
            if not row:
                continue
            time = row.get("time", row.get("Time", 0.0))
            fyy = row.get("F11", 1.0)
            length_x = row.get("length_x", params["domain_width"])
            length_z = row.get("length_z", 1.0)
            area_y = max(abs(length_x * length_z), 1.0e-30)
            ry_minus = row.get("Ry-", 0.0)
            ry_plus = row.get("Ry+", 0.0)
            # Positive compression convention for the plotted verification metric.
            stress_y = 0.5 * (abs(ry_minus) + abs(ry_plus)) / area_y
            displacement = params["domain_height"] * (1.0 - fyy)
            ideal = ideal_reaction_stress(displacement, params)
            rows.append(
                {
                    "time": time,
                    "Fyy": fyy,
                    "displacement": displacement,
                    "reaction_stress_y": stress_y,
                    "ideal_reaction_stress_y": ideal,
                    "reaction_error_y": stress_y - ideal,
                    "Ry_minus": ry_minus,
                    "Ry_plus": ry_plus,
                    "area_y": area_y,
                }
            )
    # Drop duplicate/restart rows while preserving monotonic time.
    out: list[dict[str, float]] = []
    max_time = -math.inf
    for row in rows:
        if row["time"] > max_time:
            out.append(row)
            max_time = row["time"]
    return out


def ideal_reaction_stress(displacement: float, params: dict) -> float:
    gap = float(params["initial_gap"])
    if displacement <= gap:
        return 0.0
    solid_height = max(float(params["domain_height"]) - gap, 1.0e-30)
    compression = max(displacement - gap, 0.0)
    lambda_y = max(1.0 - compression / solid_height, 1.0e-8)
    young = float(params["young_modulus"])
    nu = float(params["poisson_ratio"])
    lame_lambda = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    shear = young / (2.0 * (1.0 + nu))
    sigma_yy = lame_lambda * math.log(lambda_y) / lambda_y + shear * (lambda_y * lambda_y - 1.0) / lambda_y
    return max(-sigma_yy, 0.0)


def summarize(history: list[dict[str, float]], params: dict) -> dict[str, float]:
    if not history:
        return {"num_rows": 0}
    gap = params["initial_gap"]
    post = [row for row in history if row["displacement"] > gap]
    pre = [row for row in history if row["displacement"] <= gap]
    final = history[-1]
    errors = [row["reaction_error_y"] for row in post] if post else []
    ideal_scale = max(max((row["ideal_reaction_stress_y"] for row in history), default=0.0), 1.0e-30)
    rmse = math.sqrt(sum(e * e for e in errors) / len(errors)) if errors else 0.0
    return {
        "num_rows": len(history),
        "final_time": final["time"],
        "final_displacement": final["displacement"],
        "final_reaction_stress_y": final["reaction_stress_y"],
        "final_ideal_reaction_stress_y": final["ideal_reaction_stress_y"],
        "final_reaction_error_y": final["reaction_error_y"],
        "max_abs_preclosure_reaction_stress_y": max((abs(row["reaction_stress_y"]) for row in pre), default=0.0),
        "postclosure_rmse_reaction_stress_y": rmse,
        "postclosure_relative_rmse_reaction_stress_y": rmse / ideal_scale,
        "max_abs_reaction_error_y": max((abs(row["reaction_error_y"]) for row in history), default=0.0),
    }


def params_from_manifest(manifest: dict) -> dict:
    params = dict(DEFAULT_PARAMETERS)
    params.update({k: v for k, v in manifest.get("parameters", {}).items() if k in params})
    # Mirror the pfw_materials.contactGapClosureHyperelastic entry.  Keep these
    # local so the post-processor can run on an output bundle without importing PFW.
    params.setdefault("young_modulus", 10.0)
    params.setdefault("poisson_ratio", 0.25)
    params.setdefault("density", 1.0)
    return params


def collect_results(source_dir: Path, output_dir: Path, case_id: str) -> tuple[list[dict], dict]:
    manifest = load_manifest(source_dir, output_dir, case_id)
    params = params_from_manifest(manifest)
    results = []
    for subcase in manifest.get("subcases", []):
        name = subcase.get("name", "")
        variant = next((v for v in VARIANTS if v["name"] == name), None)
        if variant is None:
            continue
        run_dir = run_dir_from_subcase(source_dir, subcase)
        result = {**variant, "run_dir": str(run_dir) if run_dir else "", "history": [], "summary": {}}
        if run_dir is None:
            result["error"] = "run directory not found"
            results.append(result)
            continue
        csv_path = find_reaction_csv(run_dir)
        if not csv_path:
            result["error"] = f"{REACTION_CSV} not found"
            results.append(result)
            continue
        try:
            history = read_reaction_csv(csv_path, params)
            result["history"] = history
            result["summary"] = summarize(history, params)
            result["reaction_csv"] = str(csv_path)
        except Exception as exc:
            result["error"] = str(exc)
        results.append(result)
    return results, params


def write_combined_csv(output_dir: Path, results: list[dict]) -> None:
    fieldnames = [
        "variant",
        "label",
        "field_mode",
        "surface_mode",
        "gap_correction",
        "time",
        "Fyy",
        "displacement",
        "reaction_stress_y",
        "ideal_reaction_stress_y",
        "reaction_error_y",
        "Ry_minus",
        "Ry_plus",
        "area_y",
    ]
    with (output_dir / "contact_surface_gap_closure_reaction_history.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for result in results:
            for row in result.get("history", []):
                out = {
                    "variant": result["name"],
                    "label": result["label"],
                    "field_mode": result["field_mode"],
                    "surface_mode": result["surface_mode"],
                    "gap_correction": result["gap_correction"],
                }
                out.update(row)
                writer.writerow(out)


def make_plots(output_dir: Path, results: list[dict], params: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Reaction stress vs displacement, grouped by field/surface mode with gap correction style in labels.
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    any_data = False
    max_disp = params["initial_gap"] + params["extra_compression"]
    disp_grid = [max_disp * i / 300.0 for i in range(301)]
    ideal_grid = [ideal_reaction_stress(d, params) for d in disp_grid]
    ax.plot(disp_grid, ideal_grid, linewidth=2.0, label="ideal zero-then-uniaxial-strain response")
    for result in results:
        history = result.get("history", [])
        if not history:
            continue
        any_data = True
        disp = [row["displacement"] for row in history]
        stress = [row["reaction_stress_y"] for row in history]
        ax.plot(disp, stress, linewidth=1.0, label=result["label"])
    if any_data:
        ax.axvline(params["initial_gap"], linewidth=0.8, linestyle="--")
        ax.set_xlabel("imposed y-displacement")
        ax.set_ylabel("compressive y-reaction stress")
        ax.set_title("Curved contact gap closure reaction")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        fig.savefig(output_dir / "contact_gap_reaction_stress.png", dpi=180)
    plt.close(fig)

    # Relative error after closure.
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    any_data = False
    for result in results:
        history = [row for row in result.get("history", []) if row["displacement"] > params["initial_gap"]]
        if not history:
            continue
        any_data = True
        disp = [row["displacement"] for row in history]
        err = [row["reaction_error_y"] for row in history]
        ax.plot(disp, err, linewidth=1.0, label=result["label"])
    if any_data:
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("imposed y-displacement")
        ax.set_ylabel("reaction stress error")
        ax.set_title("Post-closure reaction-stress error")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        fig.savefig(output_dir / "contact_gap_reaction_error.png", dpi=180)
    plt.close(fig)


def renderer_path(source_dir: Path) -> Path:
    for parent in [source_dir] + list(source_dir.parents):
        candidate = parent / "pfw_visit_render.py"
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("could not find pfw_visit_render.py")


def find_visit_command(explicit: str) -> str:
    if explicit:
        return explicit
    return shutil.which("visit") or ""


def frame_state_label(path: Path) -> str:
    stem = path.stem.lower()
    for label in ("initial", "pct33", "pct66", "final", "middle"):
        if f"_{label}_" in stem or stem.endswith(f"_{label}"):
            return label
    return ""


def representative_result(results: list[dict]) -> dict | None:
    preferred_order = [
        "dfg_explicit_Softened",
        "dfg_explicit_Implicit",
        "multiField_explicit_Softened",
        "multiField_explicit_Implicit",
    ]
    by_name = {r["name"]: r for r in results if r.get("run_dir")}
    for name in preferred_order:
        if name in by_name:
            return by_name[name]
    return next((r for r in results if r.get("run_dir")), None)


def existing_visit_frames(output_dir: Path, variant_name: str, tokens: tuple[str, ...]) -> dict[str, Path]:
    frame_dir = output_dir / "visit_frames" / variant_name
    out: dict[str, Path] = {}
    if not frame_dir.is_dir():
        return out
    compact_tokens = tuple(re.sub(r"[^a-z0-9]+", "", token.lower()) for token in tokens)
    for png in sorted(frame_dir.rglob("*.png")):
        low = png.name.lower()
        compact_name = re.sub(r"[^a-z0-9]+", "", low)
        if variant_name.lower() not in low:
            continue
        if compact_tokens and not any(token in compact_name for token in compact_tokens):
            continue
        state = frame_state_label(png)
        if state and state not in out:
            out[state] = png
    return out


def existing_vm_frames(output_dir: Path, variant_name: str) -> dict[str, Path]:
    return existing_visit_frames(output_dir, variant_name, ("vonMises", "particleVonMisesStress"))


def existing_surface_frames(output_dir: Path, variant_name: str) -> dict[str, Path]:
    return existing_visit_frames(output_dir, variant_name, ("surfaceDiagnostics", "particleSurfaceFlag"))


def run_visit_render(args: argparse.Namespace, source_dir: Path, output_dir: Path, results: list[dict], params: dict) -> None:
    result = representative_result(results)
    if result is None:
        return
    required = {"initial", "pct33", "pct66", "final"}
    existing_vm = existing_vm_frames(output_dir, result["name"])
    existing_surface = existing_surface_frames(output_dir, result["name"])
    if required.issubset(existing_vm.keys()) and required.issubset(existing_surface.keys()) and not args.rerender_visit:
        print(f"Using existing VisIt frames for {result['name']}; pass --rerender-visit to regenerate.")
        return
    if args.no_visit:
        return
    visit_cmd = find_visit_command(args.visit_cmd)
    if not visit_cmd:
        (output_dir / "visit_render_skipped.txt").write_text("VisIt command not found; set VISIT_COMMAND or pass --visit-cmd.\n")
        return
    run_dir = Path(result.get("run_dir") or "")
    if not run_dir.is_dir():
        return
    frame_dir = output_dir / "visit_frames" / result["name"]
    frame_dir.mkdir(parents=True, exist_ok=True)
    script = renderer_path(source_dir)
    final_ideal = ideal_reaction_stress(params["initial_gap"] + params["extra_compression"], params)
    stress_color_max = max(1.25 * final_ideal, 1.0e-6)

    render_specs = [
        {
            "name": "particleVonMisesStress",
            "variable": "particleVonMisesStress",
            "case_suffix": "particleVonMisesStress",
            "color_min": 0.0,
            "color_max": stress_color_max,
            "extra": [],
        },
        {
            "name": "surfaceDiagnostics",
            "variable": "particleSurfaceFlag",
            "case_suffix": "surfaceDiagnostics",
            "color_min": 0.0,
            "color_max": 2.0,
            "extra": [
                "--colortable", "Default",
                "--vector-variable", "particleSurfaceNormal",
                "--vector-variable", "particleSurfacePosition",
                "--vector-scale", "0.25",
                "--vector-count", "500",
            ],
        },
    ]

    for spec in render_specs:
        if spec["name"] == "particleVonMisesStress" and required.issubset(existing_vm.keys()) and not args.rerender_visit:
            continue
        if spec["name"] == "surfaceDiagnostics" and required.issubset(existing_surface.keys()) and not args.rerender_visit:
            continue
        cmd = [
            visit_cmd,
            "-nowin",
            "-cli",
            "-s",
            str(script),
            "--run-dir",
            str(run_dir),
            "--output-dir",
            str(frame_dir),
            "--case-name",
            f"{result['case_name']}_{spec['case_suffix']}",
            "--variable",
            spec["variable"],
            "--strict-variable",
            "--states",
            "initial,33%,66%,final",
            "--view",
            "xy",
            "--width",
            "1024",
            "--height",
            "1024",
            "--range-mode",
            "explicit",
            "--color-min",
            f"{spec['color_min']:.12g}",
            "--color-max",
            f"{spec['color_max']:.12g}",
        ] + spec.get("extra", [])
        proc = subprocess.run(cmd, cwd=source_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_path = output_dir / f"visit_{result['name']}_{spec['name']}.log"
        log_path.write_text(proc.stdout)
        if proc.returncode != 0:
            print(f"VisIt render failed for {result['name']} {spec['name']}; see {log_path}")
        else:
            print(f"VisIt render complete for {result['name']} {spec['name']}")


def write_summary(output_dir: Path, results: list[dict], params: dict) -> None:
    summary = {result["name"]: result.get("summary", {}) for result in results}
    for result in results:
        if result.get("error"):
            summary.setdefault(result["name"], {})["error"] = result["error"]
    summary["parameters"] = params
    (output_dir / "contact_surface_gap_closure_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))


def existing_rel(output_dir: Path, rel: str) -> bool:
    return (output_dir / rel).is_file()


def append_frame_panel(lines: list[str], output_dir: Path, frames: dict[str, Path], caption: str, labels: dict[str, str]) -> None:
    frame_order = [("initial", r"initial"), ("pct33", r"near closure"), ("pct66", r"post-closure"), ("final", r"final")]
    if not any(key in frames for key, _ in frame_order):
        return
    lines.extend([r"\begin{figure}[htbp]", r"\centering"])
    idx = 0
    for key, default_label in frame_order:
        frame = frames.get(key)
        if frame is None:
            continue
        try:
            rel = frame.relative_to(output_dir)
        except ValueError:
            continue
        if idx > 0 and idx % 2 == 0:
            lines.append(r"\par\vspace{0.5em}")
        label = labels.get(key, default_label)
        lines.extend([
            r"\begin{minipage}{0.48\textwidth}",
            r"\centering",
            r"\includegraphics[width=\linewidth]{\CaseOutputDir/" + str(rel).replace("\\", "/") + r"}",
            r"\par\small " + latex_escape(label),
            r"\end{minipage}" + (r"\hfill" if idx % 2 == 0 else ""),
        ])
        idx += 1
    lines.extend([r"\caption{" + caption + r"}", r"\end{figure}", ""])


def write_latex(output_dir: Path, results: list[dict], params: dict) -> None:
    lines = [
        r"\paragraph{Generated reaction metric.}",
        r"The primary metric compares the compressive y-reaction stress against an ideal response that is zero until the curved gap closes and then follows the plane-strain uniaxial-compression response of the shared hyperelastic material.",
        r"",
        r"\begin{center}",
        r"\scriptsize",
        r"\begin{tabular}{llrrrr}",
        r"\hline",
        r"Field/surface & Gap & Samples & Final stress & Final ideal & Post-closure rel. RMSE\\",
        r"\hline",
    ]
    for result in results:
        summary = result.get("summary", {})
        short = result["label"].replace(" / ", ", ")
        if not summary or summary.get("num_rows", 0) == 0:
            lines.append(r"{} & {} & \multicolumn{{4}}{{c}}{{{}}}\\".format(latex_escape(short), latex_escape(result.get("gap_correction", "")), latex_escape(result.get("error", "no reaction history"))))
            continue
        lines.append(
            r"{} & {} & {} & {} & {} & {}\\".format(
                latex_escape(result["field_mode"] + ", " + result["surface_mode"]),
                latex_escape(result["gap_correction"]),
                int(summary.get("num_rows", 0)),
                compact_float(summary.get("final_reaction_stress_y", 0.0)),
                compact_float(summary.get("final_ideal_reaction_stress_y", 0.0)),
                compact_float(summary.get("postclosure_relative_rmse_reaction_stress_y", 0.0)),
            )
        )
    lines += [
        r"\hline",
        r"\end{tabular}",
        r"\end{center}",
        r"\normalsize",
        r"",
    ]
    if existing_rel(output_dir, "contact_gap_reaction_stress.png"):
        lines += [
            r"\begin{figure}[htbp]",
            r"\centering",
            r"\includegraphics[width=0.80\textwidth]{\CaseOutputDir/contact_gap_reaction_stress.png}",
            r"\caption{Compressive y-reaction stress versus imposed displacement for all contact surface/gap-closure variants. The ideal curve is zero before gap closure and then follows the uniaxial-strain hyperelastic response.}",
            r"\end{figure}",
            r"",
        ]
    if existing_rel(output_dir, "contact_gap_reaction_error.png"):
        lines += [
            r"\begin{figure}[htbp]",
            r"\centering",
            r"\includegraphics[width=0.80\textwidth]{\CaseOutputDir/contact_gap_reaction_error.png}",
            r"\caption{Post-closure reaction-stress error relative to the ideal mated-block compression response.}",
            r"\end{figure}",
            r"",
        ]

    rep = representative_result(results)
    frame_labels = {
        "initial": r"initial separated state",
        "pct33": r"near curved-surface closure",
        "pct66": r"post-closure compression",
        "final": r"final compression",
    }
    if rep:
        vm_frames = existing_vm_frames(output_dir, rep["name"])
        surface_frames = existing_surface_frames(output_dir, rep["name"])
        if vm_frames:
            append_frame_panel(
                lines,
                output_dir,
                vm_frames,
                r"Particle von-Mises stress overlaid on the background grid for a representative contact surface/gap-closure subcase. The frames show the zero-stress separated state, contact near gap closure, post-closure load transfer, and final compression.",
                frame_labels,
            )
        else:
            lines += [
                r"\paragraph{VisIt von-Mises frames.}",
                r"No valid particle von-Mises stress frames were found in \texttt{\CaseOutputDir/visit\_frames}. Re-run the post-processor with VisIt enabled, for example \texttt{--visit-cmd visit --rerender-visit}.",
                r"",
            ]
        if surface_frames:
            append_frame_panel(
                lines,
                output_dir,
                surface_frames,
                r"Particle surface flags with background grid and vector overlays for particle surface normal and particle surface position. These frames verify that the explicit surface description follows the curved mating geometry before and during gap closure.",
                frame_labels,
            )
        else:
            lines += [
                r"\paragraph{VisIt surface-diagnostic frames.}",
                r"No valid particle surface-flag/normal/position frames were found in \texttt{\CaseOutputDir/visit\_frames}. Re-run the post-processor with VisIt enabled, for example \texttt{--visit-cmd visit --rerender-visit}.",
                r"",
            ]
    else:
        lines += [
            r"\paragraph{VisIt frames.}",
            r"No representative completed subcase was available for VisIt rendering.",
            r"",
        ]
    (output_dir / "contact_surface_gap_closure_results.tex").write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    results, params = collect_results(source_dir, output_dir, args.case_id)
    write_combined_csv(output_dir, results)
    write_summary(output_dir, results, params)
    make_plots(output_dir, results, params)
    run_visit_render(args, source_dir, output_dir, results, params)
    write_latex(output_dir, results, params)

    errors = [r for r in results if r.get("error")]
    if errors:
        print("Contact surface/gap-closure post-processing completed with missing data for: " + ", ".join(r["name"] for r in errors))
    else:
        print("Contact surface/gap-closure post-processing complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
