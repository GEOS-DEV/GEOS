#!/usr/bin/env python3
"""Aggregate post-processing for the PBC compaction momentum test.

The script reduces each subcase momentum-history CSV to a common table, plots the
total x-momentum versus time, optionally renders VisIt x-velocity contours, and
writes the LaTeX fragment included by test.tex.
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

VARIANTS = [
    {"name": "elastic", "label": "Elastic", "case_name": "pbcCompactionMomentum_elastic"},
    {"name": "plastic", "label": "Elastic-plastic", "case_name": "pbcCompactionMomentum_plastic"},
]
MOMENTUM_CSV = "mpm_momentumHistory.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Post-process all PBC compaction momentum variants")
    parser.add_argument("--suite", default="verification")
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--case-id", default="pbcCompactionMomentumTest")
    parser.add_argument("--python", dest="python_cmd", default=sys.executable)
    parser.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", ""))
    parser.add_argument("--no-visit", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def compact_float(value: float) -> str:
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
    manifest_path = output_dir / "pbcCompactionMomentumTest_jobs.json"
    if manifest_path.is_file():
        try:
            return json.loads(manifest_path.read_text())
        except Exception:
            pass
    subcases = []
    for variant in VARIANTS:
        cid = variant_case_id(case_id, variant)
        jobs_path = source_dir / "output" / cid / f"{variant['case_name']}_jobs.json"
        jobs = {}
        if jobs_path.is_file():
            try:
                jobs = json.loads(jobs_path.read_text())
            except Exception:
                jobs = {}
        subcases.append({"name": variant["name"], "label": variant["label"], "case_id": cid, "case_name": variant["case_name"], "jobs": jobs})
    return {"case_id": case_id, "subcases": subcases}


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


def find_momentum_csv(run_dir: Path) -> Path | None:
    direct = run_dir / MOMENTUM_CSV
    if direct.is_file():
        return direct
    matches = sorted(run_dir.glob(f"**/{MOMENTUM_CSV}"))
    return matches[0] if matches else None


def float_value(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default))
    except Exception:
        return default


def read_momentum_csv(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return rows
        for raw in reader:
            row = {key: float_value(raw, key, 0.0) for key in reader.fieldnames}
            row["cycle"] = int(float_value(raw, "cycle", len(rows)))
            rows.append(row)
    return rows


def parse_log_momentum(run_dir: Path) -> list[dict[str, float]]:
    """Fallback for older executables that only print momentum to text logs.

    The printed log does not include simulation time, so this fallback is useful for
    diagnostics but not for a publication-quality verification metric.
    """
    patterns = ["*.out", "*.log", "slurm*.out", "*.txt"]
    files = []
    for pattern in patterns:
        files.extend(sorted(run_dir.glob(pattern)))
    if not files:
        return []
    particle_re = re.compile(r"GLOBAL particle momentum sum = \[\s*([-+0-9.eE]+),\s*([-+0-9.eE]+),\s*([-+0-9.eE]+)\]")
    small_re = re.compile(r"GLOBAL nodal internal force sum \(small mass nodes only\) = \[\s*([-+0-9.eE]+),\s*([-+0-9.eE]+),\s*([-+0-9.eE]+)\]")
    rows: list[dict[str, float]] = []
    pending: dict[str, float] | None = None
    for path in files:
        try:
            lines = path.read_text(errors="replace").splitlines()
        except Exception:
            continue
        for line in lines:
            match = particle_re.search(line)
            if match:
                pending = {
                    "cycle": len(rows),
                    "time": float(len(rows)),
                    "dt": 0.0,
                    "particle_momentum_x": float(match.group(1)),
                    "particle_momentum_y": float(match.group(2)),
                    "particle_momentum_z": float(match.group(3)),
                }
                rows.append(pending)
                continue
            match = small_re.search(line)
            if match and rows:
                rows[-1]["nodal_internal_force_small_mass_x"] = float(match.group(1))
                rows[-1]["nodal_internal_force_small_mass_y"] = float(match.group(2))
                rows[-1]["nodal_internal_force_small_mass_z"] = float(match.group(3))
    return rows


def summarize(history: list[dict[str, float]]) -> dict[str, float]:
    if not history:
        return {"num_rows": 0}
    px = [row.get("particle_momentum_x", 0.0) for row in history]
    times = [row.get("time", float(i)) for i, row in enumerate(history)]
    small_fx = [row.get("nodal_internal_force_small_mass_x", 0.0) for row in history]
    active_fx = [row.get("nodal_internal_force_active_x", row.get("nodal_internal_force_x", 0.0)) for row in history]
    total_fx = [a + b for a, b in zip(active_fx, small_fx)]
    return {
        "num_rows": len(history),
        "initial_time": times[0],
        "final_time": times[-1],
        "initial_particle_momentum_x": px[0],
        "final_particle_momentum_x": px[-1],
        "max_abs_particle_momentum_x": max(abs(v) for v in px),
        "max_abs_small_mass_internal_force_x": max(abs(v) for v in small_fx) if small_fx else 0.0,
        "max_abs_total_internal_force_x": max(abs(v) for v in total_fx) if total_fx else 0.0,
    }


def collect_results(source_dir: Path, output_dir: Path, case_id: str) -> list[dict]:
    manifest = load_manifest(source_dir, output_dir, case_id)
    results = []
    for subcase in manifest.get("subcases", []):
        name = subcase.get("name", "")
        variant = next((v for v in VARIANTS if v["name"] == name), None)
        if variant is None:
            continue
        run_dir = run_dir_from_subcase(source_dir, subcase)
        result = {"name": name, "label": variant["label"], "case_name": variant["case_name"], "run_dir": str(run_dir) if run_dir else "", "history": [], "summary": {}}
        if run_dir is None:
            result["error"] = "run directory not found"
            results.append(result)
            continue
        csv_path = find_momentum_csv(run_dir)
        try:
            if csv_path:
                history = read_momentum_csv(csv_path)
            else:
                history = parse_log_momentum(run_dir)
            result["history"] = history
            result["summary"] = summarize(history)
            if csv_path:
                result["momentum_csv"] = str(csv_path)
            elif history:
                result["warning"] = "parsed momentum from text logs; time axis is sample index"
            else:
                result["error"] = f"{MOMENTUM_CSV} not found"
        except Exception as exc:
            result["error"] = str(exc)
        results.append(result)
    return results


def write_combined_csv(output_dir: Path, results: list[dict]) -> None:
    fieldnames = [
        "variant",
        "label",
        "cycle",
        "time",
        "dt",
        "particle_momentum_x",
        "particle_momentum_y",
        "particle_momentum_z",
        "inactive_particle_momentum_x",
        "inactive_particle_momentum_y",
        "inactive_particle_momentum_z",
        "nodal_momentum_x",
        "nodal_momentum_y",
        "nodal_momentum_z",
        "nodal_thresholded_momentum_x",
        "nodal_thresholded_momentum_y",
        "nodal_thresholded_momentum_z",
        "nodal_internal_force_active_x",
        "nodal_internal_force_active_y",
        "nodal_internal_force_active_z",
        "nodal_internal_force_small_mass_x",
        "nodal_internal_force_small_mass_y",
        "nodal_internal_force_small_mass_z",
        "nodal_mass",
    ]
    with (output_dir / "pbc_compaction_momentum_history.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for result in results:
            for row in result.get("history", []):
                out = {"variant": result["name"], "label": result["label"]}
                out.update(row)
                writer.writerow(out)


def make_plots(output_dir: Path, results: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    any_data = False
    for result in results:
        history = result.get("history", [])
        if not history:
            continue
        any_data = True
        times = [row.get("time", float(i)) for i, row in enumerate(history)]
        px = [row.get("particle_momentum_x", 0.0) for row in history]
        ax.plot(times, px, label=result["label"])
    if any_data:
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("time")
        ax.set_ylabel("total particle x-momentum")
        ax.set_title("PBC compaction x-momentum drift")
        ax.legend()
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(output_dir / "pbc_compaction_total_x_momentum.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    any_data = False
    for result in results:
        history = result.get("history", [])
        if not history:
            continue
        any_data = True
        times = [row.get("time", float(i)) for i, row in enumerate(history)]
        active = [row.get("nodal_internal_force_active_x", row.get("nodal_internal_force_x", 0.0)) for row in history]
        small = [row.get("nodal_internal_force_small_mass_x", 0.0) for row in history]
        total = [a + s for a, s in zip(active, small)]
        ax.plot(times, total, label=result["label"] + " total")
        ax.plot(times, small, linestyle="--", label=result["label"] + " small-mass")
    if any_data:
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("time")
        ax.set_ylabel("global nodal internal force in x")
        ax.set_title("Internal-force balance diagnostic")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(output_dir / "pbc_compaction_internal_force_balance.png", dpi=180)
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


def run_visit_render(args: argparse.Namespace, source_dir: Path, output_dir: Path, results: list[dict]) -> None:
    if args.no_visit:
        return
    visit_cmd = find_visit_command(args.visit_cmd)
    if not visit_cmd:
        (output_dir / "visit_render_skipped.txt").write_text("VisIt command not found; set VISIT_COMMAND or pass --visit-cmd.\n")
        return
    script = renderer_path(source_dir)
    for result in results:
        run_dir_text = result.get("run_dir") or ""
        if not run_dir_text:
            continue
        run_dir = Path(run_dir_text)
        frame_dir = output_dir / "visit_frames" / result["name"]
        frame_dir.mkdir(parents=True, exist_ok=True)
        for variable in ("particleVelocityX", "gridVelocityX"):
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
                f"{result['case_name']}_{variable}",
                "--variable",
                variable,
                "--states",
                "initial,quarter,middle,threequarter,final",
                "--view",
                "xy",
                "--range-mode",
                "auto",
            ]
            proc = subprocess.run(cmd, cwd=source_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            log_path = output_dir / f"visit_{result['name']}_{variable}.log"
            log_path.write_text(proc.stdout)
            if proc.returncode != 0:
                print(f"VisIt render failed for {result['name']} {variable}; see {log_path}")
            else:
                print(f"VisIt render complete for {result['name']} {variable}")


def write_summary(output_dir: Path, results: list[dict]) -> None:
    summary = {result["name"]: result.get("summary", {}) for result in results}
    for result in results:
        if result.get("error"):
            summary.setdefault(result["name"], {})["error"] = result["error"]
        if result.get("warning"):
            summary.setdefault(result["name"], {})["warning"] = result["warning"]
    (output_dir / "pbc_compaction_momentum_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))


def existing_rel(output_dir: Path, rel: str) -> bool:
    return (output_dir / rel).is_file()


def write_latex(output_dir: Path, results: list[dict]) -> None:
    lines = [
        r"\paragraph{Generated momentum metric.}",
        r"The primary verification metric is the maximum absolute total particle x-momentum. The exact value is zero because the initial x-momentum is zero and the imposed y-compaction does not prescribe net x-impulse. A nonzero drift therefore measures loss of discrete momentum conservation.",
        r"",
        r"\begin{center}",
        r"\begin{tabular}{lrrrr}",
        r"\hline",
        r"Subcase & Samples & Final time & Final $P_x$ & $\max |P_x|$\\",
        r"\hline",
    ]
    for result in results:
        summary = result.get("summary", {})
        if not summary or summary.get("num_rows", 0) == 0:
            lines.append(r"{} & \multicolumn{{4}}{{c}}{{{}}}\\".format(latex_escape(result["label"]), latex_escape(result.get("error", "no momentum history"))))
            continue
        lines.append(
            r"{} & {} & {} & {} & {}\\".format(
                latex_escape(result["label"]),
                int(summary.get("num_rows", 0)),
                compact_float(summary.get("final_time", 0.0)),
                compact_float(summary.get("final_particle_momentum_x", 0.0)),
                compact_float(summary.get("max_abs_particle_momentum_x", 0.0)),
            )
        )
    lines += [
        r"\hline",
        r"\end{tabular}",
        r"\end{center}",
        r"",
    ]
    if existing_rel(output_dir, "pbc_compaction_total_x_momentum.png"):
        lines += [
            r"\begin{figure}[htbp]",
            r"\centering",
            r"\includegraphics[width=0.72\textwidth]{\CaseOutputDir/pbc_compaction_total_x_momentum.png}",
            r"\caption{Total particle x-momentum versus time. The exact curve is identically zero.}",
            r"\end{figure}",
            r"",
        ]
    if existing_rel(output_dir, "pbc_compaction_internal_force_balance.png"):
        lines += [
            r"\begin{figure}[htbp]",
            r"\centering",
            r"\includegraphics[width=0.72\textwidth]{\CaseOutputDir/pbc_compaction_internal_force_balance.png}",
            r"\caption{Global x-component of the nodal internal-force balance, including the contribution carried by small-mass nodes.}",
            r"\end{figure}",
            r"",
        ]
    frame_root = output_dir / "visit_frames"
    if frame_root.is_dir():
        pngs = sorted(frame_root.glob("*/*.png"))[:12]
        if pngs:
            lines += [r"\paragraph{VisIt x-velocity frames.}"]
            for png in pngs:
                try:
                    rel = png.relative_to(output_dir)
                except ValueError:
                    continue
                lines += [
                    r"\begin{center}",
                    r"\includegraphics[width=0.45\textwidth]{\CaseOutputDir/" + str(rel).replace("\\", "/") + r"}",
                    r"\end{center}",
                ]
    (output_dir / "pbc_compaction_results.tex").write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    results = collect_results(source_dir, output_dir, args.case_id)
    write_combined_csv(output_dir, results)
    write_summary(output_dir, results)
    make_plots(output_dir, results)
    run_visit_render(args, source_dir, output_dir, results)
    write_latex(output_dir, results)

    errors = [r for r in results if r.get("error")]
    if errors:
        print("PBC compaction post-processing completed with missing data for: " + ", ".join(r["name"] for r in errors))
    else:
        print("PBC compaction post-processing complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
