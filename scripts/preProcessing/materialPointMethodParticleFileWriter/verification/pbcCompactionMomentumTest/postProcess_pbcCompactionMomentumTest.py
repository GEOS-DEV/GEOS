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
    parser.add_argument("--no-visit", action="store_true", help="Do not launch VisIt; still include any existing valid velocity PNGs in the LaTeX fragment")
    parser.add_argument("--rerender-visit", action="store_true", help="Re-render VisIt frames even when matching PNGs already exist")
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


TEXT_MOMENTUM_COLUMNS = {"stage", "diagnostic_velocity_field"}
FORCE_BALANCE_STAGES = {"Force balance snapshot", "10b. After active grid-field mask update"}
END_STEP_STAGE = "End explicit step"


def read_momentum_csv(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return rows
        for raw in reader:
            row: dict[str, object] = {}
            for key in reader.fieldnames:
                if key in TEXT_MOMENTUM_COLUMNS:
                    row[key] = raw.get(key, "")
                else:
                    row[key] = float_value(raw, key, 0.0)
            row["cycle"] = int(float_value(raw, "cycle", len(rows)))
            rows.append(row)
    return rows


def metric_rows(history: list[dict[str, object]]) -> list[dict[str, object]]:
    """Rows used for production momentum metrics.

    New momentum ledgers carry a string stage column.  Use only the end-of-step
    rows for the reported particle momentum.  Legacy CSV files do not have a
    stage column and are already one-row-per-step, so keep them as-is.
    """
    if not history or "stage" not in history[0]:
        return history
    selected = [row for row in history if str(row.get("stage", "")) == END_STEP_STAGE]
    return selected if selected else history


def force_balance_rows(history: list[dict[str, object]]) -> list[dict[str, object]]:
    """Rows where the grid internal-force balance is meaningful.

    The end-of-step particle momentum is meaningful after grid cleanup and
    constitutive updates, but the grid force arrays are not a force-balance
    snapshot at that point.  When a stage ledger is available, use the single
    synchronized force-assembly snapshot.  Legacy one-row-per-step CSV files did
    not carry a stage column and are treated as force-balance rows.
    """
    if not history:
        return []
    if "stage" not in history[0]:
        return history
    return [row for row in history if str(row.get("stage", "")) in FORCE_BALANCE_STAGES]


def row_float(row: dict[str, object], key: str, default: float = 0.0) -> float:
    value = row.get(key, default)
    try:
        return float(value)
    except Exception:
        return default


def internal_force_total_x(row: dict[str, object]) -> float:
    if "nodal_internal_force_all_x" in row:
        return row_float(row, "nodal_internal_force_all_x", 0.0)
    return row_float(row, "nodal_internal_force_active_x", row_float(row, "nodal_internal_force_x", 0.0)) + row_float(row, "nodal_internal_force_small_mass_x", 0.0)


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


def summarize(history: list[dict[str, object]], force_history: list[dict[str, object]] | None = None) -> dict[str, float]:
    if not history:
        return {"num_rows": 0}
    px = [row_float(row, "particle_momentum_x", 0.0) for row in history]
    times = [row_float(row, "time", float(i)) for i, row in enumerate(history)]
    force_rows_for_summary = force_history or []
    small_fx = [row_float(row, "nodal_internal_force_small_mass_x", 0.0) for row in force_rows_for_summary]
    total_fx = [internal_force_total_x(row) for row in force_rows_for_summary]
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
                raw_history = read_momentum_csv(csv_path)
            else:
                raw_history = parse_log_momentum(run_dir)
            history = metric_rows(raw_history)
            force_history = force_balance_rows(raw_history)
            result["history"] = history
            result["force_history"] = force_history
            result["raw_history_rows"] = len(raw_history)
            result["summary"] = summarize(history, force_history)
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
        "stage",
        "diagnostic_velocity_field",
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
        "nodal_internal_force_all_x",
        "nodal_internal_force_all_y",
        "nodal_internal_force_all_z",
        "nodal_external_force_all_x",
        "nodal_contact_force_all_x",
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
        times = [row_float(row, "time", float(i)) for i, row in enumerate(history)]
        px = [row_float(row, "particle_momentum_x", 0.0) for row in history]
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
        history = result.get("force_history", [])
        if not history:
            continue
        any_data = True
        times = [row_float(row, "time", float(i)) for i, row in enumerate(history)]
        small = [row_float(row, "nodal_internal_force_small_mass_x", 0.0) for row in history]
        total = [internal_force_total_x(row) for row in history]
        ax.plot(times, total, label=result["label"] + " total")
        ax.plot(times, small, linestyle="--", label=result["label"] + " small-mass")
    force_plot = output_dir / "pbc_compaction_internal_force_balance.png"
    if any_data:
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("time")
        ax.set_ylabel("global nodal internal force in x")
        ax.set_title("Internal-force balance diagnostic")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(force_plot, dpi=180)
    elif force_plot.is_file():
        # Avoid carrying a stale force-balance plot from an older diagnostic run
        # into a report generated from a production end-of-step-only CSV.
        force_plot.unlink()
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
    name = path.stem.lower()
    for label in ("initial", "quarter", "middle", "threequarter", "final"):
        if f"_{label}_" in name or name.endswith(f"_{label}"):
            return label
    return ""


def is_valid_particle_velocity_frame(path: Path) -> bool:
    stem = path.stem.lower()
    # The old renderer fell back to particleDamage while the case-name still
    # contained particleVelocityX.  Do not include those misleading images in
    # the verification report.
    return "particlevelocityx" in stem and "particledamage" not in stem


def existing_particle_velocity_frames(output_dir: Path, source_dir: Path, variant_name: str) -> list[Path]:
    search_roots = [output_dir / "visit_frames" / variant_name, source_dir / "visit_frames" / variant_name, source_dir]
    frames: list[Path] = []
    seen: set[Path] = set()
    for root in search_roots:
        if not root.is_dir():
            continue
        for png in sorted(root.rglob("*.png")):
            try:
                resolved = png.resolve()
            except Exception:
                resolved = png
            if resolved in seen:
                continue
            if variant_name.lower() not in png.name.lower():
                continue
            if is_valid_particle_velocity_frame(png):
                frames.append(png)
                seen.add(resolved)
    return frames


def copy_existing_particle_velocity_frames(source_dir: Path, output_dir: Path, results: list[dict]) -> None:
    """Copy already-rendered valid particle-x-velocity frames into output/visit_frames.

    Older versions of pfw_visit_render.py left VisIt's outputToCurrentDirectory
    flag enabled, so SaveWindow wrote PNGs into the case source folder even when
    --output-dir was supplied.  Keep those frames usable without requiring a
    rerender, but deliberately ignore the known bad
    particleVelocityX_*_particleDamage.png fallback images.
    """
    for result in results:
        variant_name = result["name"]
        frame_dir = output_dir / "visit_frames" / variant_name
        frame_dir.mkdir(parents=True, exist_ok=True)
        for png in existing_particle_velocity_frames(output_dir, source_dir, variant_name):
            try:
                if png.parent.resolve() == frame_dir.resolve():
                    continue
            except Exception:
                pass
            target = frame_dir / png.name
            if not target.is_file():
                shutil.copy2(png, target)


def particle_velocity_frames_for_report(output_dir: Path, variant_name: str) -> dict[str, Path]:
    frame_dir = output_dir / "visit_frames" / variant_name
    out: dict[str, Path] = {}
    if not frame_dir.is_dir():
        return out
    for png in sorted(frame_dir.rglob("*.png")):
        if variant_name.lower() not in png.name.lower():
            continue
        if not is_valid_particle_velocity_frame(png):
            continue
        state = frame_state_label(png)
        if state and state not in out:
            out[state] = png
    return out


def run_visit_render(args: argparse.Namespace, source_dir: Path, output_dir: Path, results: list[dict]) -> None:
    copy_existing_particle_velocity_frames(source_dir, output_dir, results)
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
        variable = "particleVelocityX"
        existing = particle_velocity_frames_for_report(output_dir, result["name"])
        if {"middle", "final"}.issubset(existing.keys()) and not args.rerender_visit:
            print(
                f"Using existing particle x-velocity frames for {result['name']}; "
                "pass --rerender-visit to regenerate."
            )
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
            f"{result['case_name']}_{variable}",
            "--variable",
            variable,
            "--strict-variable",
            "--states",
            "initial,middle,final",
            "--view",
            "xy",
            "--width",
            "1024",
            "--height",
            "1024",
            "--range-mode",
            "explicit",
            "--color-min",
            "-0.05",
            "--color-max",
            "0.05",
        ]
        proc = subprocess.run(cmd, cwd=source_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_path = output_dir / f"visit_{result['name']}_{variable}.log"
        log_path.write_text(proc.stdout)
        # Copy frames again to catch older VisIt builds that ignored
        # outputDirectory and saved into cwd despite the renderer request.
        copy_existing_particle_velocity_frames(source_dir, output_dir, [result])
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
            r"\caption{Global x-component of the nodal internal-force balance at the synchronized force-assembly snapshot, including the contribution carried by small-mass nodes.}",
            r"\end{figure}",
            r"",
        ]
    else:
        lines += [
            r"\paragraph{Internal-force balance snapshot.}",
            r"No synchronized force-assembly snapshot was available in the momentum-history CSV. The production end-of-step momentum rows remain valid for the conservation metric, but the grid force arrays are not plotted after grid cleanup because they no longer represent a force-balance snapshot.",
            r"",
        ]
    # Include a compact 2x2 panel: elastic/plastic by half-time/final.
    # These are particle x-velocity pseudocolor plots overlaid on the background
    # grid; initial frames are usually nearly uniform and are omitted from the
    # report to keep the suite section concise.
    selected_frames = []
    for result in results:
        by_state = particle_velocity_frames_for_report(output_dir, result["name"])
        for state, label in (("middle", r"$t \approx T/2$"), ("final", r"$t=T$")):
            frame = by_state.get(state)
            if frame is not None:
                selected_frames.append((result["label"], label, frame))
    if selected_frames:
        lines += [
            r"\begin{figure}[htbp]",
            r"\centering",
        ]
        for idx, (variant_label, state_label, frame) in enumerate(selected_frames):
            try:
                rel = frame.relative_to(output_dir)
            except ValueError:
                continue
            if idx > 0 and idx % 2 == 0:
                lines.append(r"\par\vspace{0.5em}")
            lines += [
                r"\begin{minipage}{0.48\textwidth}",
                r"\centering",
                r"\includegraphics[width=\linewidth]{\CaseOutputDir/" + str(rel).replace("\\", "/") + r"}",
                r"\par\small " + latex_escape(variant_label) + r", " + state_label,
                r"\end{minipage}" + (r"\hfill" if idx % 2 == 0 else ""),
            ]
        lines += [
            r"\caption{Particle x-velocity contours overlaid on the background grid for the PBC compaction momentum test.  The physically expected total x-momentum is zero, so these frames visualize local transverse velocity while the momentum plot quantifies the global conservation error.}",
            r"\end{figure}",
            r"",
        ]
    else:
        lines += [
            r"\paragraph{VisIt x-velocity frames.}",
            r"No valid particle x-velocity frames were found in \texttt{\CaseOutputDir/visit\_frames}.  Re-run the post-processor with VisIt enabled, for example \texttt{--visit-cmd visit --rerender-visit}.",
            r"",
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
