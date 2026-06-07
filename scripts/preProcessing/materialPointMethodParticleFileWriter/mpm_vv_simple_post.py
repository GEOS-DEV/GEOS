#!/usr/bin/env python3
"""Common post-processing helpers for compact GEOS-MPM verification folders."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path


def parse_common_args(description: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--suite", default="verification")
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--python", dest="python_cmd", default=sys.executable)
    parser.add_argument("--visit-cmd", default=os.environ.get("VISIT_COMMAND", ""))
    parser.add_argument("--no-visit", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def latex_escape(value: object) -> str:
    text = str(value if value is not None else "")
    repl = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(repl.get(ch, ch) for ch in text)


def compact_float(value: object, digits: int = 6) -> str:
    try:
        v = float(value)
    except Exception:
        return "--"
    if not math.isfinite(v):
        return "--"
    return f"{v:.{digits}g}"


def variant_case_id(folder_case_id: str, variant: dict) -> str:
    return f"{folder_case_id}__{variant['name']}"


def load_manifest(source_dir: Path, output_dir: Path, case_id: str, default_name: str, variants: list[dict]) -> dict:
    for path in [output_dir / f"{default_name}_jobs.json", output_dir / f"{case_id}_jobs.json"]:
        if path.is_file():
            try:
                return json.loads(path.read_text())
            except Exception:
                pass
    subcases = []
    for variant in variants:
        cid = variant_case_id(case_id, variant)
        jobs_path = source_dir / "output" / cid / f"{variant['case_name']}_jobs.json"
        jobs = {}
        if jobs_path.is_file():
            try:
                jobs = json.loads(jobs_path.read_text())
            except Exception:
                jobs = {}
        subcases.append({"name": variant["name"], "label": variant.get("label", variant["name"]), "case_id": cid, "case_name": variant["case_name"], "jobs": jobs})
    return {"case_id": case_id, "subcases": subcases}


def run_dir_from_subcase(source_dir: Path, subcase: dict) -> Path | None:
    jobs = subcase.get("jobs", {}) if isinstance(subcase.get("jobs", {}), dict) else {}
    run_dir = jobs.get("run_dir")
    if run_dir:
        path = Path(run_dir)
        if path.is_dir():
            return path
    for candidate in [source_dir / "output" / str(subcase.get("case_id", "")) / str(subcase.get("case_name", "")), source_dir / "output" / str(subcase.get("case_id", ""))]:
        if candidate.is_dir():
            return candidate
    return None


def read_csv_numeric(path: Path) -> list[dict[str, float | str]]:
    if not path.is_file():
        return []
    rows: list[dict[str, float | str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        for raw in reader:
            row: dict[str, float | str] = {}
            for key, value in raw.items():
                if value is None or value == "":
                    row[key] = ""
                    continue
                try:
                    row[key] = float(value)
                except Exception:
                    row[key] = value
            rows.append(row)
    return rows


def find_first(run_dir: Path, names: list[str]) -> Path | None:
    for name in names:
        path = run_dir / name
        if path.is_file():
            return path
    for name in names:
        matches = sorted(run_dir.glob(f"**/{name}"))
        if matches:
            return matches[0]
    return None


def tracer_files(run_dir: Path, prefix: str) -> list[Path]:
    out: list[Path] = []
    for pattern in [f"{prefix}_*.csv", "tracerPoint_*.csv"]:
        out.extend(sorted(run_dir.glob(pattern)))
    return sorted(set(out))


def read_tracer(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        aliases = {"time": "t", "x": "x", "y": "y", "z": "z"}
        if "t" not in fields and "time" in fields:
            aliases["time"] = "time"
        missing = [aliases[k] for k in ["time", "x", "y", "z"] if aliases[k] not in fields]
        if missing:
            raise ValueError(f"{path.name} missing tracer columns {missing}")
        for raw in reader:
            row = {"t": float(raw[aliases["time"]]), "x": float(raw["x"]), "y": float(raw["y"]), "z": float(raw["z"])}
            for key, value in raw.items():
                if key in row or value in (None, ""):
                    continue
                try:
                    row[key] = float(value)
                except Exception:
                    pass
            rows.append(row)
    return rows


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True, default=str))


def matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_metric(output_dir: Path, file_name: str, title: str, x_label: str, y_label: str, series: list[tuple[str, list[float], list[float]]]) -> str | None:
    series = [(label, x, y) for label, x, y in series if x and y]
    if not series:
        return None
    plt = matplotlib()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for label, x, y in series:
        ax.plot(x, y, label=label)
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    out = output_dir / file_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def detect_visit(requested: str) -> str | None:
    for candidate in [requested, os.environ.get("VISIT_CMD", ""), os.environ.get("VISIT_COMMAND", ""), "/usr/gapps/visit/bin/visit", "visit"]:
        if not candidate:
            continue
        if os.path.isabs(candidate) and os.access(candidate, os.X_OK):
            return candidate
        found = shutil.which(candidate)
        if found:
            return found
    return None


def render_visit_frames(args: argparse.Namespace, run_dir: Path, output_dir: Path, case_name: str, variable: str, states: str = "initial,middle,final", view: str = "auto", colortable: str = "hot_desaturated", range_mode: str = "auto", color_min: float | None = None, color_max: float | None = None) -> list[str]:
    if args.no_visit:
        return []
    visit = detect_visit(args.visit_cmd)
    script = run_dir / "pfw_visit_render.py"
    if not visit or not script.is_file():
        return []
    frame_dir = output_dir / "visit_frames"
    frame_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        visit, "-nowin", "-cli", "-s", str(script),
        "--run-dir", str(run_dir),
        "--output-dir", str(frame_dir),
        "--case-name", case_name,
        "--states", states,
        "--view", view,
        "--variable", variable,
        "--colortable", colortable,
        "--range-mode", range_mode,
        "--list-databases",
    ]
    if color_min is not None and color_max is not None:
        cmd += ["--color-min", str(color_min), "--color-max", str(color_max), "--range-mode", "explicit"]
    proc = subprocess.run(cmd, cwd=run_dir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output_dir / f"{case_name}_visit_render.log").write_text(proc.stdout + f"\nreturncode={proc.returncode}\n")
    return [str(p) for p in sorted(frame_dir.glob(f"{case_name}_*.png"))]


def include_graphics(paths: list[Path], report_dir: Path, width: str = "0.31\\linewidth") -> str:
    lines = []
    for path in paths:
        try:
            rel = Path(os.path.relpath(path, report_dir)).as_posix()
        except Exception:
            rel = path.as_posix()
        lines.append(r"\includegraphics[width=" + width + "]{" + rel + "}")
    return "\n".join(lines)
