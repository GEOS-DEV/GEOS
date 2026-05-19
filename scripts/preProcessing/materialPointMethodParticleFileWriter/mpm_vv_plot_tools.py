#!/usr/bin/env python3
"""Generic verification history plotting utilities for GEOS-MPM.

The routines here intentionally make useful first-pass verification plots from
reactionHistory.csv and boxAverageHistory.csv even when a legacy case-specific
plot script fails.  Case-specific scripts are still run opportunistically by the
post-processor.
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
from pathlib import Path


def _matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _as_float(value):
    try:
        if value is None:
            return None
        text = str(value).strip()
        if not text or text.lower() in {"nan", "none"}:
            return None
        return float(text)
    except Exception:
        return None


def read_table(path: Path):
    path = Path(path)
    if not path.is_file():
        return [], []
    text = path.read_text(errors="replace").splitlines()
    text = [line for line in text if line.strip() and not line.lstrip().startswith("#")]
    if not text:
        return [], []
    first = text[0]
    delimiter = "," if "," in first else None
    if delimiter:
        rows = list(csv.reader(text, delimiter=delimiter))
    else:
        rows = [re.split(r"\s+", line.strip()) for line in text]
    if not rows:
        return [], []
    header = [h.strip() for h in rows[0]]
    numeric_rows = []
    for row in rows[1:]:
        vals = [_as_float(x) for x in row]
        if any(v is not None for v in vals):
            numeric_rows.append(vals)
    return header, numeric_rows


def time_column(header):
    for i, h in enumerate(header):
        low = h.lower()
        if low in {"time", "t"} or "time" in low:
            return i
    return 0


def column_values(rows, i):
    out = []
    for row in rows:
        out.append(row[i] if i < len(row) else None)
    return out


def nonempty(vals):
    return any(v is not None and math.isfinite(v) for v in vals)


def plot_columns(path: Path, out: Path, title: str, ylabel: str, selector, max_lines: int = 12):
    header, rows = read_table(path)
    if len(rows) < 2 or not header:
        return None
    t_idx = time_column(header)
    x = column_values(rows, t_idx)
    candidates = [i for i, h in enumerate(header) if i != t_idx and selector(h, i)]
    if not candidates:
        return None
    plt = _matplotlib()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    plotted = 0
    for i in candidates[:max_lines]:
        y = column_values(rows, i)
        if nonempty(y):
            ax.plot(x, y, label=header[i] if i < len(header) else f"col{i}")
            plotted += 1
    if plotted == 0:
        plt.close(fig)
        return None
    ax.set_xlabel(header[t_idx] if t_idx < len(header) else "time")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def _contains_any(text, needles):
    low = text.lower()
    return any(n in low for n in needles)


def reaction_plot(run_dir: Path, prefix: str):
    path = run_dir / "reactionHistory.csv"
    return plot_columns(
        path,
        run_dir / f"{prefix}_reaction.png",
        f"{prefix}: reaction history",
        "reaction / force",
        lambda h, i: _contains_any(h, ["reaction", "force", "rx", "ry", "rz", "_x", "_y", "_z"]),
        max_lines=12,
    )


def box_group_plot(run_dir: Path, prefix: str, group_name: str, keywords, ylabel: str):
    path = run_dir / "boxAverageHistory.csv"
    return plot_columns(
        path,
        run_dir / f"{prefix}_box_{group_name}.png",
        f"{prefix}: box-average {group_name}",
        ylabel,
        lambda h, i: _contains_any(h, keywords),
        max_lines=10,
    )


def derived_pressure_plot(run_dir: Path, prefix: str):
    path = run_dir / "boxAverageHistory.csv"
    header, rows = read_table(path)
    if len(rows) < 2:
        return None
    names = [h.lower() for h in header]
    t_idx = time_column(header)
    def find_component(options):
        for opt in options:
            for i, name in enumerate(names):
                if opt in name and "stress" in name:
                    return i
        return None
    ixx = find_component(["stress_11", "stress_0", "sxx", "xx"])
    iyy = find_component(["stress_22", "stress_1", "syy", "yy"])
    izz = find_component(["stress_33", "stress_2", "szz", "zz"])
    if ixx is None or iyy is None or izz is None:
        return None
    x = column_values(rows, t_idx)
    p = []
    for row in rows:
        vals = []
        for i in (ixx, iyy, izz):
            vals.append(row[i] if i < len(row) else None)
        if all(v is not None for v in vals):
            p.append(-(vals[0] + vals[1] + vals[2]) / 3.0)
        else:
            p.append(None)
    if not nonempty(p):
        return None
    plt = _matplotlib()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(x, p, label="pressure = -tr(stress)/3")
    ax.set_xlabel(header[t_idx] if t_idx < len(header) else "time")
    ax.set_ylabel("pressure")
    ax.set_title(f"{prefix}: derived pressure")
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = run_dir / f"{prefix}_box_pressure.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def component_magnitude_plot(run_dir: Path, prefix: str, family: str, keywords, out_name: str, ylabel: str):
    path = run_dir / "boxAverageHistory.csv"
    header, rows = read_table(path)
    if len(rows) < 2:
        return None
    low = [h.lower() for h in header]
    t_idx = time_column(header)
    comps = []
    for axis in ("0", "1", "2", "x", "y", "z"):
        for i, h in enumerate(low):
            if i == t_idx:
                continue
            if any(k in h for k in keywords) and (h.endswith(axis) or ("_" + axis) in h or axis in h[-3:]):
                if i not in comps:
                    comps.append(i)
                    break
        if len(comps) >= 3:
            break
    if len(comps) < 2:
        return None
    x = column_values(rows, t_idx)
    mag = []
    for row in rows:
        vals = []
        for i in comps[:3]:
            vals.append(row[i] if i < len(row) else None)
        vals = [v for v in vals if v is not None]
        mag.append(math.sqrt(sum(v*v for v in vals)) if vals else None)
    if not nonempty(mag):
        return None
    plt = _matplotlib()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(x, mag, label=family + " magnitude")
    ax.set_xlabel(header[t_idx] if t_idx < len(header) else "time")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{prefix}: {family} magnitude")
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = run_dir / f"{prefix}_{out_name}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out)


def generate_standard_plots(run_dir: Path, prefix: str):
    run_dir = Path(run_dir)
    products = []
    for p in [
        reaction_plot(run_dir, prefix),
        box_group_plot(run_dir, prefix, "stress", ["stress", "sxx", "syy", "szz"], "stress"),
        box_group_plot(run_dir, prefix, "strain", ["strain", "deformation", "stretch"], "strain"),
        derived_pressure_plot(run_dir, prefix),
        box_group_plot(run_dir, prefix, "damage", ["damage", "porosity", "strengthscale"], "damage / porosity / strength"),
        box_group_plot(run_dir, prefix, "energy", ["energy", "work"], "energy / work"),
        box_group_plot(run_dir, prefix, "temperature", ["temperature", "temp"], "temperature"),
        component_magnitude_plot(run_dir, prefix, "velocity", ["velocity", "vel"], "box_velocity_magnitude", "velocity magnitude"),
        component_magnitude_plot(run_dir, prefix, "momentum", ["momentum"], "box_momentum_magnitude", "momentum magnitude"),
        component_magnitude_plot(run_dir, prefix, "plastic strain", ["plastic"], "box_plastic_strain_magnitude", "plastic strain magnitude"),
    ]:
        if p:
            products.append(p)
    return products


def analyze_plot_script(path: Path):
    text = path.read_text(errors="replace") if Path(path).is_file() else ""
    low = text.lower()
    sources = []
    for name in ["reactionHistory.csv", "boxAverageHistory.csv", "particleVelocityHistory.csv", "momentumHistory.csv"]:
        if name.lower() in low:
            sources.append(name)
    if not sources:
        if "reaction" in low:
            sources.append("reactionHistory.csv")
        if "box" in low or "average" in low:
            sources.append("boxAverageHistory.csv")
    variables = []
    for key in ["reaction", "force", "stress", "strain", "damage", "plastic", "pressure", "temperature", "velocity", "momentum", "energy", "porosity", "strength"]:
        if key in low:
            variables.append(key)
    processing = []
    if "area" in low or "stress" in low and ("force" in low or "reaction" in low):
        processing.append("derive stress/traction from reaction force divided by area/length")
    if "solidvolumefraction" in low or "volume fraction" in low or "porosity" in low:
        processing.append("scale by solid volume fraction or porosity")
    if "j2" in low or "vonmises" in low or "invariant" in low or "sqrt" in low:
        processing.append("compute stress/strain invariant or vector/tensor magnitude")
    if "pressure" in low or "hyd" in path.name.lower():
        processing.append("derive pressure or hydrostatic stress from stress components")
    if "reference" in low or "analytic" in low or "analytical" in low:
        processing.append("compare with analytical/reference curve")
    return {
        "script": Path(path).name,
        "data_sources": sources or ["reactionHistory.csv when present", "boxAverageHistory.csv when present"],
        "variables": variables or ["history columns inferred from available CSV files"],
        "processing": processing or ["direct plotting with generic unit conversion where possible"],
    }


def make_plot_plan(source_dir: Path, case_id: str, family: str):
    entries = []
    for p in sorted(Path(source_dir).glob("plot*.py")):
        entries.append(analyze_plot_script(p))
    # Also look one directory up for families such as geomechanics where scripts are shared.
    if not entries:
        for p in sorted(Path(source_dir).parent.glob("plot*.py")):
            entries.append(analyze_plot_script(p))
    if not entries:
        entries.append({
            "script": "generic-history-plots",
            "data_sources": ["reactionHistory.csv when present", "boxAverageHistory.csv when present"],
            "variables": ["reaction/force components", "stress/strain/pressure/plastic/damage/velocity/energy/momentum columns"],
            "processing": ["pressure = -tr(stress)/3 when stress components are available", "component magnitudes for velocity/plastic strain/momentum when component triplets are available"],
        })
    return {"case_id": case_id, "family": family, "plots": entries}


def write_plan(path: Path, plan) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(plan, indent=2))
