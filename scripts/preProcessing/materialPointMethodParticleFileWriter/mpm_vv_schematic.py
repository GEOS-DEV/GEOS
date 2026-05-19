#!/usr/bin/env python3
"""Generate simple geometry-aware schematic PNGs from PFW inputs."""
from __future__ import annotations

import importlib.util
import os
import re
import sys
from pathlib import Path


def _matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
    return plt, Rectangle, Circle, FancyArrowPatch


def safe_case_label(case_id: str) -> str:
    return case_id.replace("__", " / ").replace("_", " ")


def load_input(input_path: Path, pfw_root: Path):
    old_path = list(sys.path)
    sys.path.insert(0, str(input_path.parent))
    sys.path.insert(0, str(pfw_root))
    try:
        spec = importlib.util.spec_from_file_location("_vv_pfw_input", input_path)
        mod = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(mod)
        return getattr(mod, "pfw", {})
    finally:
        sys.path[:] = old_path


def obj_bounds(obj):
    vals = []
    for name in ("xMin", "xMax"):
        fn = getattr(obj, name, None)
        if callable(fn):
            try:
                vals.append(fn())
            except Exception:
                vals.append(None)
    if len(vals) == 2 and vals[0] is not None and vals[1] is not None:
        try:
            return list(vals[0]), list(vals[1])
        except Exception:
            pass
    for a0, a1 in (("x0", "x1"), ("xmin", "xmax")):
        if hasattr(obj, a0) and hasattr(obj, a1):
            try:
                return list(getattr(obj, a0)), list(getattr(obj, a1))
            except Exception:
                pass
    return None


def draw_object(ax, obj, xmin, xmax, ymin, ymax, Rectangle, Circle, FancyArrowPatch, index):
    b = obj_bounds(obj)
    name = getattr(obj, "name", obj.__class__.__name__)
    cls = obj.__class__.__name__.lower()
    if b:
        lo, hi = b
        x0, y0 = lo[0], lo[1]
        w, h = hi[0] - lo[0], hi[1] - lo[1]
        if "cylinder" in cls or "sphere" in cls or abs(w - h) < 0.15 * max(abs(w), abs(h), 1e-12):
            circ = Circle((x0 + 0.5*w, y0 + 0.5*h), 0.5*max(abs(w), abs(h)), fill=False, lw=1.4)
            ax.add_patch(circ)
        else:
            ax.add_patch(Rectangle((x0, y0), w, h, fill=False, lw=1.4))
        if index < 5:
            ax.text(x0 + 0.5*w, y0 + 0.5*h, str(name), ha="center", va="center", fontsize=8)
    else:
        ax.text(0.5*(xmin+xmax), 0.5*(ymin+ymax), str(name), ha="center", va="center", fontsize=8)
    vel = getattr(obj, "vel", None)
    if vel is not None:
        try:
            vx, vy = float(vel[0]), float(vel[1])
            scale = 0.08 * max(xmax - xmin, ymax - ymin)
            if abs(vx) + abs(vy) > 0:
                cx = 0.5*(b[0][0] + b[1][0]) if b else 0.5*(xmin+xmax)
                cy = 0.5*(b[0][1] + b[1][1]) if b else 0.5*(ymin+ymax)
                norm = max((vx*vx+vy*vy)**0.5, 1e-12)
                ax.add_patch(FancyArrowPatch((cx, cy), (cx + scale*vx/norm, cy + scale*vy/norm), arrowstyle="->", mutation_scale=10, lw=1.0))
        except Exception:
            pass


def generate_schematic(input_path: Path, output_path: Path, case_id: str, pfw_root: Path) -> str:
    output_path = Path(output_path)
    try:
        pfw = load_input(Path(input_path), Path(pfw_root))
    except Exception as exc:
        pfw = {}
        error = str(exc)
    else:
        error = ""
    plt, Rectangle, Circle, FancyArrowPatch = _matplotlib()
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    xmin = float(pfw.get("xmin", 0.0)); xmax = float(pfw.get("xmax", 1.0))
    ymin = float(pfw.get("ymin", 0.0)); ymax = float(pfw.get("ymax", 1.0))
    ax.add_patch(Rectangle((xmin, ymin), xmax-xmin, ymax-ymin, fill=False, lw=1.6))
    nI = int(float(pfw.get("nI", 0) or 0)); nJ = int(float(pfw.get("nJ", 0) or 0))
    if 0 < nI <= 64:
        for i in range(1, nI):
            x = xmin + (xmax - xmin) * i / nI
            ax.plot([x, x], [ymin, ymax], lw=0.25, alpha=0.25)
    if 0 < nJ <= 64:
        for j in range(1, nJ):
            y = ymin + (ymax - ymin) * j / nJ
            ax.plot([xmin, xmax], [y, y], lw=0.25, alpha=0.25)
    objs = pfw.get("objects", []) or []
    for i, obj in enumerate(objs[:12]):
        draw_object(ax, obj, xmin, xmax, ymin, ymax, Rectangle, Circle, FancyArrowPatch, i)
    bc = pfw.get("boundaryConditionTypes", [])
    if bc:
        ax.text(xmin, ymax, "boundary types: " + str(bc), va="bottom", ha="left", fontsize=7)
    ax.set_aspect("equal", adjustable="box")
    pad = 0.08 * max(xmax-xmin, ymax-ymin, 1e-9)
    ax.set_xlim(xmin-pad, xmax+pad)
    ax.set_ylim(ymin-pad, ymax+pad)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    title = safe_case_label(case_id)
    if error:
        title += " (schematic fallback)"
        ax.text(0.5*(xmin+xmax), ymin-pad*0.6, error[:120], ha="center", va="top", fontsize=6)
    ax.set_title(title)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return str(output_path)


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--case-id", required=True)
    p.add_argument("--pfw-root", default=Path(__file__).resolve().parent)
    args = p.parse_args(argv)
    print(generate_schematic(Path(args.input), Path(args.output), args.case_id, Path(args.pfw_root)))


if __name__ == "__main__":
    main()
