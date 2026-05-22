#!/usr/bin/env python3
"""Generate publication-style schematics for GEOS-MPM example and verification reports.

This script is intentionally report-side only: it does not run GEOS and it does not
modify PFW inputs. It inspects PFW input files, geometry-object extents where they
can be evaluated safely, and case-name/code-feature patterns to create high-quality
schematic PNG/SVG assets for the examples and verification suite reports.

Typical usage from the GEOS repository root:
  python scripts/preProcessing/materialPointMethodParticleFileWriter/mpm_publication_schematics.py \
    --suite examples --suite verification --compile

The script overwrites the generated report schematic files in:
  scripts/preProcessing/materialPointMethodParticleFileWriter/examples/examples_suite_report/schematics
  scripts/preProcessing/materialPointMethodParticleFileWriter/verification/verification_suite_report/schematics

Those schematic directories are generated report products and should not be tracked
in git. Commit this source script instead.
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import importlib
import io
import json
import math
import os
import re
import runpy
import shutil
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.collections import PatchCollection
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch


@dataclasses.dataclass
class CaseInfo:
    suite: str
    case_id: str
    family: str
    input_path: Path | None
    suite_root: Path
    pfw_root: Path
    pfw: dict[str, Any]
    namespace: dict[str, Any]
    source_text: str
    report_image: Path | None = None
    diagnostics: list[str] = dataclasses.field(default_factory=list)


def sanitize_id(text: str) -> str:
    text = text.replace(" ", "_").replace("/", "__")
    text = re.sub(r"[^A-Za-z0-9_.+=@-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "case"


def human(text: str) -> str:
    text = text.replace("__", " / ").replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def short_num(x: float, unit: str = "") -> str:
    if not np.isfinite(x):
        return ""
    ax = abs(x)
    if ax == 0:
        s = "0"
    elif ax < 1e-3 or ax >= 1e4:
        s = f"{x:.2e}"
    elif ax < 0.1:
        s = f"{x:.3g}"
    else:
        s = f"{x:.4g}"
    return s + unit


def repo_root_from(start: Path) -> Path:
    p = start.resolve()
    for q in [p, *p.parents]:
        if (q / "scripts/preProcessing/materialPointMethodParticleFileWriter").is_dir():
            return q
    return p


def pfw_root_from(repo):
    """Return the PFW root, recovering from wrapper path mistakes.

    Preferred input is the repository root.  For convenience, this function also
    searches ancestors of the supplied path, the current working directory, and a
    common child named 'geos'.
    """
    from pathlib import Path
    candidates = []
    raw = Path(repo).expanduser().resolve()
    candidates.append(raw)
    candidates.extend(raw.parents)
    cwd = Path.cwd().resolve()
    candidates.append(cwd)
    candidates.extend(cwd.parents)
    # Common mistake from old wrappers: repo was supplied as the parent of ~/geos.
    candidates.extend([raw / "geos", raw / "GEOS", raw / "geosx"])

    seen = set()
    for base in candidates:
        if base in seen:
            continue
        seen.add(base)
        root = base / "scripts" / "preProcessing" / "materialPointMethodParticleFileWriter"
        if root.exists():
            return root
    root = raw / "scripts" / "preProcessing" / "materialPointMethodParticleFileWriter"
    raise FileNotFoundError(f"Could not find PFW root. Tried repo={raw}; expected {root}")

def safe_eval_input(path: Path, pfw_root: Path) -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    ns: dict[str, Any] = {
        "__file__": str(path),
        "__name__": "__mpm_schematic_input__",
    }
    diagnostics: list[str] = []
    text = path.read_text(errors="replace")
    old_cwd = Path.cwd()
    old_path = list(sys.path)
    sys.path.insert(0, str(pfw_root))
    sys.path.insert(0, str(path.parent))
    try:
        os.chdir(str(path.parent))
        try:
            code = compile(text, str(path), "exec")
            exec(code, ns)
        except SystemExit as exc:
            diagnostics.append(f"input raised SystemExit: {exc}")
        except Exception as exc:
            diagnostics.append(f"input execution warning: {type(exc).__name__}: {exc}")
            diagnostics.append(traceback.format_exc(limit=4))
    finally:
        os.chdir(str(old_cwd))
        sys.path[:] = old_path
    pfw = ns.get("pfw")
    if not isinstance(pfw, dict):
        pfw = {}
    return pfw, ns, diagnostics


def find_pfw_inputs(root: Path) -> list[Path]:
    paths: list[Path] = []
    for p in sorted(root.rglob("pfw_input*.py")):
        parts = set(p.parts)
        if "output" in parts or "_suite_runs" in parts or "suite_report" in p.name:
            continue
        stem = p.stem
        if any(tok in stem for tok in ("XXXX", "YYYY", "ZZZZ")):
            continue
        paths.append(p)
    return paths


def strip_input_prefix(stem: str) -> str:
    for prefix in ("pfw_input_", "pfw_input"):
        if stem.startswith(prefix):
            stem = stem[len(prefix):]
    return sanitize_id(stem)


def discover_cases(suite: str, pfw_root: Path) -> dict[str, CaseInfo]:
    suite_root = pfw_root / suite
    cases: dict[str, CaseInfo] = {}
    for inp in find_pfw_inputs(suite_root):
        source_text = inp.read_text(errors="replace")
        pfw, ns, diags = safe_eval_input(inp, pfw_root)
        rel_dir = inp.parent.relative_to(suite_root)
        rel_parts = [sanitize_id(x) for x in rel_dir.parts if x and x != "."]
        leaf = strip_input_prefix(inp.stem)
        if suite == "examples":
            parent = sanitize_id(inp.parent.name)
            if parent == "borehole":
                case_id = "borehole_Ghareb" if "Ghareb" in leaf or "ghareb" in source_text.lower() else leaf
            elif parent == "pdc":
                case_id = "pdc" if "r4" not in leaf.lower() else "pdc_r4"
            elif len(list(inp.parent.glob("pfw_input*.py"))) == 1:
                case_id = parent
            else:
                case_id = leaf
        else:
            case_id = "__".join(rel_parts + [leaf]) if rel_parts else leaf
        family = classify_family(case_id, inp, source_text)
        ci = CaseInfo(suite, case_id, family, inp, suite_root, pfw_root, pfw, ns, source_text, diagnostics=diags)
        for key in case_keys(ci):
            cases.setdefault(key, ci)
    return cases


def case_keys(ci: CaseInfo) -> list[str]:
    keys = {ci.case_id, sanitize_id(ci.case_id), ci.case_id.replace("__", "_")}
    if ci.input_path is not None:
        leaf = strip_input_prefix(ci.input_path.stem)
        keys.add(leaf)
        keys.add(sanitize_id(ci.input_path.parent.name))
        try:
            rel_dir = ci.input_path.parent.relative_to(ci.suite_root)
            keys.add("__".join([sanitize_id(x) for x in rel_dir.parts] + [leaf]))
        except Exception:
            pass
    return [k for k in keys if k]


def classify_family(case_id: str, path: Path | None, text: str) -> str:
    s = f"{case_id} {path or ''} {text[:5000]}".lower()
    if "braziliandisk" in s or "brazilian" in s:
        return "brazilian_disk"
    if "benchy" in s or "stl" in s:
        return "benchy"
    if "borehole" in s:
        return "borehole"
    if "pdc" in s or "cutter" in s:
        return "pdc"
    if "colliding" in s:
        return "colliding_disks"
    if "elasticdisk" in s:
        return "elastic_disk"
    if "spinningdisk" in s or "spinning_disk" in s:
        return "spinning_disk"
    if "expandingring" in s or "expanding_ring" in s:
        return "expanding_ring"
    if "expandingbar" in s or "expanding_bar" in s:
        return "expanding_bar"
    if "cohesive" in s or "cz" in s or "shearpeel" in s:
        return "cohesive_zone"
    if "geomechanics" in s or "geomodel" in s or "buckling" in s:
        return "geomechanics"
    if "ceramic" in s:
        return "ceramic"
    if "contactbc" in s or "ballsmack" in s or "bouncy" in s:
        return "contact_bc"
    if "contact" in s or "interface" in s:
        return "contact"
    if "ftable" in s or "elasticblock" in s:
        return "ftable"
    if "anneal" in s:
        return "anneal"
    if "cubic" in s:
        return "cubic"
    if "diskthrupartitions" in s or "partition" in s:
        return "partition"
    if "periodic" in s or "pbc" in s:
        return "periodic"
    if "momentum" in s:
        return "momentum"
    if "gas" in s:
        return "gas"
    if "vortex" in s:
        return "vortex"
    if "implicitfluid" in s or "fluidpressure" in s or "fluid" in s:
        return "implicit_fluid"
    if "materialswap" in s:
        return "material_swap"
    if "polymerheal" in s:
        return "polymer_heal"
    if "polymer" in s:
        return "polymer_cz"
    if "shrink" in s or "calcination" in s:
        return "shrinkage"
    if "sizeeffect" in s or "notched" in s:
        return "notched_bar"
    if "stresscontrol" in s:
        return "stress_control"
    if "temperature" in s or "temptable" in s:
        return "temperature_table"
    if "vonmises" in s:
        return "von_mises"
    if "triply" in s:
        return "triply_periodic"
    return "generic"


def get_num(pfw: dict[str, Any], *keys: str, default: float | None = None) -> float | None:
    for key in keys:
        if key in pfw:
            try:
                return float(pfw[key])
            except Exception:
                pass
    return default


def domain_bounds(ci: CaseInfo) -> tuple[float, float, float, float, float, float]:
    pfw = ci.pfw
    xmin = get_num(pfw, "xmin", "xMin", "x_min", "Xmin", default=None)
    xmax = get_num(pfw, "xmax", "xMax", "x_max", "Xmax", default=None)
    ymin = get_num(pfw, "ymin", "yMin", "y_min", "Ymin", default=None)
    ymax = get_num(pfw, "ymax", "yMax", "y_max", "Ymax", default=None)
    zmin = get_num(pfw, "zmin", "zMin", "z_min", "Zmin", default=0.0)
    zmax = get_num(pfw, "zmax", "zMax", "z_max", "Zmax", default=1.0)
    obj_ext = []
    for obj in geometry_objects(ci):
        ext = object_extents(obj)
        if ext is not None:
            obj_ext.append(ext)
    if xmin is None or xmax is None:
        if obj_ext:
            xmin = min(e[0] for e in obj_ext) if xmin is None else xmin
            xmax = max(e[1] for e in obj_ext) if xmax is None else xmax
        else:
            xmin, xmax = -0.75, 0.75
    if ymin is None or ymax is None:
        if obj_ext:
            ymin = min(e[2] for e in obj_ext) if ymin is None else ymin
            ymax = max(e[3] for e in obj_ext) if ymax is None else ymax
        else:
            ymin, ymax = -0.75, 0.75
    if zmin is None:
        zmin = 0.0
    if zmax is None:
        zmax = 1.0
    # Expand degenerate bounds.
    if abs(xmax - xmin) < 1e-12:
        xmin -= 0.5; xmax += 0.5
    if abs(ymax - ymin) < 1e-12:
        ymin -= 0.5; ymax += 0.5
    return float(xmin), float(xmax), float(ymin), float(ymax), float(zmin), float(zmax)


def is_geometry(obj: Any) -> bool:
    return hasattr(obj, "isInterior") and (hasattr(obj, "xMin") or hasattr(obj, "getSubregions") or hasattr(obj, "name"))


def geometry_objects(ci: CaseInfo) -> list[Any]:
    seen: set[int] = set()
    out: list[Any] = []
    def rec(x: Any, depth: int = 0) -> None:
        if depth > 5:
            return
        oid = id(x)
        if oid in seen:
            return
        seen.add(oid)
        if is_geometry(x):
            out.append(x)
            return
        if isinstance(x, dict):
            for v in x.values():
                rec(v, depth + 1)
        elif isinstance(x, (list, tuple, set)):
            for v in x:
                rec(v, depth + 1)
    rec(ci.pfw)
    for k, v in ci.namespace.items():
        if not k.startswith("__"):
            rec(v)
    return out


def call0(obj: Any, name: str) -> float | None:
    if not hasattr(obj, name):
        return None
    try:
        val = getattr(obj, name)
        val = val() if callable(val) else val
        return float(val)
    except Exception:
        return None


def object_extents(obj: Any) -> tuple[float, float, float, float] | None:
    xmin = call0(obj, "xMin")
    xmax = call0(obj, "xMax")
    ymin = call0(obj, "yMin")
    ymax = call0(obj, "yMax")
    if None not in (xmin, xmax, ymin, ymax):
        return xmin, xmax, ymin, ymax
    # fallback using common attributes
    center = None
    for attr in ("center", "c", "origin"):
        if hasattr(obj, attr):
            with contextlib.suppress(Exception):
                arr = np.array(getattr(obj, attr), dtype=float).ravel()
                if arr.size >= 2:
                    center = arr
                    break
    radius = None
    for attr in ("radius", "r", "R"):
        if hasattr(obj, attr):
            with contextlib.suppress(Exception):
                radius = float(getattr(obj, attr))
                break
    if center is not None and radius is not None:
        return center[0]-radius, center[0]+radius, center[1]-radius, center[1]+radius
    return None


def object_velocity(obj: Any) -> tuple[float, float] | None:
    for attr in ("vel", "velocity", "v", "initialVelocity"):
        if hasattr(obj, attr):
            try:
                val = getattr(obj, attr)
                val = val() if callable(val) else val
                arr = np.array(val, dtype=float).ravel()
                if arr.size >= 2 and np.linalg.norm(arr[:2]) > 0:
                    return float(arr[0]), float(arr[1])
            except Exception:
                pass
    return None


def add_title(ax, ci: CaseInfo, subtitle: str | None = None) -> None:
    title = human(ci.case_id)
    ax.text(0.0, 1.035, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=13, fontweight="bold")
    if subtitle:
        ax.text(0.0, 1.005, subtitle, transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5, color="0.25")


def add_scale_bar(ax, xmin, xmax, ymin, ymax, label: str | None = None) -> None:
    w = xmax - xmin; h = ymax - ymin
    L = nice_length(w * 0.25)
    x0 = xmin + 0.06*w; y0 = ymin + 0.075*h
    ax.plot([x0, x0+L], [y0, y0], lw=2.2, color="0.05", solid_capstyle="butt")
    ax.text(x0+L/2, y0-0.035*h, label or short_num(L), ha="center", va="top", fontsize=8)


def nice_length(x: float) -> float:
    if x <= 0 or not np.isfinite(x):
        return 1.0
    e = math.floor(math.log10(x))
    m = x / (10**e)
    for val in (1, 2, 5, 10):
        if m <= val:
            return val * 10**e
    return 10**(e+1)


def setup_axes(ci: CaseInfo, figsize=(7.2, 5.1)):
    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    fig.patch.set_facecolor("white")
    return fig, ax


def draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=None, fill="#ffffff"):
    w, h = xmax - xmin, ymax - ymin
    ax.add_patch(patches.Rectangle((xmin, ymin), w, h, facecolor=fill, edgecolor="0.08", lw=1.8, zorder=1))
    if grid:
        nx = 8 if w >= h else max(4, int(round(8*w/h)))
        ny = 8 if h >= w else max(4, int(round(8*h/w)))
        for x in np.linspace(xmin, xmax, nx + 1)[1:-1]:
            ax.plot([x, x], [ymin, ymax], color="0.85", lw=0.45, zorder=2)
        for y in np.linspace(ymin, ymax, ny + 1)[1:-1]:
            ax.plot([xmin, xmax], [y, y], color="0.85", lw=0.45, zorder=2)
    if partitions:
        xpar, ypar = partitions
        if xpar and xpar > 1:
            for x in np.linspace(xmin, xmax, int(xpar)+1)[1:-1]:
                ax.plot([x, x], [ymin, ymax], color="0.5", lw=0.9, ls=(0,(2,2)), zorder=3)
        if ypar and ypar > 1:
            for y in np.linspace(ymin, ymax, int(ypar)+1)[1:-1]:
                ax.plot([xmin, xmax], [y, y], color="0.5", lw=0.9, ls=(0,(2,2)), zorder=3)


def add_dimension(ax, p0, p1, text, offset=(0,0), text_offset=(0,0), color="0.1"):
    x0,y0=p0; x1,y1=p1; ox,oy=offset
    ax.annotate("", xy=(x1+ox,y1+oy), xytext=(x0+ox,y0+oy), arrowprops=dict(arrowstyle="<->", lw=1.0, color=color, shrinkA=0, shrinkB=0))
    ax.text((x0+x1)/2+ox+text_offset[0], (y0+y1)/2+oy+text_offset[1], text, ha="center", va="center", fontsize=8.5, color=color, bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.85))


def arrow(ax, xy, dxy, text=None, color="#1b5ea8", scale=1.0, zorder=12):
    x, y = xy; dx, dy = dxy
    ax.annotate("", xy=(x+dx*scale, y+dy*scale), xytext=(x, y), arrowprops=dict(arrowstyle="-|>", lw=1.8, mutation_scale=12, color=color), zorder=zorder)
    if text:
        ax.text(x+dx*scale*1.05, y+dy*scale*1.05, text, fontsize=8.5, color=color, va="center", ha="left" if dx>=0 else "right", bbox=dict(fc="white", ec="none", pad=1, alpha=0.75), zorder=zorder+1)


def draw_voronoi_disk(ax, cx, cy, r, seed=1, n=75, cmap="turbo", zorder=5):
    rng = np.random.default_rng(seed)
    pts = []
    while len(pts) < n:
        p = rng.uniform(-r, r, size=2)
        if p[0]**2+p[1]**2 <= r*r:
            pts.append([cx+p[0], cy+p[1]])
    pts = np.array(pts)
    res = 360
    xs = np.linspace(cx-r, cx+r, res)
    ys = np.linspace(cy-r, cy+r, res)
    X, Y = np.meshgrid(xs, ys)
    mask = (X-cx)**2 + (Y-cy)**2 <= r*r
    D = ((X[...,None]-pts[:,0])**2 + (Y[...,None]-pts[:,1])**2)
    labels = np.argmin(D, axis=-1).astype(float)
    labels[~mask] = np.nan
    im = ax.imshow(labels, extent=[cx-r, cx+r, cy-r, cy+r], origin="lower", cmap=cmap, interpolation="nearest", zorder=zorder)
    circ = patches.Circle((cx,cy), r, facecolor="none", edgecolor="0.05", lw=1.3, zorder=zorder+2)
    ax.add_patch(circ)
    im.set_clip_path(circ)
    # Thin cell boundaries by contouring selected label changes on a coarser grid.
    ax.contour(X, Y, np.nan_to_num(labels, nan=-1), levels=np.arange(-0.5, n, 1), colors="white", linewidths=0.12, alpha=0.28, zorder=zorder+1)
    return im


def colorbar_no_ticks(fig, ax, mappable, label):
    cb = fig.colorbar(mappable, ax=ax, fraction=0.035, pad=0.025)
    cb.set_ticks([])
    cb.ax.tick_params(length=0)
    cb.outline.set_linewidth(0.7)
    cb.set_label(label, fontsize=8.5)
    return cb


def add_feature_box(ax, features: list[str], loc=(0.985, 0.02)):
    if not features:
        return
    text = "Code features\n" + "\n".join(f"- {f}" for f in features[:7])
    ax.text(loc[0], loc[1], text, transform=ax.transAxes, ha="right", va="bottom", fontsize=7.3,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.65", lw=0.8, alpha=0.92), zorder=30)


def features(ci: CaseInfo) -> list[str]:
    s = (ci.source_text + " " + ci.case_id).lower()
    out: list[str] = []
    method = None
    if isinstance(ci.pfw, dict):
        method = ci.pfw.get("updateMethod") or ci.pfw.get("updateMethodString")
    if method:
        out.append(f"{method} update")
    elif "fmpm" in s:
        out.append("FMPM update")
    elif "xpic" in s:
        out.append("XPIC transfer")
    elif "pic" in s:
        out.append("PIC transfer")
    elif "flip" in s:
        out.append("FLIP transfer")
    if "singlepointbspline" in s or "bspline" in s or "particletype=1" in s.replace(" ", ""):
        out.append("SinglePointBSpline particles")
    if "periodic" in s or "pbc" in s:
        out.append("periodic boundaries")
    if "cohesive" in s or "cz" in s:
        out.append("cohesive-zone fields")
    if "contact" in s or "gapcorrection" in s:
        out.append("contact/multi-field interaction")
    if "reaction" in s:
        out.append("reaction-history output")
    if "boxaverage" in s or "box average" in s:
        out.append("box-average output")
    if "geomechanics" in s or "ghareb" in s:
        out.append("Geomechanics material")
    if "elastic" in s:
        out.append("elastic constitutive response")
    if "vonmises" in s:
        out.append("Von Mises plasticity")
    if "anneal" in s:
        out.append("annealing event")
    if "borehole" in s:
        out.append("borehole/confining pressure events")
    if "weibull" in s:
        out.append("Weibull strength scaling")
    # preserve order and uniqueness
    uniq=[]
    for x in out:
        if x not in uniq:
            uniq.append(x)
    return uniq


def infer_partitions(ci: CaseInfo) -> tuple[int|None, int|None]:
    def v(k):
        try:
            return int(ci.pfw.get(k))
        except Exception:
            return None
    return v("xpar"), v("ypar")


def render_case(ci: CaseInfo, out_png: Path, out_svg: Path | None = None, force: bool = True) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    if out_png.exists() and not force:
        return
    family = ci.family
    if family == "brazilian_disk":
        fig = draw_brazilian(ci)
    elif family == "benchy":
        fig = draw_benchy(ci)
    elif family == "borehole":
        fig = draw_borehole(ci)
    elif family == "pdc":
        fig = draw_pdc(ci)
    elif family in ("colliding_disks", "elastic_disk"):
        fig = draw_disk_pair_or_elastic(ci)
    elif family == "spinning_disk":
        fig = draw_spinning_disk(ci)
    elif family == "expanding_ring":
        fig = draw_expanding_ring(ci)
    elif family == "expanding_bar":
        fig = draw_expanding_bar(ci)
    elif family == "cohesive_zone":
        fig = draw_cohesive(ci)
    elif family in ("contact", "contact_bc"):
        fig = draw_contact(ci)
    elif family in ("geomechanics", "ceramic", "ftable", "anneal", "cubic", "stress_control", "temperature_table", "von_mises"):
        fig = draw_specimen(ci)
    elif family in ("partition", "periodic", "momentum", "triply_periodic"):
        fig = draw_partition_periodic(ci)
    elif family in ("gas", "vortex", "implicit_fluid", "material_swap", "polymer_cz", "polymer_heal", "shrinkage", "notched_bar"):
        fig = draw_physics_special(ci)
    else:
        fig = draw_generic(ci)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.04)
    if out_svg:
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def draw_brazilian(ci: CaseInfo):
    xmin, xmax, ymin, ymax, zmin, zmax = domain_bounds(ci)
    # For the canonical examples, force the well-known normalized geometry when input
    # bounds are not available or have arbitrary units.
    if "brazilian" in ci.case_id.lower():
        xmin, xmax, ymin, ymax = -0.75, 0.75, -0.75, 0.75
    fig, ax = setup_axes(ci, (7.0, 5.2))
    w, h = xmax-xmin, ymax-ymin
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=infer_partitions(ci), fill="#fafafa")
    piston_h = 0.10*h
    ax.add_patch(patches.Rectangle((xmin, ymax-piston_h), w, piston_h, fc="#d9d2c3", ec="0.25", lw=1.1, zorder=4))
    ax.add_patch(patches.Rectangle((xmin, ymin), w, piston_h, fc="#d9d2c3", ec="0.25", lw=1.1, zorder=4))
    r = min(w, h) / 3.0
    cx, cy = (xmin+xmax)/2, (ymin+ymax)/2
    im = draw_voronoi_disk(ax, cx, cy, r, seed=7, n=85, zorder=6)
    colorbar_no_ticks(fig, ax, im, "Weibull strength scale")
    add_dimension(ax, (cx-r, cy), (cx+r, cy), "disk diameter", offset=(0, -0.18*r), text_offset=(0, -0.04*r))
    add_dimension(ax, (xmin, ymin), (xmax, ymin), f"box width = {short_num(w, ' mm') if abs(w-1.5)<1e-9 else 'domain width'}", offset=(0, -0.12*h))
    add_dimension(ax, (xmax, ymin), (xmax, ymax), f"box height = {short_num(h, ' mm') if abs(h-1.5)<1e-9 else 'domain height'}", offset=(0.08*w, 0))
    arrow(ax, (cx-0.05*w, cy+0.04*h), (0.23*w, 0), "initial velocity: vx = 1", color="#1f5aa6")
    ax.text(xmin-0.035*w, cy, "periodic", rotation=90, ha="center", va="center", fontsize=8.5)
    ax.text(xmax+0.035*w, cy, "periodic", rotation=90, ha="center", va="center", fontsize=8.5)
    ax.text(cx, ymax+0.035*h, "frictionless piston via F-table BC", ha="center", va="bottom", fontsize=8.2)
    ax.text(cx, ymin-0.065*h, "frictionless piston via F-table BC", ha="center", va="top", fontsize=8.2)
    add_title(ax, ci, "Plane-strain Brazilian disk compression; particles colored by initial strength scale")
    add_feature_box(ax, features(ci))
    ax.set_xlim(xmin-0.25*w, xmax+0.34*w); ax.set_ylim(ymin-0.22*h, ymax+0.18*h)
    return fig


def draw_benchy(ci: CaseInfo):
    fig, ax = setup_axes(ci, (7.6, 5.2))
    xmin, xmax, ymin, ymax = -60.0, 34.0, -30.0, 40.0
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=(4, 3), fill="#fbfbfb")

    # A publication-style 3DBenchy side silhouette. The boat is drawn as a
    # copper STL body rather than a generic disk; small panel lines indicate the
    # imported triangulated surface used by the STL geometry object.
    hull = np.array([
        [-39.0, -9.0], [-30.0, -15.0], [-9.0, -16.5], [13.0, -10.5],
        [22.0, -4.0], [11.0, 0.5], [-29.0, 1.8], [-43.0, -3.5],
    ])
    hull_patch = patches.Polygon(hull, closed=True, fc="#c77c35", ec="0.10", lw=1.35, zorder=6)
    ax.add_patch(hull_patch)
    deck = np.array([[-30.0, 1.3], [13.0, 0.5], [8.0, 3.6], [-25.0, 4.4]])
    ax.add_patch(patches.Polygon(deck, closed=True, fc="#d99a57", ec="0.12", lw=0.9, zorder=7))
    cabin = np.array([[-17.0, 4.0], [2.5, 3.8], [5.2, 19.0], [-10.8, 22.0], [-20.0, 13.2]])
    ax.add_patch(patches.Polygon(cabin, closed=True, fc="#dfa965", ec="0.12", lw=1.0, zorder=8))
    bow_stack = np.array([[3.5, 5.0], [9.5, 5.0], [10.5, 12.0], [4.5, 12.0]])
    ax.add_patch(patches.Polygon(bow_stack, closed=True, fc="#b86b2e", ec="0.12", lw=0.85, zorder=8))
    ax.add_patch(patches.Rectangle((-24.0, 2.8), 8.5, 11.2, fc="#b86b2e", ec="0.12", lw=0.85, zorder=8))
    for x, y, r in [(-12.0, 10.2, 2.7), (-20.0, -6.5, 2.2), (-2.0, -7.4, 2.2)]:
        ax.add_patch(patches.Circle((x, y), r, fc="white", ec="0.18", lw=0.75, zorder=9))

    # STL facet hints clipped to the hull/deck region.
    rng = np.random.default_rng(13)
    for _ in range(48):
        x0 = rng.uniform(-38, 18)
        y0 = rng.uniform(-14, 3)
        x1 = x0 + rng.uniform(3, 8)
        y1 = y0 + rng.uniform(-1.8, 1.8)
        line, = ax.plot([x0, x1], [y0, y1], color="#8a4c24", lw=0.28, alpha=0.35, zorder=9)
        line.set_clip_path(hull_patch)

    # Steel impactor in its own particle region/contact group.
    ball_center = (-49.5, 8.5)
    ball_r = 6.5
    ax.add_patch(patches.Circle(ball_center, ball_r, fc="#7f8b96", ec="0.08", lw=1.2, zorder=10))
    ax.add_patch(patches.Circle((ball_center[0] - 2.2, ball_center[1] + 2.1), 1.7, fc="#b7c0c8", ec="none", alpha=0.9, zorder=11))
    arrow(ax, (ball_center[0] - 10.0, ball_center[1]), (18.0, 0.0), "steel ball initial velocity", color="#1f5aa6", scale=1.0, zorder=12)

    # Contact and region labels.
    ax.plot([-40, 22], [-17.8, -17.8], color="0.25", lw=0.8, ls=(0, (3, 2)), zorder=5)
    ax.text(-9, -22.0, "copper 3DBenchy STL\nParticleRegion1, contact group 0", ha="center", va="top", fontsize=8.3)
    ax.text(ball_center[0], ball_center[1] + ball_r + 4.0, "steel sphere\nParticleRegion2, contact group 1", ha="center", va="bottom", fontsize=8.3)
    ax.text(xmin + 2, ymax - 3.0, "background grid", fontsize=8.2, ha="left", va="top", color="0.25")
    ax.text(17.0, -3.0, "prescribed\nmulti-field contact", ha="center", va="center", fontsize=8.0,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.65", lw=0.7, alpha=0.92), zorder=20)

    add_title(ax, ci, "STL import: copper Benchy boat impacted by a steel ball")
    add_feature_box(ax, features(ci))
    ax.set_xlim(xmin - 6, xmax + 6)
    ax.set_ylim(ymin - 5, ymax + 8)
    return fig

def draw_borehole(ci: CaseInfo):
    fig, ax = setup_axes(ci, (6.7, 5.2))
    xmin, xmax, ymin, ymax = -1, 1, -1, 1
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=(3,3), fill="#f8f5ef")
    r = 0.22
    rock = patches.Circle((0,0), 0.92, fc="#cfc7b6", ec="0.1", lw=1.0, alpha=0.92, zorder=4)
    ax.add_patch(rock)
    hole = patches.Circle((0,0), r, fc="white", ec="0.05", lw=1.4, zorder=6)
    ax.add_patch(hole)
    # pressure arrows
    for th in np.linspace(0, 2*np.pi, 12, endpoint=False):
        x0, y0 = 0.88*np.cos(th), 0.88*np.sin(th)
        arrow(ax, (x0,y0), (-0.13*np.cos(th), -0.13*np.sin(th)), None, color="#7b4b27", scale=1.0, zorder=9)
    for th in np.linspace(0, 2*np.pi, 8, endpoint=False):
        x0, y0 = r*np.cos(th), r*np.sin(th)
        arrow(ax, (x0,y0), (0.11*np.cos(th), 0.11*np.sin(th)), None, color="#1d6f8f", scale=1.0, zorder=9)
    ax.text(0, -0.04, "borehole", ha="center", va="center", fontsize=8)
    ax.text(0.62, 0.63, "confining\npressure", ha="center", va="center", fontsize=8, color="#7b4b27")
    ax.text(-0.48, -0.26, "borehole\npressure", ha="center", va="center", fontsize=8, color="#1d6f8f")
    ax.text(xmin+0.04, ymin+0.04, "plane strain\nGhareb material", ha="left", va="bottom", fontsize=8)
    add_title(ax, ci, "Plane-strain borehole collapse/cavity expansion")
    add_feature_box(ax, features(ci))
    ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.18, 1.18)
    return fig


def draw_pdc(ci: CaseInfo):
    fig, ax = setup_axes(ci, (7.4, 5.2))
    xmin, xmax, ymin, ymax = 0, 40, 0, 30
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=(4,3), fill="#fbfbfb")
    # substrate and cutter in perspective-style 2D
    substrate = np.array([[3,5],[36,5],[36,19],[3,19]])
    ax.add_patch(patches.Polygon(substrate, fc="#244c9e", ec="0.08", lw=1.2, zorder=4, alpha=0.96))
    ax.plot([3,36],[19,19], color="0.1", lw=1.0)
    cutter = np.array([[17,18],[26,18],[29,25],[20,26]])
    ax.add_patch(patches.Polygon(cutter, fc="#2b78c6", ec="0.08", lw=1.1, zorder=7))
    contact = np.array([[16.5,17.4],[28,17.8],[29,19.3],[17.2,19.0]])
    ax.add_patch(patches.Polygon(contact, fc="#ffba3b", ec="0.12", lw=0.6, zorder=8, alpha=0.95))
    arrow(ax, (24,25), (-3.5,-6.0), "prescribed cutter motion", color="#1f5aa6")
    ax.text(8,13, "substrate", color="white", fontsize=9, ha="center", va="center")
    ax.text(24,27, "PDC cutter", fontsize=8.5, ha="center", va="bottom")
    ax.text(1.5,15, "periodic", rotation=90, fontsize=8, ha="center", va="center")
    ax.text(38.7,15, "periodic", rotation=90, fontsize=8, ha="center", va="center")
    add_title(ax, ci, "PDC cutter/substrate contact with prescribed motion")
    add_feature_box(ax, features(ci))
    ax.set_xlim(-3, 43); ax.set_ylim(-2, 34)
    return fig


def draw_disk_pair_or_elastic(ci: CaseInfo):
    fig, ax = setup_axes(ci, (6.9, 5.2))

    if ci.family == "colliding_disks":
        xmin, xmax, ymin, ymax, _, _ = domain_bounds(ci)
        # Use a stable view window when legacy inputs do not expose bounds cleanly.
        if not all(np.isfinite([xmin, xmax, ymin, ymax])) or (xmax - xmin) <= 0 or (ymax - ymin) <= 0:
            xmin, xmax, ymin, ymax = -1.0, 1.0, -1.0, 1.0
        w, h = xmax - xmin, ymax - ymin
        if w / max(h, 1e-12) > 3.0 or h / max(w, 1e-12) > 3.0:
            xmin, xmax, ymin, ymax = -1.0, 1.0, -1.0, 1.0
            w, h = 2.0, 2.0

        draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=infer_partitions(ci), fill="#fafafa")

        # Pull disk geometry from the PFW objects when possible. Fall back to the
        # diagonal two-disk arrangement used by the example if input execution is
        # not available during report-only regeneration.
        disks: list[dict[str, Any]] = []
        for obj in geometry_objects(ci):
            ext = object_extents(obj)
            if ext is None:
                continue
            x0, x1, y0, y1 = ext
            if not all(np.isfinite([x0, x1, y0, y1])):
                continue
            rw = 0.5 * (x1 - x0)
            rh = 0.5 * (y1 - y0)
            if rw <= 0 or rh <= 0:
                continue
            # Cylinders/spheres appear as near-circular extents in the x-y plot.
            if 0.45 <= rw / max(rh, 1e-12) <= 2.2:
                disks.append({
                    "center": ((x0 + x1) / 2.0, (y0 + y1) / 2.0),
                    "radius": 0.5 * (rw + rh),
                    "velocity": object_velocity(obj),
                })
        if len(disks) < 2:
            cx = 0.5 * (xmin + xmax)
            cy = 0.5 * (ymin + ymax)
            r = 0.16 * min(w, h)
            disks = [
                {"center": (cx - 0.28 * w, cy - 0.23 * h), "radius": r, "velocity": (1.0, 0.75)},
                {"center": (cx + 0.28 * w, cy + 0.23 * h), "radius": r, "velocity": (-1.0, -0.75)},
            ]
        else:
            # Keep the two largest disks. If velocities are missing, point each
            # disk toward the other so the schematic still shows the diagonal
            # collision implied by the geometry arrangement.
            disks = sorted(disks, key=lambda d: d["radius"], reverse=True)[:2]
            if disks[0].get("velocity") is None or disks[1].get("velocity") is None:
                c0 = np.array(disks[0]["center"], dtype=float)
                c1 = np.array(disks[1]["center"], dtype=float)
                d = c1 - c0
                if np.linalg.norm(d) < 1e-12:
                    d = np.array([1.0, 0.75])
                d = d / np.linalg.norm(d)
                disks[0]["velocity"] = tuple(d)
                disks[1]["velocity"] = tuple(-d)

        colors = ["#86b6d8", "#d9957b"]
        labels = ["body 1", "body 2"]
        for i, d in enumerate(disks[:2]):
            cx, cy = d["center"]
            r = max(d["radius"], 0.11 * min(w, h))
            ax.add_patch(patches.Circle((cx, cy), r, fc=colors[i], ec="0.08", lw=1.25, zorder=6 + i, alpha=0.96))
            ax.text(cx, cy, labels[i], ha="center", va="center", fontsize=8.2, zorder=9)
            vel = d.get("velocity") or ((-1.0) ** i, 0.75 * ((-1.0) ** i))
            vx, vy = float(vel[0]), float(vel[1])
            n = math.hypot(vx, vy)
            if n < 1e-12:
                other = np.array(disks[1 - i]["center"], dtype=float)
                here = np.array(d["center"], dtype=float)
                v = other - here
                n = float(np.linalg.norm(v)) or 1.0
                vx, vy = float(v[0] / n), float(v[1] / n)
            scale = 0.28 * min(w, h) / max(n, 1e-12)
            arrow(ax, (cx - 0.5 * vx * scale, cy - 0.5 * vy * scale), (vx * scale, vy * scale), "initial velocity" if i == 0 else None, color="#1f5aa6", zorder=12)

        ax.text(xmin - 0.045 * w, 0.5 * (ymin + ymax), "periodic", rotation=90, ha="center", va="center", fontsize=8.2)
        ax.text(xmax + 0.045 * w, 0.5 * (ymin + ymax), "periodic", rotation=90, ha="center", va="center", fontsize=8.2)
        ax.text(0.5 * (xmin + xmax), ymax + 0.055 * h, "diagonal disk collision in plane strain", ha="center", va="bottom", fontsize=8.5)
        add_title(ax, ci, "Two circular particle bodies with diagonal, opposing velocities")
        add_feature_box(ax, features(ci))
        ax.set_xlim(xmin - 0.18 * w, xmax + 0.18 * w)
        ax.set_ylim(ymin - 0.18 * h, ymax + 0.18 * h)
        return fig

    # Elastic disk compression case.
    xmin, xmax, ymin, ymax = -1.2, 1.2, -0.8, 0.8
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=(2, 1), fill="#fafafa")
    ax.add_patch(patches.Circle((0, 0), 0.42, fc="#87b7d9", ec="0.1", lw=1.2, zorder=5))
    ax.add_patch(patches.Rectangle((xmin, ymax - 0.10 * (ymax - ymin)), xmax - xmin, 0.10 * (ymax - ymin), fc="#d9d2c3", ec="0.25", lw=1.0, zorder=4))
    ax.add_patch(patches.Rectangle((xmin, ymin), xmax - xmin, 0.10 * (ymax - ymin), fc="#d9d2c3", ec="0.25", lw=1.0, zorder=4))
    arrow(ax, (0, 0.71), (0, -0.28), "F-table piston", color="#7b4b27")
    arrow(ax, (0, -0.71), (0, 0.28), "F-table piston", color="#7b4b27")
    add_title(ax, ci, "Elastic disk compression with prescribed boundary motion")
    add_feature_box(ax, features(ci))
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.0, 1.0)
    return fig

def draw_spinning_disk(ci: CaseInfo):
    fig, ax = setup_axes(ci, (6.4, 5.0))
    xmin, xmax, ymin, ymax = -1,1,-1,1
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=(2,2), fill="#fafafa")
    ax.add_patch(patches.Circle((0,0),0.62,fc="#87b7d9",ec="0.1",lw=1.2,zorder=5))
    for th in np.linspace(0,2*np.pi,8,endpoint=False):
        x,y=0.42*np.cos(th),0.42*np.sin(th)
        tx,ty=-0.18*np.sin(th),0.18*np.cos(th)
        arrow(ax,(x,y),(tx,ty),None,color="#1f5aa6")
    ax.text(0,0,"omega",fontsize=11,ha="center",va="center")
    add_title(ax, ci, "Spinning disk: tangential velocity field")
    add_feature_box(ax, features(ci))
    ax.set_xlim(-1.2,1.2); ax.set_ylim(-1.15,1.15)
    return fig


def draw_expanding_ring(ci: CaseInfo):
    fig, ax = setup_axes(ci, (6.6, 5.0))
    xmin, xmax, ymin, ymax = -1.1,1.1,-1.1,1.1
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=(2,2), fill="#fafafa")
    ax.add_patch(patches.Wedge((0,0),0.72,0,360,width=0.18,fc="#88bddd",ec="0.1",lw=1.2,zorder=5))
    for th in np.linspace(0,2*np.pi,12,endpoint=False):
        arrow(ax,(0.78*np.cos(th),0.78*np.sin(th)),(0.18*np.cos(th),0.18*np.sin(th)),None,color="#1f5aa6")
    ax.text(0,0,"annular specimen",ha="center",va="center",fontsize=8.5)
    add_title(ax, ci, "Expanding ring with prescribed radial initial velocity")
    add_feature_box(ax, features(ci))
    ax.set_xlim(-1.3,1.3); ax.set_ylim(-1.25,1.25)
    return fig


def draw_expanding_bar(ci: CaseInfo):
    fig, ax = setup_axes(ci, (7.0, 3.6))
    xmin, xmax, ymin, ymax = -2.8,2.8,-0.35,0.35
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=(4,1), fill="#fafafa")
    ax.add_patch(patches.Rectangle((-2.35,-0.16),4.7,0.32,fc="#89b9dc",ec="0.1",lw=1.2,zorder=5))
    arrow(ax,(-2.2,0),(-0.45,0),"-v0",color="#1f5aa6")
    arrow(ax,(2.2,0),(0.45,0),"v0",color="#1f5aa6")
    add_dimension(ax,(-2.35,-0.22),(2.35,-0.22),"bar length")
    add_title(ax, ci, "Expanding bar analytic comparison")
    add_feature_box(ax, features(ci), loc=(0.985,0.82))
    ax.set_xlim(-3.1,3.1); ax.set_ylim(-0.65,0.65)
    return fig


def draw_cohesive(ci: CaseInfo):
    fig, ax = setup_axes(ci, (6.8, 4.9))
    xmin, xmax, ymin, ymax = -1.0,1.0,-0.65,0.65
    draw_domain(ax, xmin, xmax, ymin, ymax, grid=True, partitions=(3,2), fill="#fafafa")
    ax.add_patch(patches.Rectangle((-0.82,-0.38),1.64,0.76,fc="#d7e8f4",ec="0.1",lw=1.2,zorder=4))
    angle = 0.0
    m = re.search(r"(7p5|15|45|90|0)deg", ci.case_id)
    if m:
        val=m.group(1).replace("p5",".5")
        with contextlib.suppress(Exception): angle=float(val)
    elif "shear" in ci.case_id.lower(): angle=0.0
    elif "bicrystal" in ci.case_id.lower(): angle=35.0
    rad=math.radians(angle)
    L=1.1; dx=L*math.cos(rad); dy=L*math.sin(rad)
    ax.plot([-dx/2,dx/2],[-dy/2,dy/2],color="#c42d28",lw=3.0,zorder=7,solid_capstyle="round")
    ax.text(0.02,0.08,"cohesive zone",color="#9b1d19",fontsize=8.5,ha="left")
    arrow(ax,(-0.95,0),(0.28,0),"load/contact",color="#1f5aa6")
    arrow(ax,(0.95,0),(-0.28,0),None,color="#1f5aa6")
    if angle:
        ax.text(0.42,0.32,f"{angle:g} deg interface",fontsize=8)
    add_title(ax, ci, "Cohesive-zone/interface verification")
    add_feature_box(ax, features(ci))
    ax.set_xlim(-1.25,1.25); ax.set_ylim(-0.9,0.9)
    return fig


def draw_contact(ci: CaseInfo):
    fig, ax = setup_axes(ci, (6.8, 4.8))
    xmin,xmax,ymin,ymax=-1,1,-0.65,0.65
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(2,1),fill="#fafafa")
    if ci.family == "contact_bc":
        ax.add_patch(patches.Circle((-0.35,0.15),0.22,fc="#8ebbd8",ec="0.1",lw=1.2,zorder=5))
        ax.add_patch(patches.Rectangle((0.25,-0.45),0.16,0.9,fc="#bbbbbb",ec="0.1",lw=1.2,zorder=5))
        arrow(ax,(-0.70,0.15),(0.28,0),"incident velocity",color="#1f5aa6")
        ax.text(0.43,0,"contact boundary",rotation=90,fontsize=8,va="center")
        subtitle="Particle body interacting with a contact boundary"
    else:
        ax.add_patch(patches.Rectangle((-0.78,-0.28),0.72,0.56,fc="#8ebbd8",ec="0.1",lw=1.2,zorder=5))
        ax.add_patch(patches.Rectangle((0.06,-0.28),0.72,0.56,fc="#d9957b",ec="0.1",lw=1.2,zorder=5))
        ax.plot([0,0],[-0.35,0.35],color="#c42d28",lw=2.5,zorder=7)
        ax.text(0.0,0.42,"contact interface",ha="center",fontsize=8.5)
        arrow(ax,(-0.95,0),(0.25,0),"v",color="#1f5aa6")
        arrow(ax,(0.95,0),(-0.25,0),None,color="#1f5aa6")
        subtitle="Contact/interface verification"
    add_title(ax,ci,subtitle)
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.2,1.2); ax.set_ylim(-0.82,0.82)
    return fig


def draw_specimen(ci: CaseInfo):
    fig, ax = setup_axes(ci, (6.8, 4.8))
    xmin,xmax,ymin,ymax=-1,1,-0.75,0.75
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(2,2),fill="#fafafa")
    specimen = patches.Rectangle((-0.45,-0.45),0.9,0.9,fc="#d7e8f4",ec="0.1",lw=1.2,zorder=5)
    ax.add_patch(specimen)
    subtitle="Material-point specimen"
    if ci.family == "anneal":
        ax.add_patch(patches.Rectangle((-0.45,-0.45),0.9,0.9,fc="#f5a142",ec="0.1",lw=1.2,zorder=5,alpha=0.82))
        ax.text(0,0,"anneal\nregion",ha="center",va="center",fontsize=9)
        subtitle="Annealing event applied to material-point region"
    elif ci.family == "geomechanics":
        arrow(ax,(0,0.72),(0,-0.23),"sigma",color="#7b4b27")
        arrow(ax,(0,-0.72),(0,0.23),None,color="#7b4b27")
        arrow(ax,(-0.95,0),(0.25,0),None,color="#7b4b27")
        arrow(ax,(0.95,0),(-0.25,0),None,color="#7b4b27")
        ax.text(0,0,"Ghareb\ngeomechanics",ha="center",va="center",fontsize=8)
        subtitle="Constitutive verification under prescribed stress/strain path"
    elif ci.family == "ceramic":
        im = draw_voronoi_disk(ax, 0,0,0.48, seed=11, n=45, zorder=6)
        colorbar_no_ticks(fig, ax, im, "strength scale")
        arrow(ax,(0,0.72),(0,-0.23),"compression/tension",color="#7b4b27")
        subtitle="Ceramic damage specimen with strength heterogeneity"
    elif ci.family == "cubic":
        ax.text(0,0,"cubic\ncrystal",ha="center",va="center",fontsize=9)
        ax.annotate("[100]", xy=(0.48,-0.50), xytext=(0.12,-0.55), arrowprops=dict(arrowstyle="->",lw=1.1), fontsize=8)
        ax.annotate("[010]", xy=(-0.50,0.48), xytext=(-0.65,0.08), arrowprops=dict(arrowstyle="->",lw=1.1), fontsize=8)
        subtitle="Cubic-elastic orientation check"
    elif ci.family == "ftable":
        arrow(ax,(0,0.72),(0,-0.24),"F-table displacement",color="#1f5aa6")
        ax.text(0,0,"elastic block",ha="center",va="center",fontsize=9)
        subtitle="F-table interpolation boundary-value problem"
    elif ci.family == "stress_control":
        arrow(ax,(0,0.72),(0,-0.24),"stress control",color="#7b4b27")
        ax.text(0,0,"stress\npath",ha="center",va="center",fontsize=9)
        subtitle="Stress-control event verification"
    elif ci.family == "temperature_table":
        ax.add_patch(patches.Rectangle((-0.45,-0.45),0.9,0.9,fc="#f2c46d",ec="0.1",lw=1.2,zorder=5,alpha=0.9))
        ax.text(0,0,"temperature\ntable",ha="center",va="center",fontsize=9)
        subtitle="Temperature table and thermal constitutive coupling"
    elif ci.family == "von_mises":
        arrow(ax,(0,0.72),(0,-0.24),"loading",color="#7b4b27")
        ax.text(0,0,"Von Mises\nplasticity",ha="center",va="center",fontsize=9)
        subtitle="J2 plasticity verification"
    add_title(ax,ci,subtitle)
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.22,1.22); ax.set_ylim(-0.95,0.95)
    return fig


def draw_partition_periodic(ci: CaseInfo):
    fig, ax = setup_axes(ci, (7.0, 4.7))
    xmin,xmax,ymin,ymax=-1.2,1.2,-0.65,0.65
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(5 if "5x" in ci.case_id else 3, 2),fill="#fafafa")
    if "disk" in ci.case_id.lower() or "ball" in ci.case_id.lower():
        ax.add_patch(patches.Circle((-0.55,0),0.2,fc="#89b9dc",ec="0.1",lw=1.2,zorder=5))
        arrow(ax,(-0.85,0),(0.55,0),"cross-partition trajectory",color="#1f5aa6")
    else:
        ax.add_patch(patches.Rectangle((-0.65,-0.22),0.52,0.44,fc="#89b9dc",ec="0.1",lw=1.2,zorder=5))
        arrow(ax,(-0.85,0),(0.55,0),"momentum",color="#1f5aa6")
    ax.text(xmin-0.05,ymin+0.33*(ymax-ymin),"periodic",rotation=90,fontsize=8,ha="center",va="center")
    ax.text(xmax+0.05,ymin+0.33*(ymax-ymin),"periodic",rotation=90,fontsize=8,ha="center",va="center")
    ax.text(0,ymax+0.08,"processor partitions",ha="center",fontsize=8)
    add_title(ax,ci,"Partition transfer, periodicity, and conservation diagnostics")
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.42,1.42); ax.set_ylim(-0.88,0.88)
    return fig


def draw_physics_special(ci: CaseInfo):
    if ci.family == "vortex":
        return draw_vortex(ci)
    if ci.family == "notched_bar":
        return draw_notched_bar(ci)
    if ci.family in ("polymer_cz", "polymer_heal"):
        return draw_polymer(ci)
    if ci.family == "shrinkage":
        return draw_shrinkage(ci)
    if ci.family == "gas":
        return draw_gas(ci)
    if ci.family == "implicit_fluid":
        return draw_implicit_fluid(ci)
    if ci.family == "material_swap":
        return draw_material_swap(ci)
    return draw_generic(ci)


def draw_vortex(ci):
    fig, ax = setup_axes(ci,(6.6,5.1))
    xmin,xmax,ymin,ymax=-1,1,-1,1
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(3,3),fill="#fafafa")
    xs=np.linspace(-0.75,0.75,8); ys=np.linspace(-0.75,0.75,8)
    for x in xs:
        for y in ys:
            r2=x*x+y*y
            if r2<0.75**2 and r2>0.05:
                u=-y*0.08; v=x*0.08
                arrow(ax,(x,y),(u,v),None,color="#1f5aa6",scale=1.0,zorder=5)
    ax.text(0,0,"manufactured\nvortex field",ha="center",va="center",fontsize=9)
    add_title(ax,ci,"Manufactured vortex solution")
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.18,1.18); ax.set_ylim(-1.14,1.14)
    return fig


def draw_notched_bar(ci):
    fig, ax = setup_axes(ci,(7.0,3.8))
    xmin,xmax,ymin,ymax=-2.2,2.2,-0.55,0.55
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(4,1),fill="#fafafa")
    verts=[(-1.75,-0.3),(-0.18,-0.3),(0,-0.12),(0.18,-0.3),(1.75,-0.3),(1.75,0.3),(0.18,0.3),(0,0.12),(-0.18,0.3),(-1.75,0.3)]
    ax.add_patch(patches.Polygon(verts,fc="#d7e8f4",ec="0.1",lw=1.2,zorder=5))
    arrow(ax,(-2.0,0),(-0.35,0),"tension",color="#1f5aa6")
    arrow(ax,(2.0,0),(0.35,0),"tension",color="#1f5aa6")
    ax.text(0,0.42,"symmetric notch / ligament",ha="center",fontsize=8.5)
    add_title(ax,ci,"Notched-bar size-effect verification")
    add_feature_box(ax,features(ci),loc=(0.985,0.80))
    ax.set_xlim(-2.55,2.55); ax.set_ylim(-0.85,0.85)
    return fig


def draw_polymer(ci):
    fig, ax = setup_axes(ci,(6.8,4.8))
    xmin,xmax,ymin,ymax=-1,1,-0.7,0.7
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(2,1),fill="#fafafa")
    ax.add_patch(patches.Rectangle((-0.75,-0.35),1.5,0.7,fc="#efe2f1",ec="0.1",lw=1.2,zorder=4))
    rng=np.random.default_rng(3)
    nodes=rng.uniform([-0.65,-0.25],[0.65,0.25],size=(16,2))
    for i in range(len(nodes)):
        for j in range(i+1,len(nodes)):
            if np.linalg.norm(nodes[i]-nodes[j])<0.32:
                ax.plot([nodes[i,0],nodes[j,0]],[nodes[i,1],nodes[j,1]],color="#7d4a8a",lw=0.7,zorder=5,alpha=0.75)
    ax.scatter(nodes[:,0],nodes[:,1],s=8,c="#7d4a8a",zorder=6)
    ax.plot([0,0],[-0.34,0.34],color="#c42d28",lw=2.5,zorder=7)
    if ci.family=="polymer_heal":
        ax.text(0.06,0,"healing\nzone",fontsize=8,ha="left",va="center",color="#1b7f45")
    else:
        ax.text(0.06,0,"cohesive\nzone",fontsize=8,ha="left",va="center",color="#9b1d19")
    add_title(ax,ci,"Polymer cohesive/healing verification")
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.2,1.2); ax.set_ylim(-0.9,0.9)
    return fig


def draw_shrinkage(ci):
    fig, ax = setup_axes(ci,(6.4,4.8))
    xmin,xmax,ymin,ymax=-1,1,-0.8,0.8
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(2,2),fill="#fafafa")
    body=patches.Circle((0,0),0.45,fc="#c9b19a",ec="0.1",lw=1.2,zorder=5)
    ax.add_patch(body)
    for th in np.linspace(0,2*np.pi,10,endpoint=False):
        arrow(ax,(0.62*np.cos(th),0.62*np.sin(th)),(-0.18*np.cos(th),-0.18*np.sin(th)),None,color="#7b4b27")
    ax.text(0,0,"shrinkage",ha="center",va="center",fontsize=9)
    add_title(ax,ci,"Shrinkage/calcination volume-change verification")
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.2,1.2); ax.set_ylim(-1.0,1.0)
    return fig


def draw_gas(ci):
    fig, ax=setup_axes(ci,(6.5,4.8))
    xmin,xmax,ymin,ymax=-1,1,-0.75,0.75
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(2,2),fill="#f8fbff")
    rng=np.random.default_rng(8)
    pts=rng.uniform([-0.65,-0.45],[0.65,0.45],size=(30,2))
    ax.scatter(pts[:,0],pts[:,1],s=12,c="#2b78c6",zorder=5)
    for th in np.linspace(0,2*np.pi,8,endpoint=False):
        arrow(ax,(0.35*np.cos(th),0.25*np.sin(th)),(0.22*np.cos(th),0.16*np.sin(th)),None,color="#1f5aa6")
    ax.text(0,0,"gas EOS",ha="center",va="center",fontsize=9)
    add_title(ax,ci,"Gas expansion / dilating-domain verification")
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.2,1.2); ax.set_ylim(-0.95,0.95)
    return fig


def draw_implicit_fluid(ci):
    fig, ax=setup_axes(ci,(6.5,4.8))
    xmin,xmax,ymin,ymax=-1,1,-0.7,0.7
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(2,1),fill="#fbfbfb")
    ax.add_patch(patches.Rectangle((-0.65,-0.32),1.3,0.64,fc="#c9dcec",ec="0.1",lw=1.2,zorder=5))
    for x in np.linspace(-0.5,0.5,5):
        ax.add_patch(patches.Circle((x,0.0),0.055,fc="white",ec="#2b78c6",lw=0.8,zorder=6))
    arrow(ax,(-0.9,0),(0.35,0),"pressure gradient",color="#1f5aa6")
    ax.text(0,0.43,"implicit fluid pressure",ha="center",fontsize=8.5)
    add_title(ax,ci,"Implicit pore/fluid pressure coupling")
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.2,1.2); ax.set_ylim(-0.9,0.9)
    return fig


def draw_material_swap(ci):
    fig, ax=setup_axes(ci,(6.5,4.8))
    xmin,xmax,ymin,ymax=-1,1,-0.7,0.7
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=(2,1),fill="#fafafa")
    ax.add_patch(patches.Rectangle((-0.7,-0.3),0.7,0.6,fc="#89b9dc",ec="0.1",lw=1.1,zorder=5))
    ax.add_patch(patches.Rectangle((0.0,-0.3),0.7,0.6,fc="#d9957b",ec="0.1",lw=1.1,zorder=5))
    arrow(ax,(-0.15,0.45),(0.3,0),"material swap event",color="#1f5aa6")
    ax.text(-0.35,0,"material A",ha="center",va="center",fontsize=8)
    ax.text(0.35,0,"material B",ha="center",va="center",fontsize=8)
    add_title(ax,ci,"Material-swap event verification")
    add_feature_box(ax,features(ci))
    ax.set_xlim(-1.2,1.2); ax.set_ylim(-0.9,0.9)
    return fig


def draw_generic(ci: CaseInfo):
    fig, ax = setup_axes(ci, (6.8, 4.8))
    xmin,xmax,ymin,ymax,zmin,zmax = domain_bounds(ci)
    w,h=xmax-xmin,ymax-ymin
    padx,pady=0.12*w,0.12*h
    draw_domain(ax,xmin,xmax,ymin,ymax,grid=True,partitions=infer_partitions(ci),fill="#fafafa")
    objs=geometry_objects(ci)
    rendered=False
    colors=["#8ebbd8","#d9957b","#98c091","#b69bd6","#e0c772"]
    for i,obj in enumerate(objs[:6]):
        ext=object_extents(obj)
        if ext is None:
            continue
        x0,x1,y0,y1=ext
        if not all(np.isfinite([x0,x1,y0,y1])):
            continue
        if abs(x1-x0)<1e-12 or abs(y1-y0)<1e-12:
            continue
        # Generic extent patch.
        ax.add_patch(patches.Rectangle((x0,y0),x1-x0,y1-y0,fc=colors[i%len(colors)],ec="0.1",lw=1.0,alpha=0.75,zorder=5+i))
        vel=object_velocity(obj)
        if vel:
            scale=0.18*max(w,h)/(np.linalg.norm(vel)+1e-12)
            arrow(ax,((x0+x1)/2,(y0+y1)/2),(vel[0]*scale,vel[1]*scale),"v0",color="#1f5aa6")
        rendered=True
    if not rendered:
        ax.add_patch(patches.Rectangle((xmin+0.25*w,ymin+0.25*h),0.5*w,0.5*h,fc="#d7e8f4",ec="0.1",lw=1.2,zorder=5))
    add_title(ax,ci,"PFW geometry, mesh, and boundary-condition schematic")
    add_feature_box(ax,features(ci))
    ax.set_xlim(xmin-padx,xmax+padx); ax.set_ylim(ymin-pady,ymax+pady)
    return fig


def report_tex_files(report_dir: Path) -> list[Path]:
    return sorted(report_dir.glob("*.tex"))


def images_from_tex(tex: Path) -> list[Path]:
    text = tex.read_text(errors="replace")
    paths = []
    for m in re.finditer(r"(?:includegraphics|IfFileExists)(?:\[[^\]]*\])?\{([^}]*schematics/[^}]+\.png)\}", text):
        p = (tex.parent / m.group(1)).resolve()
        paths.append(p)
    return sorted(set(paths))


def best_case_for_image(img: Path, cases: dict[str, CaseInfo], suite: str, pfw_root: Path) -> CaseInfo:
    stem = img.stem
    if stem.endswith("_schematic"):
        stem = stem[:-10]

    def unique_cases() -> list[CaseInfo]:
        out: list[CaseInfo] = []
        seen: set[int] = set()
        for ci in cases.values():
            if id(ci) not in seen:
                seen.add(id(ci))
                out.append(ci)
        return out

    norm_stem = sanitize_id(stem)
    candidates = [stem, norm_stem, stem.replace("_", "__"), stem.replace("__", "_")]
    for c in candidates:
        if c in cases:
            return cases[c]

    # Hard family matches keep report schematics from falling back to the first
    # Brazilian-disk case when a report image is present but the source input is
    # not currently discoverable. This is especially important for Benchy and
    # colliding-disks schematics.
    stem_lower = stem.lower()
    preferred_family = classify_family(stem_lower, None, stem_lower)
    if preferred_family != "generic":
        family_matches = [ci for ci in unique_cases() if ci.family == preferred_family]
        if family_matches:
            # Prefer the case whose id shares the longest token overlap with the
            # requested filename. For Brazilian-disk variants this preserves the
            # FMPM/FLIP/PIC/XPIC/BSpline distinction.
            def score(ci: CaseInfo) -> tuple[int, int]:
                cid = sanitize_id(ci.case_id).lower()
                common = 0
                for tok in re.split(r"[_\W]+", stem_lower):
                    if tok and tok in cid:
                        common += len(tok)
                return common, len(cid)
            return sorted(family_matches, key=score, reverse=True)[0]
        suite_root = pfw_root / suite
        return CaseInfo(
            suite=suite,
            case_id=norm_stem,
            family=preferred_family,
            input_path=None,
            suite_root=suite_root,
            pfw_root=pfw_root,
            pfw={},
            namespace={},
            source_text=stem,
        )

    # Fuzzy: choose case whose key is contained in stem or vice versa.
    st = stem.lower()
    best = None
    best_len = -1
    for key, ci in cases.items():
        k = key.lower()
        if (k and k in st) or (st in k):
            if len(k) > best_len:
                best = ci
                best_len = len(k)
    if best is not None:
        return best

    suite_root = pfw_root / suite
    return CaseInfo(suite=suite, case_id=stem, family=classify_family(stem, None, ""), input_path=None, suite_root=suite_root, pfw_root=pfw_root, pfw={}, namespace={}, source_text="")

def default_report_dir(pfw_root: Path, suite: str) -> Path:
    if suite == "examples":
        return pfw_root / "examples/examples_suite_report"
    if suite == "verification":
        return pfw_root / "verification/verification_suite_report"
    if suite == "validation":
        return pfw_root / "validation/validation_suite_report"
    raise ValueError(suite)


def standard_output_name(report_dir: Path, ci: CaseInfo) -> Path:
    return report_dir / "schematics" / f"{sanitize_id(ci.case_id)}_schematic.png"


def compile_report(report_dir: Path) -> None:
    texs = report_tex_files(report_dir)
    if not texs:
        return
    pdflatex = shutil.which("pdflatex")
    if not pdflatex:
        print(f"[schematics] pdflatex not found; skipped compile for {report_dir}")
        return
    tex = texs[0]
    for _ in range(2):
        subprocess.run([pdflatex, "-interaction=nonstopmode", tex.name], cwd=report_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    print(f"[schematics] refreshed PDF: {report_dir / tex.with_suffix('.pdf').name}")


def generate_suite(repo: Path, suite: str, compile_pdf: bool = False, svg: bool = False, force: bool = True) -> dict[str, Any]:
    pfw_root = pfw_root_from(repo)
    report_dir = default_report_dir(pfw_root, suite)
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "schematics").mkdir(parents=True, exist_ok=True)
    cases = discover_cases(suite, pfw_root)
    # If a report already exists, use exactly the schematic file names referenced by the report.
    targets: list[tuple[CaseInfo, Path]] = []
    for tex in report_tex_files(report_dir):
        for img in images_from_tex(tex):
            ci = best_case_for_image(img, cases, suite, pfw_root)
            targets.append((ci, img))
    if not targets:
        # Otherwise generate one schematic per discovered case.
        for ci in sorted({id(c): c for c in cases.values()}.values(), key=lambda c: c.case_id):
            targets.append((ci, standard_output_name(report_dir, ci)))
    # Deduplicate targets.
    dedup: dict[Path, CaseInfo] = {}
    for ci, img in targets:
        dedup[img] = ci
    manifest=[]
    for img, ci in sorted(dedup.items(), key=lambda kv: kv[0].name):
        out_svg = img.with_suffix(".svg") if svg else None
        try:
            render_case(ci, img, out_svg=out_svg, force=force)
            print(f"[schematics] wrote {img.relative_to(repo) if img.is_relative_to(repo) else img}")
            manifest.append({"suite": suite, "case_id": ci.case_id, "family": ci.family, "image": str(img), "input": str(ci.input_path) if ci.input_path else None, "diagnostics": ci.diagnostics[:3]})
        except Exception as exc:
            print(f"[schematics] ERROR rendering {ci.case_id}: {exc}")
            traceback.print_exc(limit=3)
            manifest.append({"suite": suite, "case_id": ci.case_id, "family": ci.family, "image": str(img), "input": str(ci.input_path) if ci.input_path else None, "error": repr(exc)})
    manifest_path = report_dir / "schematics" / "publication_schematics_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    if compile_pdf:
        compile_report(report_dir)
    return {"suite": suite, "count": len(manifest), "manifest": str(manifest_path)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo-root", default=None, help="GEOS repository root. Defaults to auto-detection from cwd.")
    parser.add_argument("--suite", action="append", choices=["examples", "verification", "validation"], help="Suite to render. May be repeated. Default: examples and verification.")
    parser.add_argument("--compile", action="store_true", help="Run pdflatex twice after replacing schematic PNGs.")
    parser.add_argument("--svg", action="store_true", help="Also write vector SVG copies beside PNGs.")
    parser.add_argument("--no-force", action="store_true", help="Do not overwrite existing schematic files.")
    args = parser.parse_args(argv)
    repo = Path(args.repo_root).expanduser().resolve() if args.repo_root else repo_root_from(Path.cwd())
    suites = args.suite or ["examples", "verification"]
    summaries=[]
    for suite in suites:
        summaries.append(generate_suite(repo, suite, compile_pdf=args.compile, svg=args.svg, force=not args.no_force))
    print(json.dumps(summaries, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
