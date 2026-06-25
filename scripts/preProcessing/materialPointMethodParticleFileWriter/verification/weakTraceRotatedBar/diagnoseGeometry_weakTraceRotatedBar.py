#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import os
from pathlib import Path
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PFW_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(PFW_ROOT))


# Geometry diagnostics are often run from a lightweight Python environment.
# Provide a minimal mpi4py stub so importing pfw_geometryObjects works even
# outside the MPI-enabled PFW runtime used for actual particle generation.
try:
    import mpi4py  # noqa: F401
except Exception:
    import types
    class _Comm:
        def Get_rank(self): return 0
        def Get_size(self): return 1
        def Barrier(self): return None
    mpi4py_stub = types.ModuleType("mpi4py")
    mpi4py_stub.MPI = types.SimpleNamespace(COMM_WORLD=_Comm())
    sys.modules.setdefault("mpi4py", mpi4py_stub)


def load_input(variant: str):
    os.environ["WEAK_TRACE_VARIANT"] = variant
    os.environ.setdefault("WEAK_TRACE_CASE_NAME", f"weakTraceRotatedBar_{variant}_geometry")
    spec = importlib.util.spec_from_file_location("_weak_trace_input", SCRIPT_DIR / "pfw_input_weakTraceRotatedBar.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod.pfw


def object_expected_normal(name: str):
    a = np.array([1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0), 0.0])
    low = name.lower()
    if "lower" in low or "compliant" in low:
        return a
    if "upper" in low or "stiff" in low:
        return -a
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Diagnose weakTraceRotatedBar PFW surface geometry")
    parser.add_argument("--variant", default=os.environ.get("WEAK_TRACE_VARIANT", "traceContactGroups"))
    parser.add_argument("--output-dir", default=str(SCRIPT_DIR / "output" / "geometryDiagnostics"))
    parser.add_argument("--sample-stride", type=int, default=int(os.environ.get("WEAK_TRACE_GEOMETRY_SAMPLE_STRIDE", "1")))
    args = parser.parse_args()

    pfw = load_input(args.variant)
    out_dir = Path(args.output_dir) / args.variant
    out_dir.mkdir(parents=True, exist_ok=True)

    surface_depth = max((pfw["xmax"] - pfw["xmin"]) / (pfw["nI"] * pfw["ppc"]), (pfw["ymax"] - pfw["ymin"]) / (pfw["nJ"] * pfw["ppc"]))
    ni = int(pfw["nI"] * pfw["ppc"])
    nj = int(pfw["nJ"] * pfw["ppc"])
    dx = (pfw["xmax"] - pfw["xmin"]) / float(ni)
    dy = (pfw["ymax"] - pfw["ymin"]) / float(nj)
    z = 0.0

    rows = []
    bad = []
    counts = {}
    tangent = np.array([-1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0), 0.0])

    for ii in range(1, ni + 1, max(1, args.sample_stride)):
        x = pfw["xmin"] + (ii - 0.5) * dx
        for jj in range(1, nj + 1, max(1, args.sample_stride)):
            y = pfw["ymin"] + (jj - 0.5) * dy
            pt = np.array([x, y, z])
            for obj in pfw["objects"]:
                flag = int(obj.isInterior(pt, surface_depth))
                if flag < 0:
                    continue
                mat = int(obj.getMat(pt) if hasattr(obj, "getMat") else obj.mat)
                group = int(obj.getGroup(pt) if hasattr(obj, "getGroup") else getattr(obj, "group", 0))
                n = np.zeros(3)
                s = np.zeros(3)
                if flag != 0:
                    n = np.array(obj.getSurfaceNormal(pt), dtype=float)
                    s = np.array(obj.getSurfacePosition(pt), dtype=float)
                nn = float(np.linalg.norm(n))
                ss = float(np.linalg.norm(s))
                expn = object_expected_normal(getattr(obj, "name", ""))
                dot = ""
                if flag != 0 and expn is not None:
                    dot_val = float(np.dot(n, expn))
                    dot = dot_val
                    tangential_error = float(np.dot(s, tangent))
                    is_bad = nn < 0.9 or abs(dot_val - 1.0) > 1.0e-8 or abs(tangential_error) > max(1.0e-10, 1.0e-8 * max(1.0, ss))
                    if is_bad:
                        bad.append((x, y, getattr(obj, "name", ""), flag, nn, dot_val, tangential_error, n.tolist(), s.tolist()))
                else:
                    tangential_error = ""
                key = (getattr(obj, "name", ""), mat, group, flag)
                counts[key] = counts.get(key, 0) + 1
                rows.append({
                    "x": x, "y": y, "z": z, "object": getattr(obj, "name", ""), "material": mat, "group": group, "surfaceFlag": flag,
                    "surfaceNormalX": n[0], "surfaceNormalY": n[1], "surfaceNormalZ": n[2], "surfaceNormalNorm": nn,
                    "surfacePositionX": s[0], "surfacePositionY": s[1], "surfacePositionZ": s[2], "surfacePositionNorm": ss,
                    "normalDotExpected": dot, "surfacePositionTangentialError": tangential_error,
                })
                break

    csv_path = out_dir / f"weakTraceRotatedBar_{args.variant}_geometry_particles.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else ["x"])
        writer.writeheader()
        writer.writerows(rows)

    surface_rows = [r for r in rows if int(r["surfaceFlag"]) != 0]
    log = []
    log.append(f"variant: {args.variant}")
    log.append(f"particleLikeSampleCount: {len(rows)}")
    log.append(f"surfaceFlaggedParticleCount: {len(surface_rows)}")
    log.append(f"badSurfaceGeometryCount: {len(bad)}")
    log.append("countsByObjectMaterialGroupFlag:")
    for key, val in sorted(counts.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3])):
        log.append(f"  {key}: {val}")
    if surface_rows:
        norms = [float(r["surfaceNormalNorm"]) for r in surface_rows]
        dots = [float(r["normalDotExpected"]) for r in surface_rows if r["normalDotExpected"] != ""]
        tang = [abs(float(r["surfacePositionTangentialError"])) for r in surface_rows if r["surfacePositionTangentialError"] != ""]
        log.append(f"normalNorm: min={min(norms):.8e} mean={sum(norms)/len(norms):.8e} max={max(norms):.8e}")
        if dots:
            log.append(f"normalDotExpected: min={min(dots):.8e} mean={sum(dots)/len(dots):.8e} max={max(dots):.8e}")
        if tang:
            log.append(f"surfacePositionTangentialErrorAbs: max={max(tang):.8e}")
    if bad:
        log.append("firstBadRows:")
        for b in bad[:20]:
            log.append("  " + repr(b))

    (out_dir / f"weakTraceRotatedBar_{args.variant}_geometry_summary.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xs = [r["x"] for r in rows]
        ys = [r["y"] for r in rows]
        flags = [r["surfaceFlag"] for r in rows]
        fig, ax = plt.subplots(figsize=(6.0, 5.0))
        sc = ax.scatter(xs, ys, c=flags, s=8, marker="s")
        fig.colorbar(sc, ax=ax, label="surfaceFlag")
        ax.set_aspect("equal")
        ax.set_title(f"{args.variant}: PFW surface flags")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.tight_layout()
        fig.savefig(out_dir / f"weakTraceRotatedBar_{args.variant}_geometry_flags.png", dpi=180)
        plt.close(fig)
    except Exception as exc:
        (out_dir / f"weakTraceRotatedBar_{args.variant}_geometry_plot_error.txt").write_text(str(exc) + "\n")

    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
