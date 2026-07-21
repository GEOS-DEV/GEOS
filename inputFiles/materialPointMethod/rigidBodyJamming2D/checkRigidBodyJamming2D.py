#!/usr/bin/env python3
"""Verification checks for the rigid-body MPM jamming example.

By default this script performs deterministic pre-run checks on the generated
particle file and performs post-run checks only when GEOS output files exist.
Pass --require-run-output to make missing GEOS history files fail the test.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PARTICLE_FILE = ROOT / "mpmParticleFile_rigidBodyJamming2D"
MANIFEST_FILE = ROOT / "rigidBodyJamming2D_manifest.json"
RIGID_HISTORY = ROOT / "rigidBodyHistory.csv"
BOX_HISTORY = ROOT / "boxAverageHistory.csv"
YIELD_STRESS = 7.0e7
TARGET_STRESS = -2.0 * YIELD_STRESS


def load_manifest() -> dict:
    if not MANIFEST_FILE.exists():
        raise FileNotFoundError(f"Missing {MANIFEST_FILE.name}; run generateRigidBodyJamming2DParticles.py first")
    return json.loads(MANIFEST_FILE.read_text(encoding="utf-8"))


def read_particles():
    if not PARTICLE_FILE.exists():
        raise FileNotFoundError(f"Missing {PARTICLE_FILE.name}; run generateRigidBodyJamming2DParticles.py first")
    with PARTICLE_FILE.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        indices = {name: i for i, name in enumerate(header)}
        rows = []
        for line_number, line in enumerate(f, start=2):
            if not line.strip():
                continue
            values = line.split()
            if len(values) != len(header):
                raise AssertionError(f"line {line_number} has {len(values)} values; expected {len(header)}")
            rows.append({name: float(values[i]) for name, i in indices.items()})
    return rows


def check_geometry(manifest: dict, rows: list[dict]):
    counts = {int(k): int(v) for k, v in manifest["particleCountsByColor"].items()}
    colors = sorted(counts)
    observed = {c: 0 for c in colors}
    surface = {c: 0 for c in colors}
    for row in rows:
        color = int(row["ParticleColor"])
        if color not in observed:
            raise AssertionError(f"unexpected ParticleColor {color}")
        observed[color] += 1
        if int(row["ContactGroup"]) != 0:
            raise AssertionError("verification expects all particles in ContactGroup=0")
        if int(row["MaterialType"]) != 0:
            raise AssertionError("verification expects all particles to use material index 0/copper")
        normal = (row["SurfaceNormalX"], row["SurfaceNormalY"], row["SurfaceNormalZ"])
        spos = (row["SurfacePositionX"], row["SurfacePositionY"], row["SurfacePositionZ"])
        normal_norm = math.sqrt(sum(x * x for x in normal))
        if int(row["SurfaceFlag"]) == 2:
            surface[color] += 1
            if abs(normal_norm - 1.0) > 1.0e-10:
                raise AssertionError(f"surface particle has non-unit normal: {normal_norm}")
            if normal[0] * spos[0] + normal[1] * spos[1] + normal[2] * spos[2] < -1.0e-12:
                raise AssertionError("surface-position vector points inward relative to the surface normal")
        else:
            if normal_norm > 1.0e-14 or math.sqrt(sum(x * x for x in spos)) > 1.0e-14:
                raise AssertionError("interior particle should not contribute explicit contact normal/position")
    if observed != counts:
        raise AssertionError(f"particle counts mismatch: observed {observed}, manifest {counts}")

    for diagnostic in manifest["nonOverlapDiagnostics"]:
        if diagnostic["conservativeCircularMargin"] <= 2.0 * manifest["dx"]:
            raise AssertionError(f"insufficient initial non-overlap margin: {diagnostic}")

    if any(value == 0 for value in surface.values()):
        raise AssertionError(f"each color must have surface particles; observed {surface}")

    print(f"geometry: {len(rows)} particles across colors {colors}; all initial non-overlap and surface checks passed")


def check_rigid_history(require: bool):
    if not RIGID_HISTORY.exists():
        if require:
            raise AssertionError(f"missing required {RIGID_HISTORY.name}")
        print(f"post-run: {RIGID_HISTORY.name} not present; skipping rigid history checks")
        return
    with RIGID_HISTORY.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise AssertionError("rigidBodyHistory.csv is empty")
    colors = {int(float(row["color"])) for row in rows}
    if colors != set(range(6)):
        raise AssertionError(f"rigidBodyHistory.csv should contain colors 0..5, observed {sorted(colors)}")
    max_force = 0.0
    max_ke = 0.0
    for row in rows:
        fx = float(row["forceX"])
        fy = float(row["forceY"])
        fz = float(row["forceZ"])
        max_force = max(max_force, math.sqrt(fx * fx + fy * fy + fz * fz))
        max_ke = max(max_ke, float(row["kineticEnergy"]))
    if max_force <= 0.0:
        raise AssertionError("rigid history did not record a nonzero resultant force")
    print(f"post-run: rigid history present; max resultant body force = {max_force:.6e}, max KE = {max_ke:.6e}")


def parse_box_history():
    with BOX_HISTORY.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        rows = list(reader)
    if not rows:
        raise AssertionError("boxAverageHistory.csv is empty")
    return rows


def check_box_history(require: bool):
    if not BOX_HISTORY.exists():
        if require:
            raise AssertionError(f"missing required {BOX_HISTORY.name}")
        print(f"post-run: {BOX_HISTORY.name} not present; skipping continuum stress checks")
        return
    rows = parse_box_history()
    final = rows[-1]
    sxx = float(final["Sxx"])
    syy = float(final["Syy"])
    tolerance = 0.30 * abs(TARGET_STRESS)
    if abs(sxx - TARGET_STRESS) > tolerance or abs(syy - TARGET_STRESS) > tolerance:
        raise AssertionError(
            f"final box stresses are not within 30% of target {TARGET_STRESS:.6e}: Sxx={sxx:.6e}, Syy={syy:.6e}"
        )
    print(f"post-run: final continuum stress check passed; Sxx={sxx:.6e}, Syy={syy:.6e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-run-output", action="store_true", help="fail if GEOS history files are absent")
    args = parser.parse_args()
    manifest = load_manifest()
    rows = read_particles()
    check_geometry(manifest, rows)
    check_rigid_history(args.require_run_output)
    check_box_history(args.require_run_output)
    print("rigid-body MPM jamming verification checks passed")


if __name__ == "__main__":
    main()
