#!/usr/bin/env python3
"""Lightweight checks for a generated rigidBodyJamming2D PFW particle file."""
from __future__ import annotations

import argparse
import math
from pathlib import Path

EXPECTED_COLORS = set(range(6))


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("particle_file", nargs="?", default="mpmParticleFile_rigidBodyJamming2D")
    return p.parse_args()


def rows_from_whitespace_table(path: Path):
    lines = path.read_text().splitlines()
    if not lines:
        return [], []
    names = lines[0].split()
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        values = line.split()
        if len(values) != len(names):
            raise SystemExit(
                f"particle row has {len(values)} values but header has {len(names)} columns: {line[:100]}"
            )
        rows.append(dict(zip(names, values)))
    return names, rows


def as_float(row, name):
    return float(row.get(name, 0.0))


def main() -> int:
    args = parse_args()
    path = Path(args.particle_file)
    if not path.is_file():
        raise SystemExit(f"missing particle file: {path}")

    names, rows = rows_from_whitespace_table(path)
    for required in ("ParticleColor", "ContactGroup", "SurfaceFlag", "SurfaceNormalX", "SurfacePositionX"):
        if required not in names:
            raise SystemExit(f"missing required particle column {required}")

    colors = set()
    counts = {c: 0 for c in EXPECTED_COLORS}
    surface_counts = {c: 0 for c in EXPECTED_COLORS}
    bad_normals = 0
    bad_interior_surface_data = 0

    for row in rows:
        c = int(float(row["ParticleColor"]))
        colors.add(c)
        if c in counts:
            counts[c] += 1
        if int(float(row["ContactGroup"])) != 0:
            raise SystemExit("all particles should share ContactGroup=0")
        flag = int(float(row["SurfaceFlag"]))
        normal = [as_float(row, f"SurfaceNormal{axis}") for axis in "XYZ"]
        spos = [as_float(row, f"SurfacePosition{axis}") for axis in "XYZ"]
        nrm = math.sqrt(sum(v * v for v in normal))
        if flag != 0:
            surface_counts[c] = surface_counts.get(c, 0) + 1
            if abs(nrm - 1.0) > 1.0e-8:
                bad_normals += 1
        elif nrm > 1.0e-14 or math.sqrt(sum(v * v for v in spos)) > 1.0e-14:
            bad_interior_surface_data += 1

    if colors != EXPECTED_COLORS:
        raise SystemExit(f"expected colors {sorted(EXPECTED_COLORS)}, observed {sorted(colors)}")
    if any(counts[c] <= 0 for c in EXPECTED_COLORS):
        raise SystemExit(f"empty object color found: {counts}")
    if any(surface_counts.get(c, 0) <= 0 for c in EXPECTED_COLORS):
        raise SystemExit(f"missing surface particles for at least one color: {surface_counts}")
    if bad_normals:
        raise SystemExit(f"surface particles with non-unit normals: {bad_normals}")
    if bad_interior_surface_data:
        raise SystemExit(f"interior particles with nonzero explicit surface data: {bad_interior_surface_data}")

    print("rigidBodyJamming2D particle-file checks passed")
    print("particle counts by color:", counts)
    print("surface counts by color:", surface_counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
