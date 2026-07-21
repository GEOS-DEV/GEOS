#!/usr/bin/env python3
"""Generate a 2D rigid-body MPM jamming verification particle file.

The geometry intentionally mixes analytic disks and axis-aligned blocks.  Object
colors are distinct from contact group; all objects share contact group 0 so the
rigid-body event must partition local nodal fields from ParticleColor.

Surface particles carry exact outward normals and exact surface-position vectors
from the particle point to the analytic object boundary.  Interior particles keep
zero normals/positions so they do not contribute to explicit contact geometry.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

DX = 0.0125
Z = 0.0
THICKNESS = 0.1
VOLUME = DX * DX * THICKNESS
SURFACE_BAND = 0.75 * DX
OUT = Path(__file__).resolve().parent
PARTICLE_FILE = OUT / "mpmParticleFile_rigidBodyJamming2D"
MANIFEST_FILE = OUT / "rigidBodyJamming2D_manifest.json"

SURFACE_FLAG_INTERIOR = 0
SURFACE_FLAG_SURFACE = 2
PARTICLE_TYPE_SINGLE_POINT = 0
CONTACT_GROUP = 0
MATERIAL_COPPER = 0

OBJECTS = [
    {"kind": "disk", "color": 0, "center": (-0.34, -0.20), "radius": 0.105},
    {"kind": "block", "color": 1, "center": (-0.06, -0.19), "width": 0.18, "height": 0.115},
    {"kind": "disk", "color": 2, "center": (0.24, -0.18), "radius": 0.090},
    {"kind": "block", "color": 3, "center": (-0.26, 0.14), "width": 0.16, "height": 0.155},
    {"kind": "disk", "color": 4, "center": (0.05, 0.12), "radius": 0.115},
    {"kind": "block", "color": 5, "center": (0.34, 0.14), "width": 0.145, "height": 0.135},
]

HEADER = [
    "ID", "PositionX", "PositionY", "PositionZ",
    "VelocityX", "VelocityY", "VelocityZ",
    "MaterialType", "ParticleType", "ContactGroup", "ParticleColor",
    "SurfaceFlag", "Damage", "Porosity", "Temperature", "TemperatureRate",
    "StrengthScale", "CZTag",
    "RVectorXX", "RVectorXY", "RVectorXZ",
    "RVectorYX", "RVectorYY", "RVectorYZ",
    "RVectorZX", "RVectorZY", "RVectorZZ",
    "MaterialDirectionXX", "MaterialDirectionXY", "MaterialDirectionXZ",
    "MaterialDirectionYX", "MaterialDirectionYY", "MaterialDirectionYZ",
    "MaterialDirectionZX", "MaterialDirectionZY", "MaterialDirectionZZ",
    "SurfaceNormalX", "SurfaceNormalY", "SurfaceNormalZ",
    "SurfacePositionX", "SurfacePositionY", "SurfacePositionZ",
    "SurfaceTractionX", "SurfaceTractionY", "SurfaceTractionZ",
]


def disk_contains(obj: dict, x: float, y: float) -> bool:
    cx, cy = obj["center"]
    return (x - cx) ** 2 + (y - cy) ** 2 <= obj["radius"] ** 2 + 1.0e-14


def block_contains(obj: dict, x: float, y: float) -> bool:
    cx, cy = obj["center"]
    return abs(x - cx) <= 0.5 * obj["width"] + 1.0e-14 and abs(y - cy) <= 0.5 * obj["height"] + 1.0e-14


def disk_surface(obj: dict, x: float, y: float):
    cx, cy = obj["center"]
    radius = obj["radius"]
    rx = x - cx
    ry = y - cy
    r = math.hypot(rx, ry)
    if r < 1.0e-14:
        nx, ny = 1.0, 0.0
        dist_to_surface = radius
    else:
        nx, ny = rx / r, ry / r
        dist_to_surface = radius - r
    sx = dist_to_surface * nx
    sy = dist_to_surface * ny
    is_surface = dist_to_surface <= SURFACE_BAND
    return is_surface, nx, ny, sx, sy


def block_surface(obj: dict, x: float, y: float):
    cx, cy = obj["center"]
    half_w = 0.5 * obj["width"]
    half_h = 0.5 * obj["height"]
    dx_right = cx + half_w - x
    dx_left = x - (cx - half_w)
    dy_top = cy + half_h - y
    dy_bottom = y - (cy - half_h)
    candidates = [
        (dx_right, 1.0, 0.0, dx_right, 0.0),
        (dx_left, -1.0, 0.0, -dx_left, 0.0),
        (dy_top, 0.0, 1.0, 0.0, dy_top),
        (dy_bottom, 0.0, -1.0, 0.0, -dy_bottom),
    ]
    dist, nx, ny, sx, sy = min(candidates, key=lambda item: item[0])
    is_surface = dist <= SURFACE_BAND
    return is_surface, nx, ny, sx, sy


def object_contains(obj: dict, x: float, y: float) -> bool:
    if obj["kind"] == "disk":
        return disk_contains(obj, x, y)
    return block_contains(obj, x, y)


def object_surface(obj: dict, x: float, y: float):
    if obj["kind"] == "disk":
        return disk_surface(obj, x, y)
    return block_surface(obj, x, y)


def representative_radius(obj: dict) -> float:
    if obj["kind"] == "disk":
        return obj["radius"]
    return 0.5 * math.hypot(obj["width"], obj["height"])


def bounding_box(obj: dict):
    cx, cy = obj["center"]
    if obj["kind"] == "disk":
        r = obj["radius"]
        return cx - r, cx + r, cy - r, cy + r
    return cx - 0.5 * obj["width"], cx + 0.5 * obj["width"], cy - 0.5 * obj["height"], cy + 0.5 * obj["height"]


def check_non_overlap() -> list[dict]:
    diagnostics = []
    for i, a in enumerate(OBJECTS):
        for b in OBJECTS[i + 1:]:
            ax, ay = a["center"]
            bx, by = b["center"]
            center_distance = math.hypot(ax - bx, ay - by)
            # Conservative separating-circle check for a cheap documented margin.
            margin = center_distance - representative_radius(a) - representative_radius(b)
            diagnostics.append({
                "colorA": a["color"],
                "colorB": b["color"],
                "conservativeCircularMargin": margin,
            })
            if margin <= 2.0 * DX:
                raise RuntimeError(
                    f"Objects {a['color']} and {b['color']} have insufficient conservative margin: {margin:g}"
                )
    return diagnostics


def generate_particles():
    rows = []
    particle_id = 1
    counts = {}
    surface_counts = {}
    for obj in OBJECTS:
        xmin, xmax, ymin, ymax = bounding_box(obj)
        ix0 = math.floor(xmin / DX) - 1
        ix1 = math.ceil(xmax / DX) + 1
        iy0 = math.floor(ymin / DX) - 1
        iy1 = math.ceil(ymax / DX) + 1
        counts[obj["color"]] = 0
        surface_counts[obj["color"]] = 0
        for ix in range(ix0, ix1 + 1):
            x = (ix + 0.5) * DX
            for iy in range(iy0, iy1 + 1):
                y = (iy + 0.5) * DX
                if not object_contains(obj, x, y):
                    continue
                is_surface, nx, ny, sx, sy = object_surface(obj, x, y)
                if not is_surface:
                    nx = ny = sx = sy = 0.0
                surface_flag = SURFACE_FLAG_SURFACE if is_surface else SURFACE_FLAG_INTERIOR
                counts[obj["color"]] += 1
                surface_counts[obj["color"]] += int(is_surface)
                row = [
                    particle_id, x, y, Z,
                    0.0, 0.0, 0.0,
                    MATERIAL_COPPER, PARTICLE_TYPE_SINGLE_POINT, CONTACT_GROUP, obj["color"],
                    surface_flag, 0.0, 0.0, 300.0, 0.0,
                    1.0, 0,
                    VOLUME, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0,
                    nx, ny, 0.0,
                    sx, sy, 0.0,
                    0.0, 0.0, 0.0,
                ]
                rows.append(row)
                particle_id += 1
    return rows, counts, surface_counts


def main():
    margins = check_non_overlap()
    rows, counts, surface_counts = generate_particles()
    with PARTICLE_FILE.open("w", encoding="utf-8") as f:
        f.write("\t".join(HEADER) + "\n")
        for row in rows:
            f.write(" ".join(f"{value:.17g}" if isinstance(value, float) else str(value) for value in row) + "\n")
    manifest = {
        "particleFile": PARTICLE_FILE.name,
        "dx": DX,
        "thickness": THICKNESS,
        "particleVolume": VOLUME,
        "objects": OBJECTS,
        "particleCountsByColor": counts,
        "surfaceParticleCountsByColor": surface_counts,
        "nonOverlapDiagnostics": margins,
        "notes": [
            "SurfaceFlag=2 marks analytic surface particles; interior particles carry zero normal/position.",
            "All bodies use ContactGroup=0 and distinct ParticleColor values to exercise color partitioning.",
            "Copper yield stress in the XML is 70 MPa; the continuum stress-control target is -140 MPa in x and y.",
        ],
    }
    MANIFEST_FILE.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"wrote {len(rows)} particles to {PARTICLE_FILE}")
    print(f"wrote manifest to {MANIFEST_FILE}")


if __name__ == "__main__":
    main()
