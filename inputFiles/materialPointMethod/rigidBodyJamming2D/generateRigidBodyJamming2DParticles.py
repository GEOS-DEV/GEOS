#!/usr/bin/env python3
"""Generate the MPI/periodic rigidBodyJamming2D verification particle file.

The case uses a 2 x 3 plane-strain decomposition, periodic x, and moving y walls.
Six colors define six rigid bodies while every particle uses catch-all
ContactGroup=0. One disk straddles the x-periodic seam, so the generated particle
set explicitly exercises minimum-image body moments and torque.

The physical grid has 60 periodic x cells and 28 nonperiodic y cells. Thus
DY=1.5*DX. Two particles per x cell and three per y cell give equal particle
spacing and square particle domains. All particles use SinglePointBSpline.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

# Physical domain and requested grid/particle layout.
XMIN = -0.5
XMAX = 0.5
YMIN = -0.35
YMAX = 0.35
PERIOD_X = XMAX - XMIN
GRID_NX = 60
GRID_NY = 28
GRID_DX = PERIOD_X / GRID_NX
GRID_DY = (YMAX - YMIN) / GRID_NY
PPCX = 2
PPCY = 3
PARTICLE_DX = GRID_DX / PPCX
PARTICLE_DY = GRID_DY / PPCY
THICKNESS = GRID_DX
Z = 0.0
VOLUME = PARTICLE_DX * PARTICLE_DY * THICKNESS
SURFACE_BAND = 0.75 * PARTICLE_DX

assert math.isclose(GRID_DY, 1.5 * GRID_DX, rel_tol=0.0, abs_tol=1.0e-15)
assert math.isclose(PARTICLE_DX, PARTICLE_DY, rel_tol=0.0, abs_tol=1.0e-15)

OUT = Path(__file__).resolve().parent
PARTICLE_FILE = OUT / "mpmParticleFile_rigidBodyJamming2D"
MANIFEST_FILE = OUT / "rigidBodyJamming2D_manifest.json"

SURFACE_FLAG_INTERIOR = 0
SURFACE_FLAG_SURFACE = 2
PARTICLE_TYPE_SINGLE_POINT_BSPLINE = 1
CONTACT_GROUP = 0
MATERIAL_COPPER = 0

OBJECTS = [
    {"kind": "disk", "color": 0, "center": (-0.28, -0.20), "radius": 0.105},
    {"kind": "block", "color": 1, "center": (0.02, -0.19), "width": 0.18, "height": 0.115},
    {
        "kind": "disk",
        "color": 2,
        "center": (0.47, -0.18),
        "radius": 0.075,
        "periodicX": True,
    },
    {"kind": "block", "color": 3, "center": (-0.28, 0.14), "width": 0.16, "height": 0.155},
    {"kind": "disk", "color": 4, "center": (0.03, 0.12), "radius": 0.115},
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


def minimum_image_delta(delta: float) -> float:
    return delta - PERIOD_X * round(delta / PERIOD_X)


def object_dx(obj: dict, x: float) -> float:
    dx = x - obj["center"][0]
    return minimum_image_delta(dx) if obj.get("periodicX", False) else dx


def disk_contains(obj: dict, x: float, y: float) -> bool:
    _, cy = obj["center"]
    dx = object_dx(obj, x)
    dy = y - cy
    return dx * dx + dy * dy <= obj["radius"] ** 2 + 1.0e-14


def block_contains(obj: dict, x: float, y: float) -> bool:
    cx, cy = obj["center"]
    return abs(x - cx) <= 0.5 * obj["width"] + 1.0e-14 and abs(y - cy) <= 0.5 * obj["height"] + 1.0e-14


def disk_surface(obj: dict, x: float, y: float):
    _, cy = obj["center"]
    radius = obj["radius"]
    rx = object_dx(obj, x)
    ry = y - cy
    radial_distance = math.hypot(rx, ry)
    if radial_distance < 1.0e-14:
        nx, ny = 1.0, 0.0
        distance_to_surface = radius
    else:
        nx, ny = rx / radial_distance, ry / radial_distance
        distance_to_surface = radius - radial_distance
    sx = distance_to_surface * nx
    sy = distance_to_surface * ny
    return distance_to_surface <= SURFACE_BAND, nx, ny, sx, sy


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
    distance, nx, ny, sx, sy = min(candidates, key=lambda item: item[0])
    return distance <= SURFACE_BAND, nx, ny, sx, sy


def object_contains(obj: dict, x: float, y: float) -> bool:
    return disk_contains(obj, x, y) if obj["kind"] == "disk" else block_contains(obj, x, y)


def object_surface(obj: dict, x: float, y: float):
    return disk_surface(obj, x, y) if obj["kind"] == "disk" else block_surface(obj, x, y)


def representative_radius(obj: dict) -> float:
    if obj["kind"] == "disk":
        return obj["radius"]
    return 0.5 * math.hypot(obj["width"], obj["height"])


def conservative_center_distance(a: dict, b: dict) -> float:
    dx = a["center"][0] - b["center"][0]
    # The domain is periodic in x for every body, whether or not an individual
    # body crosses the seam at t=0.
    dx = minimum_image_delta(dx)
    dy = a["center"][1] - b["center"][1]
    return math.hypot(dx, dy)


def check_non_overlap() -> list[dict]:
    diagnostics = []
    for index, first in enumerate(OBJECTS):
        for second in OBJECTS[index + 1:]:
            center_distance = conservative_center_distance(first, second)
            margin = center_distance - representative_radius(first) - representative_radius(second)
            diagnostics.append({
                "colorA": first["color"],
                "colorB": second["color"],
                "periodicConservativeCircularMargin": margin,
            })
            if margin <= 2.0 * PARTICLE_DX:
                raise RuntimeError(
                    f"Objects {first['color']} and {second['color']} have insufficient periodic margin: {margin:g}"
                )
    return diagnostics


def candidate_indices(obj: dict):
    if obj.get("periodicX", False):
        x_indices = range(GRID_NX * PPCX)
    else:
        cx, _ = obj["center"]
        half_width = obj["radius"] if obj["kind"] == "disk" else 0.5 * obj["width"]
        ix0 = max(0, math.floor((cx - half_width - XMIN) / PARTICLE_DX) - 1)
        ix1 = min(GRID_NX * PPCX - 1, math.ceil((cx + half_width - XMIN) / PARTICLE_DX) + 1)
        x_indices = range(ix0, ix1 + 1)

    _, cy = obj["center"]
    half_height = obj["radius"] if obj["kind"] == "disk" else 0.5 * obj["height"]
    iy0 = max(0, math.floor((cy - half_height - YMIN) / PARTICLE_DY) - 1)
    iy1 = min(GRID_NY * PPCY - 1, math.ceil((cy + half_height - YMIN) / PARTICLE_DY) + 1)
    return x_indices, range(iy0, iy1 + 1)


def generate_particles():
    rows = []
    particle_id = 1
    counts = {}
    surface_counts = {}
    x_half = 0.5 * PARTICLE_DX
    y_half = 0.5 * PARTICLE_DY
    z_half = 0.5 * THICKNESS

    for obj in OBJECTS:
        x_indices, y_indices = candidate_indices(obj)
        counts[obj["color"]] = 0
        surface_counts[obj["color"]] = 0
        for ix in x_indices:
            x = XMIN + (ix + 0.5) * PARTICLE_DX
            for iy in y_indices:
                y = YMIN + (iy + 0.5) * PARTICLE_DY
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
                    MATERIAL_COPPER, PARTICLE_TYPE_SINGLE_POINT_BSPLINE, CONTACT_GROUP, obj["color"],
                    surface_flag, 0.0, 0.0, 300.0, 0.0,
                    1.0, 0,
                    x_half, 0.0, 0.0,
                    0.0, y_half, 0.0,
                    0.0, 0.0, z_half,
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
    with PARTICLE_FILE.open("w", encoding="utf-8") as particle_file:
        particle_file.write("\t".join(HEADER) + "\n")
        for row in rows:
            particle_file.write(
                " ".join(f"{value:.17g}" if isinstance(value, float) else str(value) for value in row) + "\n"
            )

    periodic_color_x = [row[1] for row in rows if int(row[10]) == 2]
    manifest = {
        "particleFile": PARTICLE_FILE.name,
        "physicalDomain": {"xmin": XMIN, "xmax": XMAX, "ymin": YMIN, "ymax": YMAX},
        "mpiPartitions": [2, 3, 1],
        "periodic": [1, 0, 0],
        "gridCells": {"xPhysical": GRID_NX, "yPhysical": GRID_NY, "yTotalWithBuffers": GRID_NY + 2},
        "gridSpacing": {"dx": GRID_DX, "dy": GRID_DY, "dyOverDx": GRID_DY / GRID_DX},
        "particlesPerCell": {"x": PPCX, "y": PPCY, "z": 1},
        "particleSpacing": {"dx": PARTICLE_DX, "dy": PARTICLE_DY},
        "particleType": PARTICLE_TYPE_SINGLE_POINT_BSPLINE,
        "thickness": THICKNESS,
        "particleVolume": VOLUME,
        "objects": OBJECTS,
        "particleCountsByColor": counts,
        "surfaceParticleCountsByColor": surface_counts,
        "periodicColor2XExtent": [min(periodic_color_x), max(periodic_color_x)],
        "nonOverlapDiagnostics": margins,
        "notes": [
            "SurfaceFlag=2 marks analytic surface particles; interior particles carry zero normal/position.",
            "All colors use ContactGroup=0; maxGridFields=2 is independent of the six global colors.",
            "ParticleColor=2 crosses the x-periodic seam and is solved as one unwrapped rigid body.",
            "Copper yield stress is 70 MPa; continuum stress control targets Syy=-140 MPa.",
        ],
    }
    MANIFEST_FILE.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"wrote {len(rows)} particles to {PARTICLE_FILE}")
    print(f"wrote manifest to {MANIFEST_FILE}")


if __name__ == "__main__":
    main()
