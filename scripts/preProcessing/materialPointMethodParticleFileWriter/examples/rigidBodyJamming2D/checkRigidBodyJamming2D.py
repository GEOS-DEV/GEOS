#!/usr/bin/env python3
"""Lightweight checks for the generated periodic rigidBodyJamming2D PFW particle file."""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET

EXPECTED_COLORS = set(range(6))
PERIODIC_COLOR = 2


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("particle_file", nargs="?", default="mpmParticleFile_rigidBodyJamming2D")
    return parser.parse_args()


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




def check_generated_xml(path: Path) -> None:
    if not path.is_file():
        return
    root = ET.parse(path).getroot()
    mesh = root.find("./Mesh/InternalMesh")
    solver = root.find("./Solvers/SolidMechanics_MPM")
    rigid = root.find("./Solvers/SolidMechanics_MPM/MPMEvents/RigidBodyMPM")
    continuum = root.find("./Solvers/SolidMechanics_MPM/MPMEvents/DeformationUpdate[@name='continuumCompression']")
    if mesh is None or solver is None or rigid is None or continuum is None:
        raise SystemExit("generated XML is missing the periodic grid, solver, or staged events")
    if mesh.get("periodic", "").replace(" ", "") != "{1,0,0}":
        raise SystemExit("generated XML should be periodic in x only")
    if (mesh.get("nx"), mesh.get("ny"), mesh.get("nz")) != ("{ 60 }", "{ 30 }", "{ 3 }"):
        raise SystemExit("generated XML should use nx=60, ny=30, nz=3")
    if solver.get("planeStrain") != "1" or solver.get("damageFieldPartitioning") != "1":
        raise SystemExit("generated XML should use plane strain and continuum DFG auto-contact")
    if solver.get("rigidBodyMaxGridFields") != "2" or rigid.get("maxGridFields") != "2":
        raise SystemExit("generated XML should reserve exactly two rigid/DFG velocity fields")
    if continuum.get("stressControl", "").replace(" ", "") != "{0,1,0}":
        raise SystemExit("generated XML should compact the continuum in y only")

def main() -> int:
    args = parse_args()
    path = Path(args.particle_file)
    if not path.is_file():
        raise SystemExit(f"missing particle file: {path}")

    names, rows = rows_from_whitespace_table(path)
    required_columns = (
        "ParticleType",
        "ParticleColor",
        "ContactGroup",
        "SurfaceFlag",
        "SurfaceNormalX",
        "SurfacePositionX",
        "RVectorXX",
        "RVectorYY",
    )
    for required in required_columns:
        if required not in names:
            raise SystemExit(f"missing required particle column {required}")

    colors = set()
    counts = {color: 0 for color in EXPECTED_COLORS}
    surface_counts = {color: 0 for color in EXPECTED_COLORS}
    periodic_x_values = []
    bad_normals = 0
    bad_interior_surface_data = 0
    bad_particle_type = 0
    nonsquare_particle_domains = 0

    for row in rows:
        color = int(float(row["ParticleColor"]))
        colors.add(color)
        if color in counts:
            counts[color] += 1
        if color == PERIODIC_COLOR:
            periodic_x_values.append(as_float(row, "PositionX"))

        if int(float(row["ParticleType"])) != 1:
            bad_particle_type += 1
        if int(float(row["ContactGroup"])) != 0:
            raise SystemExit("all colors should share the catch-all ContactGroup=0")

        half_x = abs(as_float(row, "RVectorXX"))
        half_y = abs(as_float(row, "RVectorYY"))
        if not math.isclose(half_x, half_y, rel_tol=1.0e-12, abs_tol=1.0e-14):
            nonsquare_particle_domains += 1

        flag = int(float(row["SurfaceFlag"]))
        normal = [as_float(row, f"SurfaceNormal{axis}") for axis in "XYZ"]
        surface_position = [as_float(row, f"SurfacePosition{axis}") for axis in "XYZ"]
        normal_norm = math.sqrt(sum(value * value for value in normal))
        if flag != 0:
            surface_counts[color] = surface_counts.get(color, 0) + 1
            if abs(normal_norm - 1.0) > 1.0e-8:
                bad_normals += 1
        elif normal_norm > 1.0e-14 or math.sqrt(sum(value * value for value in surface_position)) > 1.0e-14:
            bad_interior_surface_data += 1

    if colors != EXPECTED_COLORS:
        raise SystemExit(f"expected colors {sorted(EXPECTED_COLORS)}, observed {sorted(colors)}")
    if any(counts[color] <= 0 for color in EXPECTED_COLORS):
        raise SystemExit(f"empty object color found: {counts}")
    if any(surface_counts.get(color, 0) <= 0 for color in EXPECTED_COLORS):
        raise SystemExit(f"missing surface particles for at least one color: {surface_counts}")
    if bad_particle_type:
        raise SystemExit(f"expected all particles to use SinglePointBSpline ParticleType=1; bad rows: {bad_particle_type}")
    if nonsquare_particle_domains:
        raise SystemExit(f"particle r-vector x/y half-widths are not square in {nonsquare_particle_domains} rows")
    if bad_normals:
        raise SystemExit(f"surface particles with non-unit normals: {bad_normals}")
    if bad_interior_surface_data:
        raise SystemExit(f"interior particles with nonzero explicit surface data: {bad_interior_surface_data}")

    if not periodic_x_values or min(periodic_x_values) >= -0.45 or max(periodic_x_values) <= 0.45:
        raise SystemExit(
            "ParticleColor=2 should straddle the x-periodic seam with particles on both x<-0.45 and x>0.45"
        )

    check_generated_xml(path.with_name("mpm_rigidBodyJamming2D.xml"))

    print("rigidBodyJamming2D PFW particle-file/XML checks passed")
    print("particle counts by color:", counts)
    print("surface counts by color:", surface_counts)
    print(
        "periodic color-2 x extent:",
        f"[{min(periodic_x_values):.6g}, {max(periodic_x_values):.6g}]",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
