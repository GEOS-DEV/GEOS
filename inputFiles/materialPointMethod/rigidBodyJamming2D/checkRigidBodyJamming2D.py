#!/usr/bin/env python3
"""Verification checks for the MPI/periodic rigidBodyJamming2D case.

The script always checks the committed/generated particle geometry and XML. It
checks GEOS history files when present, or requires them with
``--require-run-output``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
PARTICLE_FILE = ROOT / "mpmParticleFile_rigidBodyJamming2D"
MANIFEST_FILE = ROOT / "rigidBodyJamming2D_manifest.json"
XML_FILE = ROOT / "mpm_rigidBodyJamming2D.xml"
RIGID_HISTORY = ROOT / "rigidBodyHistory.csv"
BOX_HISTORY = ROOT / "boxAverageHistory.csv"
YIELD_STRESS = 7.0e7
TARGET_STRESS_YY = -2.0 * YIELD_STRESS


def load_manifest() -> dict:
    if not MANIFEST_FILE.exists():
        raise FileNotFoundError(f"Missing {MANIFEST_FILE.name}; run generateRigidBodyJamming2DParticles.py first")
    return json.loads(MANIFEST_FILE.read_text(encoding="utf-8"))


def read_particles():
    if not PARTICLE_FILE.exists():
        raise FileNotFoundError(f"Missing {PARTICLE_FILE.name}; run generateRigidBodyJamming2DParticles.py first")
    with PARTICLE_FILE.open("r", encoding="utf-8") as particle_file:
        header = particle_file.readline().strip().split("\t")
        indices = {name: index for index, name in enumerate(header)}
        rows = []
        for line_number, line in enumerate(particle_file, start=2):
            if not line.strip():
                continue
            values = line.split()
            if len(values) != len(header):
                raise AssertionError(f"line {line_number} has {len(values)} values; expected {len(header)}")
            rows.append({name: float(values[index]) for name, index in indices.items()})
    return rows


def check_xml():
    root = ET.parse(XML_FILE).getroot()
    mesh = root.find("./Mesh/InternalMesh")
    particle_mesh = root.find("./Mesh/ParticleMesh")
    solver = root.find("./Solvers/SolidMechanics_MPM")
    rigid = root.find("./Solvers/SolidMechanics_MPM/MPMEvents/RigidBodyMPM")
    continuum = root.find("./Solvers/SolidMechanics_MPM/MPMEvents/DeformationUpdate[@name='continuumCompression']")
    if None in (mesh, particle_mesh, solver, rigid, continuum):
        raise AssertionError("missing required mesh, solver, or event node in XML")

    if mesh.get("periodic", "").replace(" ", "") != "{1,0,0}":
        raise AssertionError(f"expected periodic x only, observed {mesh.get('periodic')}")
    if mesh.get("nx") != "{60}" or mesh.get("ny") != "{30}" or mesh.get("nz") != "{3}":
        raise AssertionError("expected grid counts nx=60, ny=30, nz=3")
    if particle_mesh.get("particleTypes", "").replace(" ", "") != "{SinglePointBSpline}":
        raise AssertionError("expected SinglePointBSpline particle mesh")
    if solver.get("planeStrain") != "1":
        raise AssertionError("expected planeStrain=1")
    if solver.get("boundaryConditionTypes", "").replace(" ", "") != "{0,0,3,3,1,1}":
        raise AssertionError("expected periodic/outflow x, Contact y, and symmetry z boundary entries")
    if solver.get("rigidBodyMaxGridFields") != "2" or rigid.get("maxGridFields") != "2":
        raise AssertionError("expected a two-field rigid/contact allocation")
    if solver.get("damageFieldPartitioning") != "1":
        raise AssertionError("continuum handoff should restore two-field DFG auto-contact")
    if continuum.get("stressControl", "").replace(" ", "") != "{0,1,0}":
        raise AssertionError("continuum verification should stress-control y only")

    x0, x1 = [float(value) for value in mesh.get("xCoords").strip("{}").split(",")]
    y0, y1 = [float(value) for value in mesh.get("yCoords").strip("{}").split(",")]
    dx = (x1 - x0) / 60.0
    dy = (y1 - y0) / 30.0
    if not math.isclose(dy, 1.5 * dx, rel_tol=0.0, abs_tol=1.0e-14):
        raise AssertionError(f"expected DY=1.5*DX, observed DX={dx}, DY={dy}")

    print("xml: periodic x, 2-field rigid event, y loading, and rectangular grid checks passed")


def check_geometry(manifest: dict, rows: list[dict]):
    counts = {int(key): int(value) for key, value in manifest["particleCountsByColor"].items()}
    colors = sorted(counts)
    observed = {color: 0 for color in colors}
    surface = {color: 0 for color in colors}
    periodic_color_x = []

    for row in rows:
        color = int(row["ParticleColor"])
        if color not in observed:
            raise AssertionError(f"unexpected ParticleColor {color}")
        observed[color] += 1
        if color == 2:
            periodic_color_x.append(row["PositionX"])
        if int(row["ContactGroup"]) != 0:
            raise AssertionError("all colors should share catch-all ContactGroup=0")
        if int(row["MaterialType"]) != 0:
            raise AssertionError("verification expects material index 0/copper")
        if int(row["ParticleType"]) != 1:
            raise AssertionError("verification expects SinglePointBSpline ParticleType=1")

        half_x = abs(row["RVectorXX"])
        half_y = abs(row["RVectorYY"])
        if not math.isclose(half_x, half_y, rel_tol=1.0e-12, abs_tol=1.0e-14):
            raise AssertionError("x/y particle half-widths should be equal (square particles)")

        normal = (row["SurfaceNormalX"], row["SurfaceNormalY"], row["SurfaceNormalZ"])
        surface_position = (row["SurfacePositionX"], row["SurfacePositionY"], row["SurfacePositionZ"])
        normal_norm = math.sqrt(sum(value * value for value in normal))
        if int(row["SurfaceFlag"]) == 2:
            surface[color] += 1
            if abs(normal_norm - 1.0) > 1.0e-10:
                raise AssertionError(f"surface particle has non-unit normal: {normal_norm}")
            if sum(normal[index] * surface_position[index] for index in range(3)) < -1.0e-12:
                raise AssertionError("surface-position vector points inward relative to the surface normal")
        elif normal_norm > 1.0e-14 or math.sqrt(sum(value * value for value in surface_position)) > 1.0e-14:
            raise AssertionError("interior particle should not contribute explicit contact normal/position")

    if observed != counts:
        raise AssertionError(f"particle counts mismatch: observed {observed}, manifest {counts}")
    if any(value == 0 for value in surface.values()):
        raise AssertionError(f"each color must have surface particles; observed {surface}")

    for diagnostic in manifest["nonOverlapDiagnostics"]:
        if diagnostic["periodicConservativeCircularMargin"] <= 2.0 * manifest["particleSpacing"]["dx"]:
            raise AssertionError(f"insufficient initial periodic non-overlap margin: {diagnostic}")

    if min(periodic_color_x) >= -0.45 or max(periodic_color_x) <= 0.45:
        raise AssertionError("ParticleColor=2 should have particles on both sides of the x-periodic seam")

    grid = manifest["gridSpacing"]
    spacing = manifest["particleSpacing"]
    if not math.isclose(grid["dy"], 1.5 * grid["dx"], rel_tol=0.0, abs_tol=1.0e-14):
        raise AssertionError("manifest grid spacing does not satisfy DY=1.5*DX")
    if not math.isclose(spacing["dx"], spacing["dy"], rel_tol=0.0, abs_tol=1.0e-14):
        raise AssertionError("manifest particle spacing is not square")

    print(
        f"geometry: {len(rows)} SinglePointBSpline particles across colors {colors}; "
        "periodic non-overlap, seam crossing, and exact surface checks passed"
    )


def check_rigid_history(require: bool):
    if not RIGID_HISTORY.exists():
        if require:
            raise AssertionError(f"missing required {RIGID_HISTORY.name}")
        print(f"post-run: {RIGID_HISTORY.name} not present; skipping rigid history checks")
        return
    with RIGID_HISTORY.open("r", encoding="utf-8") as history_file:
        rows = list(csv.DictReader(history_file))
    if not rows:
        raise AssertionError("rigidBodyHistory.csv is empty")
    colors = {int(float(row["color"])) for row in rows}
    if colors != set(range(6)):
        raise AssertionError(f"rigid history should contain colors 0..5, observed {sorted(colors)}")
    max_force = 0.0
    max_ke = 0.0
    for row in rows:
        force = [float(row[name]) for name in ("forceX", "forceY", "forceZ")]
        max_force = max(max_force, math.sqrt(sum(value * value for value in force)))
        max_ke = max(max_ke, float(row["kineticEnergy"]))
    if max_force <= 0.0:
        raise AssertionError("rigid history did not record a nonzero resultant force")
    print(f"post-run: rigid history present; max body force={max_force:.6e}, max KE={max_ke:.6e}")


def check_box_history(require: bool):
    if not BOX_HISTORY.exists():
        if require:
            raise AssertionError(f"missing required {BOX_HISTORY.name}")
        print(f"post-run: {BOX_HISTORY.name} not present; skipping continuum stress checks")
        return
    with BOX_HISTORY.open("r", encoding="utf-8") as history_file:
        rows = list(csv.DictReader(history_file, skipinitialspace=True))
    if not rows:
        raise AssertionError("boxAverageHistory.csv is empty")
    syy = float(rows[-1]["Syy"])
    tolerance = 0.30 * abs(TARGET_STRESS_YY)
    if abs(syy - TARGET_STRESS_YY) > tolerance:
        raise AssertionError(
            f"final Syy is not within 30% of target {TARGET_STRESS_YY:.6e}: Syy={syy:.6e}"
        )
    print(f"post-run: final continuum y-stress check passed; Syy={syy:.6e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-run-output", action="store_true", help="fail if GEOS history files are absent")
    args = parser.parse_args()
    check_xml()
    manifest = load_manifest()
    rows = read_particles()
    check_geometry(manifest, rows)
    check_rigid_history(args.require_run_output)
    check_box_history(args.require_run_output)
    print("rigid-body MPM MPI/periodic jamming verification checks passed")


if __name__ == "__main__":
    main()
