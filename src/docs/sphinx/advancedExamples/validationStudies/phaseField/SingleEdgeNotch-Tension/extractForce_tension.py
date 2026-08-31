#!/usr/bin/env python3
"""Extract top-boundary force vs. time from GEOS VTK output."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


HEX_FACES = np.array(
    [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ],
    dtype=np.int64,
)

TIME_TO_DISPLACEMENT_MM = 1.0e-5
# The model is 0.01 mm thick in z; scale the integrated force to the
# 1 mm-thick benchmark load before converting N to kN.
MODEL_THICKNESS_MM = 0.01
REFERENCE_THICKNESS_MM = 1.0
FORCE_TO_LOAD_KN = REFERENCE_THICKNESS_MM / MODEL_THICKNESS_MM * 1.0e-3

# Digitized from the right-hand reference plot for L = 0.0075 mm, eta = 0.
# Columns are displacement u [mm] and load F [kN].
REFERENCE_L_0P0075_ETA_0 = np.array(
    [
        [0.00000, 0.000],
        [0.00070, 0.090],
        [0.00140, 0.183],
        [0.00210, 0.276],
        [0.00280, 0.370],
        [0.00350, 0.463],
        [0.00420, 0.556],
        [0.00490, 0.650],
        [0.00555, 0.737],
        [0.00562, 0.704],
        [0.00569, 0.633],
        [0.00576, 0.532],
        [0.00582, 0.405],
        [0.00586, 0.282],
        [0.00589, 0.158],
        [0.00590, 0.050],
        [0.00590, 0.000],
    ],
    dtype=float,
)


def find_default_pvds(name: str) -> list[Path]:
    """Look for the GEOS output in the current directory, then in the ATS output tree."""
    for pattern in (name, f"results/{name}", f"results/run_*/{name}", f"*/{name}"):
        matches = sorted(Path.cwd().glob(pattern))
        if matches:
            return matches
    raise RuntimeError(f"No {name} found under {Path.cwd()}; pass the path explicitly")


def parse_pvd(pvd_path: Path) -> list[tuple[float, Path]]:
    root = ET.parse(pvd_path).getroot()
    datasets = []
    for dataset in root.iter("DataSet"):
        time = float(dataset.attrib["timestep"])
        vtm_path = pvd_path.parent / dataset.attrib["file"]
        datasets.append((time, vtm_path))
    return datasets


def region_piece_paths(vtm_path: Path, region_name: str) -> list[Path]:
    root = ET.parse(vtm_path).getroot()
    pieces = []
    marker = f"/{region_name}/"
    for dataset in root.iter("DataSet"):
        rel = dataset.attrib.get("file", "")
        if marker in f"/{rel}":
            pieces.append(vtm_path.parent / rel)
    return sorted(pieces)


def read_grid(vtu_path: Path) -> vtk.vtkUnstructuredGrid:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    return reader.GetOutput()


def global_y_bounds(piece_paths: list[Path]) -> tuple[float, float]:
    ymin = math.inf
    ymax = -math.inf
    for path in piece_paths:
        grid = read_grid(path)
        points = vtk_to_numpy(grid.GetPoints().GetData())
        ymin = min(ymin, float(points[:, 1].min()))
        ymax = max(ymax, float(points[:, 1].max()))
    return ymin, ymax


def quad_areas(points: np.ndarray) -> np.ndarray:
    area_1 = 0.5 * np.linalg.norm(np.cross(points[:, 1] - points[:, 0], points[:, 2] - points[:, 0]), axis=1)
    area_2 = 0.5 * np.linalg.norm(np.cross(points[:, 2] - points[:, 0], points[:, 3] - points[:, 0]), axis=1)
    return area_1 + area_2


def cell_connectivity(grid: vtk.vtkUnstructuredGrid) -> np.ndarray:
    cell_types = vtk.vtkCellTypes()
    grid.GetCellTypes(cell_types)
    unique_types = [cell_types.GetCellType(index) for index in range(cell_types.GetNumberOfTypes())]
    if unique_types != [vtk.VTK_HEXAHEDRON]:
        found = ", ".join(str(int(value)) for value in sorted(unique_types))
        raise RuntimeError(f"Only VTK_HEXAHEDRON cells are supported; found cell types {found}")

    connectivity = vtk_to_numpy(grid.GetCells().GetConnectivityArray())
    if connectivity.size != grid.GetNumberOfCells() * 8:
        raise RuntimeError("Expected 8 connectivity entries per hexahedron cell")
    return connectivity.reshape((-1, 8))


def top_boundary_force(piece_paths: list[Path], top_y: float, tol: float) -> tuple[float, float, float, float, int]:
    force_x = 0.0
    force_y = 0.0
    force_z = 0.0
    area_total = 0.0
    face_count = 0
    seen_global_cells: set[int] = set()

    for path in piece_paths:
        grid = read_grid(path)
        points = vtk_to_numpy(grid.GetPoints().GetData())
        stress_array = grid.GetCellData().GetArray("averageStress")
        cell_gid_array = grid.GetCellData().GetArray("localToGlobalMap")
        if stress_array is None:
            raise RuntimeError(f"{path} does not contain CellData array 'averageStress'")
        if cell_gid_array is None:
            raise RuntimeError(f"{path} does not contain CellData array 'localToGlobalMap'")

        stresses = vtk_to_numpy(stress_array)
        gids = vtk_to_numpy(cell_gid_array)
        connectivity = cell_connectivity(grid)

        if seen_global_cells:
            keep = np.array([int(gid) not in seen_global_cells for gid in gids], dtype=bool)
            connectivity = connectivity[keep]
            stresses = stresses[keep]
            gids = gids[keep]
        seen_global_cells.update(int(gid) for gid in gids)

        piece_force_x = 0.0
        piece_force_y = 0.0
        piece_force_z = 0.0
        piece_area = 0.0
        piece_faces = 0
        for face in HEX_FACES:
            face_point_ids = connectivity[:, face]
            face_points = points[face_point_ids]
            mask = np.all(np.abs(face_points[:, :, 1] - top_y) <= tol, axis=1)
            if not np.any(mask):
                continue

            areas = quad_areas(face_points[mask])
            face_stresses = stresses[mask]
            piece_force_x += float(np.dot(face_stresses[:, 5], areas))
            piece_force_y += float(np.dot(face_stresses[:, 1], areas))
            piece_force_z += float(np.dot(face_stresses[:, 3], areas))
            piece_area += float(np.sum(areas))
            piece_faces += int(np.count_nonzero(mask))

        force_x += piece_force_x
        force_y += piece_force_y
        force_z += piece_force_z
        area_total += piece_area
        face_count += piece_faces

    return force_x, force_y, force_z, area_total, face_count


def write_csv(rows: list[dict[str, float | int]], csv_path: Path) -> None:
    fieldnames = ["time", "force_x", "force_y", "force_z", "top_area", "top_faces"]
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_plot(rows: list[dict[str, float | int]], png_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if shutil.which("latex") is None:
        # No TeX installation available: fall back to matplotlib's own mathtext.
        plt.rcParams.update({
            "text.usetex": False,
            "font.family": "serif",
        })
    else:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",  # Use a serif font family (required for LaTeX rendering)
            "font.serif": ["Palatino"],  # Specify the custom LaTeX font
            "mathtext.fontset": "custom",  # Use the custom LaTeX font for math text
            "mathtext.default": "regular",  # Set the math font to regular
        })

    displacement_mm = [float(row["time"]) * TIME_TO_DISPLACEMENT_MM for row in rows]
    load_kn = [float(row["force_y"]) * FORCE_TO_LOAD_KN for row in rows]

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.plot(
        REFERENCE_L_0P0075_ETA_0[:, 0],
        REFERENCE_L_0P0075_ETA_0[:, 1],
        color="black",
        linewidth=1.8,
        linestyle="-",
        label=r"Reference",
    )
    ax.plot(displacement_mm,
            load_kn,
            color="#cd1313",
            linestyle="none",
            marker="o",
            markersize=5,
            markeredgewidth=1.8,
            markerfacecolor="none",
            label="GEOS")
    ax.set_xlabel("Displacement $u$ [mm]", fontsize=12)
    ax.set_ylabel("Applied force [kN]", fontsize=12)
    ax.grid(True, color="#d9d9d9", linewidth=0.8)
    ax.legend(frameon=True, facecolor='white', edgecolor='black', fontsize=12, loc="upper left")
    # ax.set_title("Single edge notch tension: top boundary force")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pvd",
        nargs="?",
        type=Path,
        help="Path to the GEOS PVD time series. Defaults to singleEdgeNotchTension.pvd in the "
        "current directory, or results/run_*/singleEdgeNotchTension.pvd.",
    )
    parser.add_argument("--region", default="Region1", help="Cell element region to integrate.")
    parser.add_argument("--csv", type=Path, help="Output CSV path.")
    parser.add_argument("--png", type=Path, help="Output PNG path.")
    parser.add_argument(
        "--top-y",
        type=float,
        help="Top boundary y-coordinate. Defaults to the maximum y-coordinate in the first frame.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1.0e-9,
        help="Absolute tolerance for detecting faces on the top boundary.",
    )
    args = parser.parse_args()

    pvd_path = args.pvd.resolve() if args.pvd else find_default_pvds("singleEdgeNotchTension.pvd")[0]
    out_dir = pvd_path.parent
    csv_path = args.csv.resolve() if args.csv else out_dir / "top_boundary_force.csv"
    png_path = args.png.resolve() if args.png else out_dir / "top_boundary_force.png"

    datasets = parse_pvd(pvd_path)
    if not datasets:
        raise RuntimeError(f"No timesteps found in {pvd_path}")

    first_pieces = region_piece_paths(datasets[0][1], args.region)
    if not first_pieces:
        raise RuntimeError(f"No VTU pieces for region '{args.region}' in {datasets[0][1]}")
    ymin, ymax = global_y_bounds(first_pieces)
    top_y = ymax if args.top_y is None else args.top_y
    print(f"Using top y = {top_y:g} (mesh y range {ymin:g} to {ymax:g})", flush=True)

    rows: list[dict[str, float | int]] = []
    for index, (time, vtm_path) in enumerate(datasets, start=1):
        pieces = region_piece_paths(vtm_path, args.region)
        force_x, force_y, force_z, area, faces = top_boundary_force(pieces, top_y, args.tol)
        rows.append(
            {
                "time": time,
                "force_x": force_x,
                "force_y": force_y,
                "force_z": force_z,
                "top_area": area,
                "top_faces": faces,
            }
        )
        print(
            f"[{index:03d}/{len(datasets):03d}] time={time:g} "
            f"Fy={force_y:.12g} area={area:.12g} faces={faces}",
            flush=True,
        )

    write_csv(rows, csv_path)
    write_plot(rows, png_path)
    print(f"Wrote {csv_path}", flush=True)
    print(f"Wrote {png_path}", flush=True)


if __name__ == "__main__":
    main()
