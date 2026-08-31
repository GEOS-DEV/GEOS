#!/usr/bin/env python3
"""Extract combined top-boundary force vs. time from restarted GEOS VTK output."""

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

COMPONENT_INDEX = {
    "x": "force_x",
    "y": "force_y",
    "z": "force_z",
}

TIME_TO_DISPLACEMENT_MM = 1.0e-5
# The model is 0.01 mm thick in z; scale the integrated force to the
# 1 mm-thick benchmark load before converting N to kN.
MODEL_THICKNESS_MM = 0.01
REFERENCE_THICKNESS_MM = 1.0
FORCE_TO_LOAD_KN = REFERENCE_THICKNESS_MM / MODEL_THICKNESS_MM * 1.0e-3

# Digitized from the right-hand reference plot for the shear benchmark,
# eta = 0. Columns are displacement u [mm] and load F [kN].
REFERENCE_SHEAR_ETA_0 = np.array(
    [
        [0.00000, 0.000],
        [0.00090, 0.059],
        [0.00180, 0.118],
        [0.00270, 0.177],
        [0.00360, 0.236],
        [0.00450, 0.295],
        [0.00540, 0.354],
        [0.00630, 0.414],
        [0.00720, 0.474],
        [0.00810, 0.539],
        [0.00870, 0.558],
        [0.00950, 0.572],
        [0.01030, 0.577],
        [0.01110, 0.571],
        [0.01180, 0.548],
        [0.01250, 0.508],
        [0.01320, 0.458],
        [0.01390, 0.412],
        [0.01460, 0.382],
        [0.01480, 0.372],
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


def parse_pvd(pvd_path: Path) -> list[tuple[float, Path, Path]]:
    root = ET.parse(pvd_path).getroot()
    datasets = []
    for dataset in root.iter("DataSet"):
        time = float(dataset.attrib["timestep"])
        vtm_path = pvd_path.parent / dataset.attrib["file"]
        datasets.append((time, vtm_path, pvd_path))
    return datasets


def parse_pvds(pvd_paths: list[Path]) -> list[tuple[float, Path, Path]]:
    datasets = []
    for pvd_path in pvd_paths:
        datasets.extend(parse_pvd(pvd_path))

    datasets.sort(key=lambda item: (item[0], str(item[1])))
    deduped = []
    seen_times: set[float] = set()
    for dataset in datasets:
        time = dataset[0]
        if time in seen_times:
            continue
        deduped.append(dataset)
        seen_times.add(time)
    return deduped


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

        for face in HEX_FACES:
            face_point_ids = connectivity[:, face]
            face_points = points[face_point_ids]
            mask = np.all(np.abs(face_points[:, :, 1] - top_y) <= tol, axis=1)
            if not np.any(mask):
                continue

            areas = quad_areas(face_points[mask])
            face_stresses = stresses[mask]
            force_x += float(np.dot(face_stresses[:, 5], areas))
            force_y += float(np.dot(face_stresses[:, 1], areas))
            force_z += float(np.dot(face_stresses[:, 3], areas))
            area_total += float(np.sum(areas))
            face_count += int(np.count_nonzero(mask))

    return force_x, force_y, force_z, area_total, face_count


def write_csv(rows: list[dict[str, float | int | str]], csv_path: Path) -> None:
    fieldnames = ["time", "force_x", "force_y", "force_z", "top_area", "top_faces", "source_pvd"]
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_plot(rows: list[dict[str, float | int | str]], png_path: Path, component: str) -> None:
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

    force_key = COMPONENT_INDEX[component]
    displacement_mm = [float(row["time"]) * TIME_TO_DISPLACEMENT_MM for row in rows]
    load_kn = [float(row[force_key]) * FORCE_TO_LOAD_KN for row in rows]

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    if component == "x":
        ax.plot(
            REFERENCE_SHEAR_ETA_0[:, 0],
            REFERENCE_SHEAR_ETA_0[:, 1],
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
    ax.legend(frameon=True, facecolor='white', edgecolor='black', fontsize=12)
    # ax.set_title("Single edge notch shear: top boundary force")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pvd",
        nargs="*",
        type=Path,
        help="Path(s) to GEOS PVD time series. Defaults to singleEdgeNotchShear.pvd in the "
        "current directory, or results/run_*/singleEdgeNotchShear.pvd.",
    )
    parser.add_argument("--region", default="Region1", help="Cell element region to integrate.")
    parser.add_argument("--csv", type=Path, help="Output CSV path.")
    parser.add_argument("--png", type=Path, help="Output PNG path.")
    parser.add_argument(
        "--component",
        choices=sorted(COMPONENT_INDEX),
        default="x",
        help="Force component to plot. CSV always contains all components.",
    )
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

    pvd_paths = [path.resolve() for path in args.pvd]
    if not pvd_paths:
        pvd_paths = find_default_pvds("singleEdgeNotchShear.pvd")

    out_dir = pvd_paths[0].parent
    csv_path = args.csv.resolve() if args.csv else out_dir / "top_boundary_force_combined.csv"
    png_path = args.png.resolve() if args.png else out_dir / f"top_boundary_force_{args.component}_combined.png"

    datasets = parse_pvds(pvd_paths)
    if not datasets:
        raise RuntimeError("No timesteps found in PVD input")

    first_pieces = region_piece_paths(datasets[0][1], args.region)
    if not first_pieces:
        raise RuntimeError(f"No VTU pieces for region '{args.region}' in {datasets[0][1]}")
    ymin, ymax = global_y_bounds(first_pieces)
    top_y = ymax if args.top_y is None else args.top_y
    print("PVD input:")
    for pvd_path in pvd_paths:
        print(f"  {pvd_path}")
    print(f"Using top y = {top_y:g} (mesh y range {ymin:g} to {ymax:g})", flush=True)
    print(f"Combined timesteps = {len(datasets)}", flush=True)

    rows: list[dict[str, float | int | str]] = []
    for index, (time, vtm_path, source_pvd) in enumerate(datasets, start=1):
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
                "source_pvd": source_pvd,
            }
        )
        print(
            f"[{index:03d}/{len(datasets):03d}] time={time:g} "
            f"Fx={force_x:.12g} Fy={force_y:.12g} area={area:.12g} faces={faces}",
            flush=True,
        )

    write_csv(rows, csv_path)
    write_plot(rows, png_path, args.component)
    print(f"Wrote {csv_path}", flush=True)
    print(f"Wrote {png_path}", flush=True)


if __name__ == "__main__":
    main()
