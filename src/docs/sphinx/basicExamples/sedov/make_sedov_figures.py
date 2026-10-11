"""Plot the elastic Sedov example from GEOS histories or the checked-in data.

The documentation build uses the PNG files; it does not execute this script.
Use --history to collect new results. Use --data to reproduce the figures
without GEOS or h5py. All stress values are in Pa and all times are in s in
the data file. Plot labels disclose the changes of units.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ElementTree

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="geos-sedov-matplotlib-"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

HERE = Path(__file__).resolve().parent
INPUT = (HERE / "../../../../../inputFiles/solidMechanics/sedov_basicExample.xml").resolve()
COLORS = ["#64748b", "#c97923", "#16877a", "#245b9e"]


def configuration(path):
    """Read the physical parameters from the standalone XML input."""
    problem = ElementTree.parse(path)
    material = problem.find("Constitutive/ElasticIsotropic")
    mesh = problem.find("Mesh/InternalMesh")
    events = problem.find("Events")
    return {
        "density_kg_m3": float(material.get("defaultDensity")),
        "bulk_modulus_Pa": float(material.get("defaultBulkModulus")),
        "shear_modulus_Pa": float(material.get("defaultShearModulus")),
        "cell_count": int(mesh.get("nx").strip("{} "))**3,
        "end_time_s": float(events.get("maxTime")),
        "time_step_s": float(events.find("PeriodicEvent[@name='solverApplications']").get("forceDt")),
    }


def collect(history, parameters):
    """Keep the plotted snapshots and the complete source-node history."""
    import h5py

    with h5py.File(history, "r") as handle:
        times = handle["shale_stress Time"][...].ravel()
        if times[0] != 0 or np.any(np.diff(times) <= 0):
            raise ValueError("Expected a true initial record and increasing times.")
        for field in ["averageStress", "totalDisplacement", "velocity"]:
            if not np.array_equal(times, handle[field + " Time"][...].ravel()):
                raise ValueError("The four fields must have the same sampling times.")
        requested = [0, 25 * parameters["time_step_s"], 50 * parameters["time_step_s"], times[-1]]
        indices = np.array([int(np.argmin(np.abs(times - time))) for time in requested])
        cells = handle["shale_stress elementCenter"][0]
        nodes = handle["velocity ReferencePosition"][0]
        # Coordinate meta-collectors write when the selected set changes.
        # For this fixed mesh, only record zero contains coordinates; later
        # reserved rows can contain zeros. Do not use them as moving positions.
        if not np.array_equal(cells, handle["averageStress elementCenter"][0]):
            raise ValueError("The two cell fields must use the same initial coordinates.")
        if not np.array_equal(nodes, handle["totalDisplacement ReferencePosition"][0]):
            raise ValueError("The two nodal fields must use the same reference coordinates.")
        # PackCollection flattens the last two dimensions: 8 points x 6 components.
        stress = handle["shale_stress"][indices].reshape(len(indices), len(cells), 8, 6)
        source = np.flatnonzero(np.all(np.isclose(nodes, [1, 1, 1]), axis=1))
        if len(source) != 1:
            raise ValueError("Expected one node at (1, 1, 1) m; use the serial example.")
        # The regular cubes have equal quadrature weights.
        data = {
            "times": times[indices],
            "cell_centers": cells,
            "node_positions": nodes,
            "stress": stress.mean(axis=2),
            "average_stress": handle["averageStress"][indices],
            "displacement": handle["totalDisplacement"][indices],
            "velocity": handle["velocity"][indices],
            "history_times": times,
            "source_velocity": handle["velocity"][:, source[0], :],
        }
    return data


def validate(data, parameters):
    """Check the physical initial conditions and the data used by each plot."""
    cells = data["cell_centers"]
    nodes = data["node_positions"]
    if len(cells) != parameters["cell_count"] or len(nodes) != 1331:
        raise ValueError("Expected the unmodified 10 x 10 x 10, serial C3D8 example.")
    if data["times"][0] != 0 or data["times"][-1] >= parameters["end_time_s"]:
        raise ValueError("History snapshots must start at zero and precede the final cleanup output.")
    expected_times = np.arange(round(parameters["end_time_s"] / parameters["time_step_s"])) * parameters["time_step_s"]
    if data["history_times"].shape != expected_times.shape or not np.allclose(
            data["history_times"], expected_times, rtol=0, atol=1e-15):
        raise ValueError("Expected the complete history: 100 records from 0 to 0.99 ms at 10 microsecond intervals.")
    for name, values in data.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Non-finite values in {name}.")
    source = np.all(np.isclose(cells, [.5, .5, .5]), axis=1)
    expected = np.zeros((len(cells), 6))
    expected[source, :3] = -1e6
    if source.sum() != 1 or not np.allclose(data["stress"][0], expected, rtol=0, atol=1e-7):
        raise ValueError("The first record does not contain the prescribed stress pulse.")
    if not np.allclose(data["average_stress"][0], expected, rtol=0, atol=1e-7):
        raise ValueError("The initial averageStress does not match the prescribed stress.")
    if np.any(data["displacement"][0]) or np.any(data["velocity"][0]):
        raise ValueError("The first record must precede the first solver step.")
    for component in range(3):
        on_plane = np.isclose(nodes[:, component], 0)
        if np.any(data["velocity"][:, on_plane, component]):
            raise ValueError("The normal velocity must be zero on each symmetry plane.")
    if not np.allclose(data["source_velocity"], data["source_velocity"][:, :1], rtol=0, atol=1e-13):
        raise ValueError("The three source-node velocity components must obey the cube symmetry.")


def cube_faces(low, high):
    """Return the six faces of a cube for the setup diagram."""
    x, y, z = low
    X, Y, Z = high
    return [[(x, y, z), (X, y, z), (X, Y, z), (x, Y, z)], [(x, y, Z), (X, y, Z), (X, Y, Z), (x, Y, Z)],
            [(x, y, z), (X, y, z), (X, y, Z), (x, y, Z)], [(x, Y, z), (X, Y, z), (X, Y, Z), (x, Y, Z)],
            [(x, y, z), (x, Y, z), (x, Y, Z), (x, y, Z)], [(X, y, z), (X, Y, z), (X, Y, Z), (X, y, Z)]]


def save(figure, output, name):
    figure.savefig(output / name, dpi=200, facecolor="white")
    plt.close(figure)


def setup_figure(data, output):
    figure = plt.figure(figsize=(10, 5.2), layout="constrained")
    axes = figure.add_subplot(121, projection="3d")
    axes.add_collection3d(
        Poly3DCollection(cube_faces((0, 0, 0), (10, 10, 10)),
                         facecolor="#dfe9f3",
                         edgecolor="#7a91a8",
                         alpha=.12,
                         linewidth=.8))
    for value in range(1, 10):
        axes.plot([0, 10], [value, value], [0, 0], color="#b5c5d5", lw=.5)
        axes.plot([value, value], [0, 10], [0, 0], color="#b5c5d5", lw=.5)
    axes.add_collection3d(
        Poly3DCollection(cube_faces((0, 0, 0), (1, 1, 1)), facecolor="#c34e3b", edgecolor="#8c3025", linewidth=.8))
    axes.set(xlim=(0, 10),
             ylim=(0, 10),
             zlim=(0, 10),
             xticks=[0, 5, 10],
             yticks=[0, 5, 10],
             zticks=[0, 5, 10],
             xlabel="$x$ (m)",
             ylabel="$y$ (m)",
             zlabel="$z$ (m)",
             title="Computational domain: one octant")
    axes.set_box_aspect((1, 1, 1))
    axes.view_init(elev=23, azim=-55)
    for axis in [axes.xaxis, axes.yaxis, axes.zaxis]:
        axis.pane.fill = False
    zoom = figure.add_subplot(122)
    zoom.set(xlim=(0, 4),
             ylim=(0, 4),
             xticks=range(5),
             yticks=range(5),
             aspect="equal",
             xlabel="$x$ (m)",
             ylabel="$y$ (m)",
             title="Source-cell detail at $z=0.5$ m")
    zoom.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor="#c34e3b", edgecolor="#8c3025"))
    zoom.grid(color="#cbd5e1", linewidth=.65)
    zoom.plot(.5, .5, "o", color="white", ms=5)
    zoom.annotate("One loaded cell\n$\\sigma_{xx}=\\sigma_{yy}=\\sigma_{zz}=-1$ MPa",
                  xy=(.5, .5),
                  xytext=(1.3, 2.3),
                  fontsize=10,
                  arrowprops={
                      "arrowstyle": "->",
                      "color": "#475569"
                  })
    zoom.text(.5,
              -.24,
              "Normal velocity is zero on $x=0$, $y=0$, and $z=0$.",
              transform=zoom.transAxes,
              ha="center",
              fontsize=10)
    save(figure, output, "SedovSetup.png")


def stress_maps(data, output):
    cells = data["cell_centers"]
    selected = np.flatnonzero(np.isclose(cells[:, 2], .5))
    selected = selected[np.lexsort((cells[selected, 0], cells[selected, 1]))]
    figure, axes = plt.subplots(2, 2, figsize=(9.5, 7.3), layout="constrained", sharex=True, sharey=True)
    normalizer = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for index, axis in enumerate(axes.flat):
        image = axis.imshow(data["stress"][index, selected, 0].reshape(10, 10) / 1e6,
                            origin="lower",
                            extent=(0, 10, 0, 10),
                            cmap="RdBu_r",
                            norm=normalizer,
                            interpolation="nearest")
        axis.set(title=f"$t={data['times'][index] * 1e3:.2f}$ ms",
                 xlim=(0, 4),
                 ylim=(0, 4),
                 xticks=range(5),
                 yticks=range(5))
        axis.grid(which="major", color="#d3d9df", lw=.35)
    for axis in axes[-1]:
        axis.set_xlabel("$x$ (m)")
    for axis in axes[:, 0]:
        axis.set_ylabel("$y$ (m)")
    figure.suptitle("Material stress $\\sigma_{xx}$ near the source, $z=0.5$ m", fontsize=13)
    figure.colorbar(image, ax=axes, label="$\\sigma_{xx}$ (MPa); negative values mean compression", shrink=.85)
    save(figure, output, "SedovStressMaps.png")


def profiles(data, output):
    cells = data["cell_centers"]
    nodes = data["node_positions"]
    cell_line = np.flatnonzero(np.all(np.isclose(cells[:, 1:], .5), axis=1))
    node_line = np.flatnonzero(np.all(np.isclose(nodes[:, 1:], 0), axis=1))
    cell_line = cell_line[np.argsort(cells[cell_line, 0])]
    node_line = node_line[np.argsort(nodes[node_line, 0])]
    if len(cell_line) != 10 or len(node_line) != 11:
        raise ValueError("Expected all 10 cells and 11 nodes on the selected lines.")
    figure, axes = plt.subplots(3, 1, figsize=(9.2, 8), layout="constrained", sharex=True)
    for index, color in enumerate(COLORS):
        label = f"{data['times'][index] * 1e3:.2f} ms"
        axes[0].plot(cells[cell_line, 0],
                     data["stress"][index, cell_line, 0] / 1e6,
                     color=color,
                     marker="o",
                     ms=4,
                     label=label)
        axes[1].plot(nodes[node_line, 0],
                     data["displacement"][index, node_line, 0] * 1e6,
                     color=color,
                     marker="o",
                     ms=4)
        axes[2].plot(nodes[node_line, 0], data["velocity"][index, node_line, 0], color=color, marker="o", ms=4)
    axes[0].set(ylabel="$\\sigma_{xx}$ (MPa)", title="Stress: cell centers at $y=z=0.5$ m")
    axes[0].legend(title="Time", ncol=4, loc="lower right", frameon=True, edgecolor="white", fontsize=10)
    axes[1].set(ylabel="$u_x$ ($\\mu$m)", title="Displacement: nodes at $y=z=0$ m")
    axes[2].set(ylabel="$v_x$ (m/s)", title="Velocity: nodes at $y=z=0$ m", xlabel="$x$ (m)")
    for axis in axes:
        axis.set(xlim=(0, 10), xticks=range(0, 11, 2))
        axis.grid(color="#e2e8f0", lw=.6)
        axis.spines[["top", "right"]].set_visible(False)
    save(figure, output, "SedovProfiles.png")


def velocity_history(data, output):
    figure, axis = plt.subplots(figsize=(9.2, 4.5), layout="constrained")
    time = data["history_times"] * 1e3
    velocity = data["source_velocity"]
    axis.plot(time, velocity[:, 0], color=COLORS[-1], lw=2, label="$v_x=v_y=v_z$")
    axis.plot(time, np.linalg.norm(velocity, axis=1), color=COLORS[2], lw=2, ls="--", label="Speed $|\\mathbf{v}|$")
    axis.set(xlabel="Time (ms)",
             ylabel="Velocity / speed (m/s)",
             xlim=(0, 1),
             title="Motion of the source node at $(1,1,1)$ m")
    axis.grid(color="#e2e8f0", lw=.6)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    save(figure, output, "SedovVelocityHistory.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--history", type=Path, help="GEOS sedovHistory.hdf5 from a serial run")
    inputs.add_argument("--data", type=Path, help="Curated sedovResults.npz; defaults to the checked-in file")
    parser.add_argument("--input", type=Path, default=INPUT, help="GEOS XML input used for the run")
    parser.add_argument("--run-log", type=Path, help="GEOS log with the executable version and compiler")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    parameters = configuration(args.input)
    if args.history:
        data = collect(args.history, parameters)
    else:
        with np.load(args.data or HERE / "sedovResults.npz", allow_pickle=False) as saved:
            data = dict(saved)
    validate(data, parameters)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 11,
        "axes.labelsize": 11,
        "legend.fontsize": 10,
        "lines.linewidth": 1.7,
        "font.family": "DejaVu Sans"
    })
    setup_figure(data, args.output_dir)
    stress_maps(data, args.output_dir)
    profiles(data, args.output_dir)
    velocity_history(data, args.output_dir)
    if args.history:
        np.savez_compressed(args.output_dir / "sedovResults.npz", **data)
        manifest = {
            "input_file":
            "inputFiles/solidMechanics/sedov_basicExample.xml",
            "input_sha256":
            hashlib.sha256(args.input.read_bytes()).hexdigest(),
            "history_sha256":
            hashlib.sha256(args.history.read_bytes()).hexdigest(),
            "curated_data_sha256":
            hashlib.sha256((args.output_dir / "sedovResults.npz").read_bytes()).hexdigest(),
            "parameters":
            parameters,
            "snapshots_s":
            data["times"].tolist(),
            "history_records":
            len(data["history_times"]),
            "source_node_m": [1, 1, 1],
            "stress_components": ["xx", "yy", "zz", "yz", "xz", "xy"],
            "stress_reduction":
            "Arithmetic mean of eight equally weighted quadrature points in each regular cube.",
            "spatial_coordinates":
            "Reference coordinates; all selected points are retained, with no smoothing.",
            "initial_average_stress_max_error_Pa":
            float(np.max(np.abs(data["average_stress"][0] - data["stress"][0]))),
            "average_stress_at_last_sample_equals_initial":
            np.array_equal(data["average_stress"][-1], data["average_stress"][0]),
            "unit_changes_for_plots":
            "s to ms, Pa to MPa, displacement m to micrometers; velocity remains m/s.",
        }
        if args.run_log:
            lines = args.run_log.read_text().splitlines()
            manifest["runtime"] = [
                line.strip() for line in lines if line.startswith("GEOS version:") or "C++ compiler:" in line
            ]
        (args.output_dir / "sedovResults.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Created four figures in {args.output_dir}.")
    print("Snapshot times (ms):", data["times"] * 1e3)
    print("Initial stress and motion, sampling times, finite values, and symmetry checks passed.")


if __name__ == "__main__":
    main()
