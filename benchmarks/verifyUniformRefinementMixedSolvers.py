#!/usr/bin/env python3
"""Check connected hex/prism FEM and conservative FVM through refinement.

The repository's extruded hybrid mesh provides independent connected geometry.
Only the correctness harness uses numpy/h5py; GEOS needs no Python VTK package.
"""

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
from xml.etree import ElementTree as ET

from verifyUniformRefinementInitialization import make_inputs


ID_BASE = 9007199254741001
MESH_BASE = "Problem/domain/MeshBodies/mesh/meshLevels/Level0/"


def read_fixture(path):
    """Read the committed ASCII fixture; reject other formats/shapes explicitly."""
    import numpy as np
    tokens = path.read_text().split()
    assert "ASCII" in tokens and "BINARY" not in tokens
    cursor = tokens.index("POINTS")
    count = int(tokens[cursor + 1])
    points = np.asarray(tokens[cursor + 3:cursor + 3 + 3 * count], dtype=float).reshape(-1, 3)
    cursor = tokens.index("CELLS")
    count, entries = map(int, tokens[cursor + 1:cursor + 3])
    cursor += 3
    start = cursor
    cells = []
    for _ in range(count):
        size = int(tokens[cursor])
        cells.append([int(p) for p in tokens[cursor + 1:cursor + 1 + size]])
        cursor += size + 1
    assert cursor - start == entries
    cursor = tokens.index("CELL_TYPES")
    assert int(tokens[cursor + 1]) == count
    types = np.asarray(tokens[cursor + 2:cursor + 2 + count], dtype=int)
    cursor = tokens.index("LOOKUP_TABLE")
    attributes = np.asarray(tokens[cursor + 2:cursor + 2 + count], dtype=float).astype(int)
    volumes = []
    for cell, kind in zip(cells, types, strict=True):
        arity = {12: 4, 15: 5, 16: 6}[int(kind)]
        assert len(cell) == 2 * arity
        bottom, top = points[cell[:arity]], points[cell[arity:]]
        assert np.all(bottom[:, 2] == bottom[0, 2]) and np.all(top[:, 2] == top[0, 2])
        np.testing.assert_array_equal(bottom[:, :2], top[:, :2])
        signed_area = sum(a[0] * b[1] - a[1] * b[0] for a, b in zip(bottom, np.roll(bottom, -1, axis=0))) / 2
        volume = abs(signed_area) * (top[0, 2] - bottom[0, 2])
        assert volume > 0
        volumes.append(volume)
    assert set(types) == {12, 15, 16}, "This oracle expects connected hexes and both prism sizes"
    return points, cells, types, attributes, np.asarray(volumes)


def write_mesh(path, fixture, pulse):
    points, cells, types, attributes, _ = fixture
    root = ET.Element("VTKFile", type="UnstructuredGrid", version="0.1", byte_order="LittleEndian")
    piece = ET.SubElement(ET.SubElement(root, "UnstructuredGrid"), "Piece",
                          NumberOfPoints=str(len(points)), NumberOfCells=str(len(cells)))

    def array(parent, values, kind="Float64", **kwargs):
        element = ET.SubElement(parent, "DataArray", type=kind, format="ascii", **kwargs)
        element.text = " ".join(str(value) for value in values)

    point_data = ET.SubElement(piece, "PointData", GlobalIds="pointIds")
    array(point_data, (ID_BASE + 7 * p for p in range(len(points))), "Int64", Name="pointIds", IdType="1")
    cell_data = ET.SubElement(piece, "CellData", GlobalIds="cellIds")
    array(cell_data, (ID_BASE + 1000 + 13 * c for c in range(len(cells))), "Int64", Name="cellIds", IdType="1")
    array(cell_data, attributes, "Int32", Name="attribute")
    array(cell_data, (100000 + (20000 * (int(attribute) + 1) if pulse else 0)
                      for attribute in attributes), Name="seedPressure")
    array(cell_data, (0.08 + 0.01 * int(attribute) + 0.000001 * c
                      for c, attribute in enumerate(attributes)), Name="seedPorosity")
    array(ET.SubElement(piece, "Points"), points.flat, NumberOfComponents="3")
    cell_arrays = ET.SubElement(piece, "Cells")
    array(cell_arrays, (p for cell in cells for p in cell), "Int64", Name="connectivity")
    total = 0
    offsets = []
    for cell in cells:
        total += len(cell)
        offsets.append(total)
    array(cell_arrays, offsets, "Int64", Name="offsets")
    array(cell_arrays, types, "UInt8", Name="types")
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def make_cases(directory, fixture, levels):
    # Reuse the established solver/material configuration, with a unit domain.
    make_inputs(directory, (1, 1, 1), (1, 1, 1))
    write_mesh(directory / "constant.vtu", fixture, False)
    write_mesh(directory / "pulse.vtu", fixture, True)
    paths = {}
    for case in ("constant", "pulse", "fem"):
        for level in levels:
            kind = "fem" if case == "fem" else "fvm"
            tree = ET.parse(directory / f"{kind}-L{level}.xml")
            root = tree.getroot()
            vtk = root.find("Mesh/VTKMesh")
            vtk.set("file", str(directory / ("pulse.vtu" if case == "pulse" else "constant.vtu")))
            root.find("ElementRegions/CellElementRegion").set("cellBlocks", "{ * }")
            if kind == "fvm":
                vtk.set("fieldsToImport", "{ seedPressure, seedPorosity }")
                vtk.set("fieldNamesInGEOS", "{ pressure, porosity_referencePorosity }")
                root.remove(root.find("FieldSpecifications"))
            path = directory / f"{case}-L{level}.xml"
            tree.write(path, encoding="utf-8", xml_declaration=True)
            paths[case, level] = path
    return paths


def locate_roots(centers, fixture):
    """Coarse cell containing each fine-cell center.

    The fixture cells are vertical extrusions, so a center lies in a cell when it
    is inside its z range and inside its bottom polygon (winding number test).
    Refinement keeps only the coarse root ID with the GEOS cells and does not
    write it to checkpoints, so the checks locate roots from geometry.
    """
    import numpy as np
    points, cells, types, _, _ = fixture
    roots = np.full(len(centers), -1, dtype=np.int64)
    for root, (cell, kind) in enumerate(zip(cells, types, strict=True)):
        arity = {12: 4, 15: 5, 16: 6}[int(kind)]
        bottom = points[cell[:arity]]
        low, high = bottom[0, 2], points[cell[arity], 2]
        x, y = centers[:, 0], centers[:, 1]
        winding = np.zeros(len(centers), dtype=np.int64)
        for a, b in zip(bottom, np.roll(bottom, -1, axis=0)):
            side = (b[0] - a[0]) * (y - a[1]) - (x - a[0]) * (b[1] - a[1])
            winding += ((a[1] <= y) & (b[1] > y) & (side > 0)).astype(np.int64)
            winding -= ((a[1] > y) & (b[1] <= y) & (side < 0)).astype(np.int64)
        inside = (winding != 0) & (centers[:, 2] > min(low, high)) & (centers[:, 2] < max(low, high))
        assert not np.any(inside & (roots >= 0)), "A fine cell lies in two coarse cells"
        roots[inside] = root
    assert np.all(roots >= 0), "Every fine cell must lie in a coarse cell"
    return roots


def values(group):
    """Read the repository wrapper's stored permutation, retaining exact integers."""
    raw = group["__values__"][:]
    if "__dimensions__" not in group:
        return raw
    dims = tuple(int(i) for i in group["__dimensions__"][:])
    perm = tuple(int(i) for i in group["__permutation__"][:])
    return raw.reshape(tuple(dims[i] for i in perm)).transpose(tuple(perm.index(i) for i in range(len(perm))))


def inspect(file, case, level, rank, fixture, records, step):
    import h5py
    import numpy as np
    points, cells, types, attributes, coarse_volumes = fixture
    result = {"owned": 0, "ghosts": 0, "volume": 0.0, "mass": 0.0,
              "initial_mass": 0.0, "pressure_change": 0.0, "affine_error": 0.0, "new_boundary_nodes": 0}
    with h5py.File(file) as data:
        assert int(values(data["Problem/Mesh/mesh/uniformRefinement"])[0]) == level
        mesh = data[MESH_BASE]
        node_ids = values(mesh["nodeManager/localToGlobalMap"])
        positions = values(mesh["nodeManager/ReferencePosition"])
        node_ghosts = values(mesh["nodeManager/ghostRank"])
        if case == "fem":
            displacement = values(mesh["nodeManager/totalDisplacement"])
            bulk, shear = 5.5556e9, 4.16667e9
            young = 9 * bulk * shear / (3 * bulk + shear)
            poisson = (3 * bulk - 2 * shear) / (2 * (3 * bulk + shear))
            # The final checkpoint is written after the ramp reaches time two.
            assert step == 2
            expected = positions * np.array([1, -poisson, -poisson]) * (2e6 / young)
            np.testing.assert_allclose(displacement, expected, rtol=2e-7, atol=2e-11, err_msg=str(file))
            result["affine_error"] = float(np.max(abs(displacement - expected), initial=0))
            on_boundary = np.zeros(len(node_ids), dtype=bool)
            for name, axis, coordinate in (("left", 0, 0), ("right", 0, 1), ("yZero", 1, 0), ("zZero", 2, 0)):
                selected = abs(positions[:, axis] - coordinate) <= 0.001
                np.testing.assert_array_equal(values(mesh["nodeManager/sets/" + name]), np.flatnonzero(selected),
                                               err_msg="Boundary sets must include the final fine nodes")
                on_boundary |= selected
            result["new_boundary_nodes"] = int(sum(on_boundary & (node_ghosts < 0) & (node_ids > ID_BASE + 7 * (len(points) - 1))))
            for i, gid in enumerate(node_ids):
                records["nodes"].setdefault(int(gid), []).append(
                    (int(node_ghosts[i]), rank, positions[i], displacement[i]))
        subregions = mesh["ElementRegions/elementRegionsGroup/region/elementSubRegions"]
        for subregion in subregions.values():
            if not isinstance(subregion, h5py.Group):
                continue
            gids = values(subregion["localToGlobalMap"])
            ghosts = values(subregion["ghostRank"])
            volumes = values(subregion["elementVolume"])
            centers = values(subregion["elementCenter"])
            assert np.all(volumes > 0)
            # Lineage is used during import only. GEOS writes a NO_WRITE wrapper as an
            # empty placeholder, so no lineage entry may carry values.
            assert not any(name.startswith("_geosUniform") and isinstance(subregion[name], h5py.Group) for name in subregion)
            if level:
                root_indices = locate_roots(centers, fixture)
            else:
                root_indices = (gids - ID_BASE - 1000) // 13
                assert np.all((root_indices >= 0) & (root_indices < len(cells)))
                np.testing.assert_array_equal(gids, ID_BASE + 1000 + 13 * root_indices)
                np.testing.assert_array_equal(root_indices, locate_roots(centers, fixture))
            field = values(subregion["pressure"]) if case != "fem" else None
            if case != "fem":
                porosity = values(subregion["porosity/referencePorosity"])
                porosity = porosity.reshape(len(gids), -1)[:, 0] if len(gids) else porosity.reshape(0)
                expected_porosity = 0.08 + 0.01 * attributes[root_indices] + 0.000001 * root_indices
                np.testing.assert_allclose(porosity, expected_porosity, rtol=0, atol=1e-15,
                                           err_msg="Imported material field lost source-cell alignment")
                initial_pressure = np.full(len(gids), 100000) + (20000 * (attributes[root_indices] + 1) if case == "pulse" else 0)
                if step == 0:
                    np.testing.assert_array_equal(values(subregion["pressure_n"]), initial_pressure)
                if case == "constant":
                    np.testing.assert_allclose(field, initial_pressure, rtol=0, atol=1e-8)
                result["pressure_change"] = max(result["pressure_change"], float(np.max(abs(field - initial_pressure), initial=0)))
                mask = ghosts < 0
                mass = values(subregion["mass"])
                result["mass"] += float(sum(mass[mask]))
                result["initial_mass"] += float(sum(volumes[mask] * porosity[mask] * (1 + 1e-9 * initial_pressure[mask]) *
                                                    1000 * np.exp(1e-9 * initial_pressure[mask])))
            result["owned"] += int(sum(ghosts < 0))
            result["ghosts"] += int(sum(ghosts >= 0))
            result["volume"] += float(sum(volumes[ghosts < 0]))
            for i, gid in enumerate(gids):
                records["cells"].setdefault(int(gid), []).append((int(ghosts[i]), rank, centers[i], field[i] if field is not None else None))
                if ghosts[i] < 0:
                    root = int(root_indices[i])
                    records["root_counts"][root] += 1
                    records["root_ranks"].setdefault(root, set()).add(rank)
                    records["root_volumes"][root] += float(volumes[i])
                    key = tuple(np.round(centers[i], 11))
                    assert key not in records["solution"]
                    if field is not None:
                        records["solution"][key] = field[i]
        # Original vertex namespaces survive unchanged, independent of partition count.
        coarse = {tuple(point): ID_BASE + 7 * i for i, point in enumerate(points)}
        for i, point in enumerate(positions):
            if tuple(point) in coarse:
                assert node_ids[i] == coarse[tuple(point)]
    return result


def verify_records(records, fixture, level):
    import numpy as np
    _, cells, types, _, volumes = fixture
    for root in range(len(cells)):
        expected = 1 if level == 0 else {12: 8, 15: 10, 16: 12}[int(types[root])] * 8**(level - 1)
        assert records["root_counts"][root] == expected, root
        # Every child keeps the owner of its coarse root.
        assert len(records["root_ranks"][root]) == 1, (root, "children of one root on several ranks")
        np.testing.assert_allclose(records["root_volumes"][root], volumes[root], rtol=1e-11, atol=1e-14)
    for collection in (records["nodes"], records["cells"]):
        for gid, replicas in collection.items():
            owners = [record for record in replicas if record[0] < 0]
            assert len(owners) == 1, (gid, "global ownership")
            _, rank, position, field = owners[0]
            for ghost_rank, replica_rank, replica_position, replica_field in replicas:
                if ghost_rank >= 0:
                    assert ghost_rank == rank and replica_rank != rank
                np.testing.assert_allclose(replica_position, position, rtol=0, atol=1e-12)
                if field is not None:
                    np.testing.assert_allclose(replica_field, field, rtol=1e-12, atol=1e-8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path, help="new directory for inputs, checkpoints and evidence")
    parser.add_argument("--launcher", default="mpiexec")
    parser.add_argument("--ranks", nargs="+", type=int, default=[1, 4])
    parser.add_argument("--levels", nargs="+", type=int, choices=[0, 1, 2], default=[0, 1, 2])
    parser.add_argument("--cases", nargs="+", choices=["constant", "pulse", "fem"], default=["constant", "pulse", "fem"])
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--verify-only", action="store_true", help="inspect completed runs in output without launching GEOS")
    args = parser.parse_args()
    if not __debug__:
        parser.error("run without Python -O; checks use assertions")
    if 1 not in args.ranks or min(args.ranks) < 1 or args.timeout <= 0:
        parser.error("ranks must include one and be positive; timeout must be positive")
    import numpy as np
    import h5py
    del h5py
    executable = args.executable.resolve(strict=True)
    launcher = shlex.split(args.launcher)
    if not launcher:
        parser.error("launcher must not be empty")
    directory = args.output.resolve()
    source = Path(__file__).resolve().parents[1] / "inputFiles/solidMechanics/hybridHexPrismMesh.vtk"
    fixture = read_fixture(source)
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    if args.verify_only:
        assert json.loads((directory / "fixture.json").read_text())["sha256"] == source_sha
        completed = {row["case"]: row["returncode"] for row in json.loads((directory / "status.json").read_text())}
        paths = {(case, level): directory / f"{case}-L{level}.xml" for case in args.cases for level in args.levels}
    else:
        directory.mkdir(parents=True, exist_ok=False)
        paths = make_cases(directory, fixture, args.levels)
        (directory / "fixture.json").write_text(json.dumps({"source": str(source), "sha256": source_sha,
                                                           "cell_types": dict(Counter(map(int, fixture[2]))), "coarse_volume": float(sum(fixture[-1]))}, indent=2) + "\n")
    checks, status, reference = [], [], {}
    for ranks in sorted(set(args.ranks)):
        for case in args.cases:
            for level in sorted(set(args.levels)):
                name = f"{case}-L{level}-p{ranks}"
                out = directory / name
                command = [str(executable), "-i", str(paths[case, level]), "-o", str(out)]
                if ranks > 1 or Path(launcher[0]).name == "srun":
                    command = launcher + ["-n", str(ranks)] + command
                log_path = directory / f"{name}.log"
                if args.verify_only:
                    returncode = completed[name]
                else:
                    out.mkdir()
                    with log_path.open("w") as log:
                        process = subprocess.Popen(command, stdout=log, stderr=log, start_new_session=True)
                        try:
                            returncode = process.wait(timeout=args.timeout)
                        except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
                            try:
                                os.killpg(process.pid, signal.SIGTERM)
                            except ProcessLookupError:
                                pass
                            try:
                                process.wait(timeout=10)
                            except subprocess.TimeoutExpired:
                                os.killpg(process.pid, signal.SIGKILL)
                                process.wait()
                            raise error
                status.append({"case": name, "returncode": returncode, "command": command})
                if not args.verify_only:
                    (directory / "status.json").write_text(json.dumps(status, indent=2) + "\n")
                assert returncode == 0, f"{name}: see {log_path}"
                if case != "constant":
                    assert "linear solver solve time" in log_path.read_text(), "The oracle must follow an actual solve"
                final_records = None
                for step in ([2] if case == "fem" else [0, 1, 2]):
                    records = {"cells": {}, "nodes": {}, "root_counts": Counter(), "root_volumes": Counter(),
                               "root_ranks": {}, "solution": {}}
                    totals = Counter()
                    for rank in range(ranks):
                        file = out / f"{case}-L{level}_restart_{step:09}" / f"rank_{rank:07}.hdf5"
                        result = inspect(file, case, level, rank, fixture, records, step)
                        for key, value in result.items():
                            totals[key] = max(totals[key], value) if key in ("pressure_change", "affine_error") else totals[key] + value
                    verify_records(records, fixture, level)
                    np.testing.assert_allclose(totals["volume"], sum(fixture[-1]), rtol=1e-11, atol=1e-14)
                    if ranks > 1:
                        assert totals["ghosts"] > 0
                    if case == "fem" and level > 0:
                        assert totals["new_boundary_nodes"] > 0
                    if case != "fem":
                        # Closed-boundary fluxes cancel globally despite nonzero internal fluxes.
                        np.testing.assert_allclose(totals["mass"], totals["initial_mass"], rtol=2e-11, atol=1e-12)
                        if case == "pulse":
                            assert totals["pressure_change"] > 1, "The conservation check must exercise a nonzero solve"
                    checks.append({"case": name, "step": step, **totals})
                    final_records = records
                if case != "fem":
                    solution = final_records["solution"]
                    if ranks == 1:
                        reference[case, level] = solution
                    else:
                        expected = reference[case, level]
                        assert solution.keys() == expected.keys()
                        np.testing.assert_allclose(list(solution.values()), [expected[key] for key in solution], rtol=2e-10, atol=2e-5)
                (directory / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
                print(f"Passed {name}", flush=True)
    print(f"Passed {len(status)} mixed-mesh application runs; evidence: {directory / 'verification.json'}")


if __name__ == "__main__":
    main()
