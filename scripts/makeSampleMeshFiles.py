#!/usr/bin/python3

import argparse
import numpy as np
import vtk


def perturb(point_array, index_set):
    # Option for later: parameterize spacing and change perturbation size
    nrows = np.size(point_array, 0)

    if not index_set:
        return np.zeros((nrows, 3))

    x_nudge = np.zeros((nrows, 1))
    y_nudge = np.zeros((nrows, 1))
    z_nudge = np.zeros((nrows, 1))

    # Make each perturbation between +/-0.1, uniform distribution
    max_nudge = 0.1
    min_nudge = -1 * max_nudge

    if "x" in index_set:
        x_nudge = np.random.uniform(min_nudge, max_nudge, size=(nrows, 1))

    if "y" in index_set:
        y_nudge = np.random.uniform(min_nudge, max_nudge, size=(nrows, 1))

    if "z" in index_set:
        z_nudge = np.random.uniform(min_nudge, max_nudge, size=(nrows, 1))

    nudges = np.concatenate((x_nudge, y_nudge, z_nudge), axis=1)

    return nudges


def get_index_map(ideal, actual):
    indmap = dict()
    ideal_tuples = [tuple(point) for point in ideal]
    actual_tuples = [tuple(point) for point in actual]

    for ind, coord in enumerate(ideal_tuples):
        indmap[ind] = actual_tuples.index(coord)

    return indmap


def get_cube_indices(ncubes, indmap):
    index_list = []
    starting_corners = np.vstack(list(map(np.ravel, np.mgrid[0:ncubes, 0:ncubes, 0:ncubes]))).T

    for x, y, z in starting_corners:
        points = [
            (x, y, z),
            (x + 1, y, z),
            (x + 1, y + 1, z),
            (x, y + 1, z),
            (x, y, z + 1),
            (x + 1, y, z + 1),
            (x + 1, y + 1, z + 1),
            (x, y + 1, z + 1),
        ]
        cube_inds = []

        for xi, yi, zi in points:
            ind = (ncubes + 1) * (ncubes + 1) * xi + (ncubes + 1) * yi + zi
            cube_inds.append(indmap[ind])
        index_list.append(cube_inds)

    return index_list


def make_nxnxn_framework(n):
    # Option for later: make axbxc rectangle instead of nxnxn cube

    # Generate all categories of points
    corner_combos = np.mgrid[0 : n + 1 : n, 0 : n + 1 : n, 0 : n + 1 : n]
    corners = np.vstack(list(map(np.ravel, corner_combos))).T
    x_edge_combos = np.mgrid[1:n, 0 : n + 1 : n, 0 : n + 1 : n]
    x_edges = np.vstack(list(map(np.ravel, x_edge_combos))).T
    y_edge_combos = np.mgrid[0 : n + 1 : n, 1:n, 0 : n + 1 : n]
    y_edges = np.vstack(list(map(np.ravel, y_edge_combos))).T
    z_edge_combos = np.mgrid[0 : n + 1 : n, 0 : n + 1 : n, 1:n]
    z_edges = np.vstack(list(map(np.ravel, z_edge_combos))).T
    xy_face_combos = np.mgrid[1:n, 1:n, 0 : n + 1 : n]
    xy_faces = np.vstack(list(map(np.ravel, xy_face_combos))).T
    yz_face_combos = np.mgrid[0 : n + 1 : n, 1:n, 1:n]
    yz_faces = np.vstack(list(map(np.ravel, yz_face_combos))).T
    xz_face_combos = np.mgrid[1:n, 0 : n + 1 : n, 1:n]
    xz_faces = np.vstack(list(map(np.ravel, xz_face_combos))).T
    center_combos = np.mgrid[1:n, 1:n, 1:n]
    center = np.vstack(list(map(np.ravel, center_combos))).T

    # Organize points into list, set corresponding axes for nudging
    points_list = [corners, x_edges, y_edges, z_edges, xy_faces, yz_faces, xz_faces, center]
    points = np.concatenate(points_list, axis=0)
    nudge_axes = [
        set(),
        set(["x"]),
        set(["y"]),
        set(["z"]),
        set(["x", "y"]),
        set(["y", "z"]),
        set(["x", "z"]),
        set(["x", "y", "z"]),
    ]

    # Get an orderly array of all points and calculate cube corner indices
    ideal_points = np.vstack(list(map(np.ravel, np.mgrid[0 : n + 1, 0 : n + 1, 0 : n + 1]))).T
    # Increments z, then y, then x
    indmap = get_index_map(ideal_points, points)
    cube_indices = get_cube_indices(n, indmap)

    # Calculate perturbations to original points
    perturbation_list = []

    for i in range(0, len(points_list)):
        perturbation_list.append(perturb(points_list[i], nudge_axes[i]))

    # print(perturbation_list)
    perturbations = np.concatenate(perturbation_list, axis=0)
    final_points = points + perturbations

    return (final_points, cube_indices)


def make_points_object(point_array):
    points = vtk.vtkPoints()

    for x, y, z in point_array:
        points.InsertNextPoint(x, y, z)

    return points


def make_cell_array(index_list):
    cells = vtk.vtkCellArray()
    n_nodes = len(index_list[0])

    for node_list in index_list:
        if n_nodes == 8:
            cell = vtk.vtkHexahedron()
        elif n_nodes == 4:
            cell = vtk.vtkTetra()
        elif n_nodes == 5:
            cell = vtk.vtkPyramid()
        elif n_nodes == 6:
            cell = vtk.vtkWedge()
        else:
            print(f"{n_nodes} nodes are not supported at this time")

        for ind, value in enumerate(node_list):
            cell.GetPointIds().SetId(ind, value)

        cells.InsertNextCell(cell)

    if n_nodes == 8:
        cell_type = vtk.VTK_HEXAHEDRON
    elif n_nodes == 4:
        cell_type = vtk.VTK_TETRA
    elif n_nodes == 5:
        cell_type = vtk.VTK_PYRAMID
    elif n_nodes == 6:
        cell_type = vtk.VTK_WEDGE

    return cells, cell_type


def make_tets_from_cubes(cube_list):
    tet_list = []

    for a, b, c, d, e, f, g, h in cube_list:
        tet_list.extend([[a, b, c, f], [a, c, d, h], [a, e, f, h], [c, f, g, h], [a, f, c, h]])

    return tet_list


def make_wedges_from_cubes(cube_list):
    wedge_list = []

    for a, b, c, d, e, f, g, h in cube_list:
        wedge_list.extend([[a, d, h, b, c, g], [a, h, e, b, g, f]])

    return wedge_list


def make_pyramids_from_cubes(cube_list):
    pyramid_list = []

    for a, b, c, d, e, f, g, h in cube_list:
        pyramid_list.extend([[a, e, f, b, h], [b, f, g, c, h], [a, b, c, d, h]])

    return pyramid_list


def write_vtu(points, inds, fname):
    vtk_points = make_points_object(points)
    vtk_cells, vtk_type = make_cell_array(inds)
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(vtk_points)
    ugrid.SetCells(vtk_type, vtk_cells)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(fname)
    writer.SetDataModeToAscii()
    writer.Update()


def write_msh(points, inds, fname):
    with open(fname, "w") as fout:
        fout.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n")
        print(f"{len(points)}", file=fout)

        for i, point in enumerate(points, start=1):
            # Note that .msh files are 1-indexed
            x, y, z = point
            print(i, x, y, z, file=fout)

        fout.write("$EndNodes\n$Elements\n")
        print(f"{len(inds)}", file=fout)

        for i, elem in enumerate(inds, start=1):
            fout.write(str(i) + " ")

            if len(elem) == 8:
                fout.write("5")
            elif len(elem) == 4:
                fout.write("4")
            elif len(elem) == 5:
                fout.write("7")
            elif len(elem) == 6:
                fout.write("6")
            else:
                print(f"{len(elem)} nodes is not supported at this time")
            fout.write(" 2 0 1")

            for index in elem:
                fout.write(" " + str(index + 1))
            fout.write("\n")
        fout.write("$EndElements")


def generate_all_types(n):
    points, inds = make_nxnxn_framework(n)
    # Write hexes
    hex_name = "hex_" + str(n)
    write_msh(points, inds, hex_name + ".msh")
    write_vtu(points, inds, hex_name + ".vtu")
    # Write tets
    tet_name = "tet_" + str(n)
    tet_inds = make_tets_from_cubes(inds)
    write_msh(points, tet_inds, tet_name + ".msh")
    write_vtu(points, tet_inds, tet_name + ".vtu")
    # Write wedges
    wedge_name = "wedge_" + str(n)
    wedge_inds = make_wedges_from_cubes(inds)
    write_msh(points, wedge_inds, wedge_name + ".msh")
    write_vtu(points, wedge_inds, wedge_name + ".vtu")
    # Write pyramids
    pyramid_name = "pyramid_" + str(n)
    pyramid_inds = make_pyramids_from_cubes(inds)
    write_msh(points, pyramid_inds, pyramid_name + ".msh")
    write_vtu(points, pyramid_inds, pyramid_name + ".vtu")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script that generates NxNxN meshes with a cubic bounding box. Non-corner nodes are perturbed from their structured locations and outputs are written as hexahedra, tetrahedra, wedges, and pyramids in both .vtu and .msh formats."
    )
    parser.add_argument(
        "sideLengths", metavar="N", type=int, nargs="+", help="an integer side length for mesh construction"
    )
    parser.add_argument(
        "--seed", type=int, nargs="?", const=5, default=5, help="Integer random number generator seed, default=5"
    )
    args = parser.parse_args()

    for ncubes in args.sideLengths:
        np.random.seed(seed)
        generate_all_types(ncubes)
