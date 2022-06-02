#!/usr/bin/python3
# Could easily parameterize size/spacing of cuboctahedra by adding a scale
# argument to make_cuboctahedron, make_octahedron, and get_spaced_centers
# and adjusting get_all_centers start points accordingly

import os
import numpy as np
import argparse


def make_cuboctahedron(center):
    cx, cy, cz = center
    # Define vertices
    p0 = (cx, cy - 1, cz - 1)
    p1 = (cx + 1, cy, cz - 1)
    p2 = (cx, cy + 1, cz - 1)
    p3 = (cx - 1, cy, cz - 1)
    p4 = (cx + 1, cy - 1, cz)
    p5 = (cx + 1, cy + 1, cz)
    p6 = (cx - 1, cy + 1, cz)
    p7 = (cx - 1, cy - 1, cz)
    p8 = (cx, cy - 1, cz + 1)
    p9 = (cx + 1, cy, cz + 1)
    p10 = (cx, cy + 1, cz + 1)
    p11 = (cx - 1, cy, cz + 1)
    vertices = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]

    # Define faces
    faces = [
        (p0, p3, p2, p1),
        (p0, p4, p8, p7),
        (p1, p5, p9, p4),
        (p2, p6, p10, p5),
        (p3, p7, p11, p6),
        (p8, p9, p10, p11),
        (p0, p1, p4),
        (p1, p2, p5),
        (p2, p3, p6),
        (p0, p7, p3),
        (p4, p9, p8),
        (p5, p10, p9),
        (p6, p11, p10),
        (p7, p8, p11),
    ]

    return vertices, faces


def make_octahedron(center):
    cx, cy, cz = center
    # Define vertices
    p0 = (cx, cy, cz - 1)
    p1 = (cx, cy - 1, cz)
    p2 = (cx + 1, cy, cz)
    p3 = (cx, cy + 1, cz)
    p4 = (cx - 1, cy, cz)
    p5 = (cx, cy, cz + 1)
    vertices = [p0, p1, p2, p3, p4, p5]

    # Define faces
    faces = [
        (p0, p2, p1),
        (p0, p3, p2),
        (p0, p4, p3),
        (p0, p1, p4),
        (p1, p2, p5),
        (p2, p3, p5),
        (p3, p4, p5),
        (p1, p5, p4),
    ]

    return vertices, faces


def get_spaced_centers(nx, ny, nz, start_point=(1, 1, 1)):
    # add check for nx, ny, nz > 0
    x0, y0, z0 = start_point
    total_points = nx * ny * nz
    grid = np.mgrid[
        x0 : (x0 + 2 * nx) : 2, y0 : (y0 + 2 * ny) : 2, z0 : (z0 + 2 * nz) : 2
    ]
    centers = np.vstack(list(map(np.ravel, grid))).T

    return tuple(map(tuple, centers))


def get_all_centers(nx, ny, nz):
    cuboctahedral = get_spaced_centers(nx, ny, nz, (1, 1, 1))
    octahedral = get_spaced_centers(nx - 1, ny - 1, nz - 1, (2, 2, 2))

    return cuboctahedral, octahedral


def separate_nodes_faces(cub_list, oct_list):
    nodes = list()
    nodes += (cubeocta[0] for cubeocta in cub_list)
    nodes += (octa[0] for octa in oct_list)
    faces = list()
    faces += (cubeocta[1] for cubeocta in cub_list)
    faces += (octa[1] for octa in oct_list)

    return nodes, faces


def map_nodes_to_inds(nodes):
    unique_points = set()

    for cell in nodes:
        unique_points.update(cell)

    index_map = {}

    for i, point in enumerate(unique_points):
        index_map[point] = i

    return index_map


def unpack_points(index_map):
    unpacked = list()

    for point, _ in sorted(index_map.items(), key=lambda x: x[1]):
        unpacked.extend(point)

    return unpacked


def translate_nodes(nodes, index_map):
    translated = list()
    offsets = list()
    count = 0

    for cell in nodes:
        count += len(cell)
        translated.extend([index_map[point] for point in cell])
        offsets.append(count)

    return translated, offsets


def translate_faces(faces, index_map):
    translated = list()
    offsets = list()
    count = 0

    for facelist in faces:
        temp_list = [len(facelist)]
        count += 1

        for face in facelist:
            temp_list.append(len(face))
            temp_list.extend([index_map[point] for point in face])
            count += len(face) + 1
        translated.extend(temp_list)
        offsets.append(count)

    return translated, offsets


def write_header(outfile):
    outfile.write('<?xml version="1.0"?>\n')


def open_node(outfile, space_count, name, **kwargs):
    arg_string = ""

    for arg, value in kwargs.items():
        arg_string += f' {arg}="{value}"'

    outfile.write(space_count * " " + f"<{name}" + arg_string + ">\n")
    space_count += 2

    return space_count


def close_node(outfile, space_count, name):
    space_count -= 2
    outfile.write(space_count * " " + f"</{name}>\n")

    return space_count


def write_data(outfile, space_count, data, **kwargs):
    space_count = open_node(outfile, space_count, "DataArray", **kwargs)
    outfile.write(space_count * " ")

    for index, value in enumerate(data, 1):
        outfile.write(f" {value}")
        # Write a newline every 15 values

        if (index % 15) == 0:
            outfile.write("\n" + space_count * " ")

    outfile.write("\n")
    space_count = close_node(outfile, space_count, "DataArray")


def print_to_vtu(fname, nodes, faces, index_map):
    # Precompute necessary values
    npoints = len(index_map)
    ncells = len(nodes)
    # Assign value 0 to cuboctahedra and 1 to octahedra for visualization
    cell_data = [0 if len(cell) == 12 else 1 for cell in nodes]
    points = unpack_points(index_map)
    cell_connectivity, cell_offsets = translate_nodes(nodes, index_map)
    cell_types = [42 for _ in range(ncells)]
    face_connectivity, face_offsets = translate_faces(faces, index_map)
    # Write to file
    with open(fname, "w") as outfile:
        write_header(outfile)
        space_count = 2
        space_count = open_node(
            outfile,
            space_count,
            "VTKFile",
            type="UnstructuredGrid",
            version="0.1",
            byte_order="LittleEndian",
        )
        space_count = open_node(outfile, space_count, "UnstructuredGrid")
        space_count = open_node(
            outfile, space_count, "Piece", NumberOfPoints=npoints, NumberOfCells=ncells
        )
        # Optional: add PointData later
        space_count = open_node(
            outfile, space_count, "CellData", Scalars="cell_scalars"
        )
        write_data(
            outfile,
            space_count,
            cell_data,
            type="Int32",
            Name="cell_scalars",
            format="ascii",
        )
        space_count = close_node(outfile, space_count, "CellData")
        space_count = open_node(outfile, space_count, "Points")
        write_data(
            outfile,
            space_count,
            points,
            type="Float32",
            NumberOfComponents=3,
            format="ascii",
        )
        space_count = close_node(outfile, space_count, "Points")
        space_count = open_node(outfile, space_count, "Cells")
        write_data(
            outfile,
            space_count,
            cell_connectivity,
            type="Int64",
            IdType=1,
            Name="connectivity",
            format="ascii",
        )
        write_data(
            outfile,
            space_count,
            cell_offsets,
            type="Int64",
            IdType=1,
            Name="offsets",
            format="ascii",
        )
        write_data(
            outfile, space_count, cell_types, type="UInt8", Name="types", format="ascii"
        )
        write_data(
            outfile,
            space_count,
            face_connectivity,
            type="Int64",
            IdType=1,
            Name="faces",
            format="ascii",
        )
        write_data(
            outfile,
            space_count,
            face_offsets,
            type="Int64",
            IdType=1,
            Name="faceoffsets",
            format="ascii",
        )
        space_count = close_node(outfile, space_count, "Cells")
        space_count = close_node(outfile, space_count, "Piece")
        space_count = close_node(outfile, space_count, "UnstructuredGrid")
        space_count = close_node(outfile, space_count, "VTKFile")


def make_cuboctahedral_mesh(nx, ny, nz, fname):
    # Generate centers, places centers in positive octant by default
    cub_centers, oct_centers = get_all_centers(nx, ny, nz)

    # Generate cuboctahedra
    cub_list = list(map(make_cuboctahedron, cub_centers))

    # Generate octahedra
    oct_list = list(map(make_octahedron, oct_centers))

    # Convert to node lists and face lists
    nodes, faces = separate_nodes_faces(cub_list, oct_list)

    # Make map of points to indices
    index_map = map_nodes_to_inds(nodes)

    # Write to file
    print_to_vtu(fname, nodes, faces, index_map)


if __name__ == "__main__":
    # Read input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "nx", help="Number of cuboctahedra in the x direction", type=int
    )
    parser.add_argument(
        "ny", help="Number of cuboctahedra in the y direction", type=int
    )
    parser.add_argument(
        "nz", help="Number of cuboctahedra in the z direction", type=int
    )
    parser.add_argument("outfile", help="Output file name")
    args = parser.parse_args()
    nx, ny, nz = (args.nx, args.ny, args.nz)
    fname = args.outfile
    make_cuboctahedral_mesh(nx, ny, nz, fname)
