import argparse
import logging
from typing import List
import sys
# import itertools

from collections import defaultdict

import numpy

import vtk
import vtk.util
import vtk.util.numpy_support


def cyclic_perm(a):
    n = len(a)
    b = [[a[i - j] for i in range(n)] for j in range(n)]
    return b


def find_neighboring_cells(grid):
    output = grid

    nodes = grid.GetPoints()
    # fracture_nodes contains all the nodes belonging to fracture cells.
    fracture_nodes = numpy.zeros(nodes.GetNumberOfPoints(), dtype=bool)
    for i in range(output.GetNumberOfCells()):
        cell = grid.GetCell(i)
        if cell.GetCellDimension() == 2:
            for j in range(cell.GetNumberOfPoints()):
                fracture_nodes[cell.GetPointId(j)] = 1

    # Mapping from the cell id to the sorted nodes
    cell_to_nodes = [set()] * output.GetNumberOfCells()
    for i in range(output.GetNumberOfCells()):
        cell = grid.GetCell(i)
        if cell.GetCellDimension() == 3:
            cell_to_nodes[i] = set(map(cell.GetPointId, range(cell.GetNumberOfPoints())))

    # Keeping the cells which touch the fracture cells with 3 nodes or more.
    cells_in_contact = []
    for cell_index, cell_nodes in enumerate(cell_to_nodes):
        num_connected_nodes = 0
        for n in cell_nodes:
            num_connected_nodes += fracture_nodes[n]
        if num_connected_nodes > 2:
            cells_in_contact.append(cell_index)

    # For each triangle, returns the points of the face
    all_triangle_points = {}
    for i in range(output.GetNumberOfCells()):
        if grid.GetCellType(i) == vtk.VTK_TRIANGLE:
            triangle = grid.GetCell(i)
            all_triangle_points[i] = triangle.GetPointId(0), triangle.GetPointId(1), triangle.GetPointId(2)

    progress = dict() # use a numpy full of -1
    # used_cells = numpy.zeros(len(cells_in_contact), dtype=int)
    for triangle, triangle_points in all_triangle_points.items():
        # print(f"Searching face {triangle}")
        triangle_points_set = set(triangle_points)
        counter = 0
        for i, cell_index in enumerate(cells_in_contact):
            # if used_cells[i] < 2:
                cell_nodes = cell_to_nodes[cell_index]
                if cell_nodes.issuperset(triangle_points_set):
                    # TODO this is good for tets, not all the types of elements
                    # print(f"    tet {cell_index} contains the face")
                    counter += 1
                    # used_cells[i] += 1
                    progress[triangle] = counter
                    if counter == 2:
                        break
    # TODO Store the tet id and find which tet is on which side of the normal
    #      Then split the nodes! (the nodes inside the fracture, not at the boundary) vtkMarkBoundaryFilter?
    # Done

# without used_cells
# 8.18 s ± 175 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
# with used_cells
# 19.5 s ± 502 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

    progress = dict() # use a numpy full of -1
    for triangle, triangle_points in all_triangle_points.items():
        print(f"Searching face {triangle}")
        triangle_points_set = set(triangle_points)
        counter = 0
        for i, nodes in enumerate(cell_to_nodes):
            if nodes.issuperset(triangle_points_set):
                print(f"    tet {i} contains the face")
                counter += 19
                progress[triangle] = counter
                if counter == 2:
                    break
    # TODO Store the tet id and find which tet is on which side of the normal
    assert set(progress.values()) == {2}


    # # Finding the tets that contain all the three points of each fracture cell
    # progress = dict()
    # for triangle, triangle_points in all_triangle_points.items():
    #     print(f"Searching face {triangle}")
    #     for i in range(output.GetNumberOfCells()):
    #         if grid.GetCellType(i) == vtk.VTK_TETRA:
    #             counter = 0
    #             tet = grid.GetCell(i)
    #             tet_points = set(map(tet.GetPointId, (0, 1, 2, 3)))
    #             if tet_points.issuperset(set(triangle_points)):
    #                 print(f"    tet {i} contains the face")
    #                 counter += 1
    #                 progress[i] = counter
    #                 if counter == 2:
    #                     break


    # Finding the faces
    progress = defaultdict(int)
    for triangle, triangle_points in all_triangle_points.items():
        print(f"Searching face {triangle}")
        for i in range(output.GetNumberOfCells()):
            if grid.GetCellType(i) == vtk.VTK_TETRA:
                tet = grid.GetCell(i)
                # Warning, do not store instances to the faces since there is one single instance.
                # You would end up working with the same unique face, not all the faces of the tet.
                tet_face_points = []
                for j in 0, 1, 2, 3:
                    face = tet.GetFace(j)
                    tet_face_points.append((face.GetPointId(0), face.GetPointId(1), face.GetPointId(2)))
                for tet_face_point in tet_face_points:
                    for cyclic_tet_face_point in cyclic_perm(tet_face_point):
                        if tuple(cyclic_tet_face_point) == tuple(triangle_points):
                            print(f"    found tet {i}.")
                            progress[triangle] += 1
                    for cyclic_tet_face_point in cyclic_perm(list(reversed(tet_face_point))):
                        if tuple(cyclic_tet_face_point) == tuple(triangle_points):
                            print(f"    found reversed tet {i}.")
                            progress[triangle] += 2


    # for i in range(output.GetNumberOfCells()):
    #     if grid.GetCellType(i) == vtk.VTK_TETRA:
    #         tet = grid.GetCell(i)
    #         # Warning, do not store instances to the faces since there is one single instance.
    #         # You would end up working with the same unique face, not all the faces of the tet.
    #         face_points = []
    #         for j in 0, 1, 2, 3:
    #             face = tet.GetFace(j)
    #             face_points.append((face.GetPointId(0), face.GetPointId(1), face.GetPointId(2)))
    #         for face_points in face_points:
    #             # face_points = face.GetPointId(0), face.GetPointId(1), face.GetPointId(2)
    #             for triangle, points in triangle_points.items():
    #                 if face_points == points: # check permutations
    #                     print(f"Found face {triangle} for tet {i}.")
    #                     # print(face_points, points)



def clean_input_mesh(vtk_input_file: str):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_input_file)
    reader.Update()
    input_ = reader.GetOutput()
    # iCellTypes = vtk.vtkCellTypes()
    # input_.GetCellTypes(iCellTypes)
    # Keeping only a specific type of cells
    # keep_cells = {vtk.VTK_TRIANGLE: defaultdict(list), vtk.VTK_TETRA: defaultdict(list)}
    # vtk type -> zone -> cell id
    keep_cells = defaultdict(lambda: defaultdict(list))
    # keep_cell_types = keep_cells.keys()
    # Getting the physical field.
    zones = input_.GetCellData().GetArray('CellEntityIds')
    assert zones.GetNumberOfValues() == input_.GetNumberOfCells()

    for i in range(input_.GetNumberOfCells()):
        zone, cell_type = zones.GetValue(i), input_.GetCellType(i)
        # if cell_type in keep_cell_types:
        keep_cells[cell_type][zone].append(i)

    # const_keep_cells = dict(keep_cells)
    keep_cells = dict.copy(keep_cells) # const version
    # We only want to keep triangles and for faults and the tets.
    faults = {1, 4, 7, 10, 11, 12}
    num_tets = len(keep_cells[vtk.VTK_TETRA][1])
    num_triangles = 0
    for fault in faults:
        num_triangles += len(keep_cells[vtk.VTK_TRIANGLE][fault])
    # The new total of cells will be num_tets + num_triangles

    output_cells = vtk.vtkCellArray()
    output_cells.AllocateExact(num_tets + num_triangles, num_tets * 4 + num_triangles * 3)

    for tet_id in keep_cells[vtk.VTK_TETRA][1]:
        tet = input_.GetCell(tet_id)
        output_cells.InsertNextCell(tet)
    attributes = [0] * num_tets
    output_cell_types = [vtk.VTK_TETRA] * num_tets

    for region, triangle_ids in keep_cells[vtk.VTK_TRIANGLE].items():
        if region in faults:
            for triangle_id in triangle_ids:
                triangle = input_.GetCell(triangle_id)
                output_cells.InsertNextCell(triangle)
            attributes += [region] * len(triangle_ids) # keep same id, don't increment.
    output_cell_types += [vtk.VTK_TRIANGLE] * num_triangles

    assert len(attributes) == num_tets + num_triangles
    assert len(output_cell_types) == num_tets + num_triangles
    assert output_cells.GetNumberOfCells() == num_tets + num_triangles

    # Let's loop over all the cells and insert them
    # In the meantime, we should be creating the attribute array.
    # attributeCellData = vtk.vtkCellData()
    att = vtk.util.numpy_support.numpy_to_vtk(attributes)
    # att = vtk.vtkIntArray()
    att.SetName("attribute")
    # att.SetNumberOfComponents(1)
    # att.SetNumberOfTuples(num_tets + num_triangles)

    # creating another mesh from some components of the input mesh.
    output = vtk.vtkUnstructuredGrid()
    output.SetPoints(input_.GetPoints())
    output.SetCells(output_cell_types, output_cells)
    output.GetCellData().AddArray(att)

    return output


def parse(cli_args: List[str]):
    """
    Parse the command line arguments and return the corresponding structure.
    :param cli_args: The list of arguments (as strings)
    :return: The struct
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--vtk-input-file',
                        type=str,
                        dest="vtk_input_file")
    args = parser.parse_args(cli_args)
    return args


if __name__ == '__main__':
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s', level=logging.INFO)
    args = parse(sys.argv[1:])
    logging.info("Playing with hi24L coarse.")
    # vtk_input_file="/tmp/tmp.EtfQgEUXLg/tmp_buffer/HI24L_coars.vtk"
    vtk_input_file = "HI24L_coars.vtk"
    output_mesh = clean_input_mesh(vtk_input_file)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName("output.vtu")
    writer.SetInputData(output_mesh)
    writer.Write()


# triangles CellEntityIds
# faults 1, 4, 7, 10, 11, 12
# BC bottom/z- 2
# BC top/z+ 3
# BC y+ 5
# BC y- 6
# BC x+ 8
# BC x- 9

# tets CellEntityIds
# body 1 (unique value)

