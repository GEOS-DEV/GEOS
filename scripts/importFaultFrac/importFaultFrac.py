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


def id_list_to_tuple(i):
    o = []
    for v in range(i.GetNumberOfIds()):
        o.append(i.GetId(v))
    return tuple(o)


def duplicate_nodes(grid, cells_in_contact, internal_fracture_nodes):
    points = grid.GetPoints()

    is_cell_in_contact = numpy.zeros(grid.GetNumberOfCells(), dtype=bool)
    for cell_index in cells_in_contact:
        is_cell_in_contact[cell_index] = True
    is_fracture_node = numpy.zeros(grid.GetNumberOfPoints(), dtype=bool)
    for node_index in internal_fracture_nodes:
        is_fracture_node[cell_index] = True

    # internal fracture nodes to cells
    if_nodes_to_cells = defaultdict(set)
    for n in internal_fracture_nodes:
        for nc in cells_in_contact:
            if n in id_list_to_tuple(grid.GetCell(nc).GetPointIds()):
                if_nodes_to_cells[n].add(nc)
    # internal fracture nodes to triangles
    if_nodes_to_triangles = defaultdict(set)
    for n in internal_fracture_nodes:
        for nc in range(grid.GetNumberOfCells()):
            cell = grid.GetCell(nc)
            if cell.GetCellDimension() == 2:
                if n in id_list_to_tuple(cell.GetPointIds()):
                    if_nodes_to_triangles[n].add(nc)

    # new_cells = vtk.vtkCellArray()
    # new_cells.DeepCopy(grid.GetCells())
    # # new_cells.AllocateCopy(grid.GetCells())  # too much
    # for i, contact in enumerate(is_cell_in_contact):
    #     cell = grid.GetCell(cell)
    #     if not contact:
    #         # copy old cell
    #         new_cells.append(cell)
    #     else:

    # # creating another mesh from some components of the input mesh.
    # output = vtk.vtkUnstructuredGrid()
    # output.SetPoints(points)
    # output.SetCells(output_cell_types, output_cells)
    # output.GetCellData().AddArray(att)
    #
    # return output


def find_neighboring_cells(grid):
    """
    Find the cells that are in contact of the fracture/fault. Computes all the nodes that need to be split.
    :param grid:
    :return: cells_in_contact, internal_fracture_nodes
    """
    # output = grid

    nodes = grid.GetPoints()
    # fracture_nodes contains all the nodes belonging to fracture cells.
    is_fracture_nodes = numpy.zeros(nodes.GetNumberOfPoints(), dtype=bool)
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        if cell.GetCellDimension() == 2:
            for j in range(cell.GetNumberOfPoints()):
                is_fracture_nodes[cell.GetPointId(j)] = 1

    # Mapping from the cell id to the sorted nodes
    cell_to_nodes = [set()] * grid.GetNumberOfCells()
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        if cell.GetCellDimension() == 3:
            cell_to_nodes[i] = set(map(cell.GetPointId, range(cell.GetNumberOfPoints())))

    # Keeping the cells which touch the fracture cells with 3 nodes or more.
    cells_in_contact = []
    for cell_index, cell_nodes in enumerate(cell_to_nodes):
        num_connected_nodes = 0
        for n in cell_nodes:
            num_connected_nodes += is_fracture_nodes[n]
        if num_connected_nodes > 2:
            cells_in_contact.append(cell_index)

    # # For each triangle, returns the points of the face
    # all_triangle_points = {}
    # for i in range(grid.GetNumberOfCells()):
    #     if grid.GetCellType(i) == vtk.VTK_TRIANGLE:
    #         triangle = grid.GetCell(i)
    #         all_triangle_points[i] = triangle.GetPointId(0), triangle.GetPointId(1), triangle.GetPointId(2)

    # # progress = dict() # use a numpy full of -1
    # triangle_cell_neighbors = defaultdict(list)
    # for triangle, triangle_points in all_triangle_points.items():
    #     print(f"Searching face {triangle}")
    #     triangle_points_set = set(triangle_points)
    #     counter = 0
    #     for cell_index in cells_in_contact:
    #         cell_nodes = cell_to_nodes[cell_index]
    #         if cell_nodes.issuperset(triangle_points_set):
    #             # TODO this is good for tets, not all the types of elements
    #             print(f"    tet {cell_index} contains the face")
    #             triangle_cell_neighbors[triangle].append(cell_index)
    #             if len(triangle_cell_neighbors[triangle]) == 2:
    #                 break
    # assert set(map(lambda s: len(s), triangle_cell_neighbors.values())) == {2}

    # Now it's necessary to find which nodes need to be split.
    # The nodes on the boundary do not need to be split.
    fracture_edges = defaultdict(int)  # TODO master the hash...
    # all_fracture_edges = []
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        if cell.GetCellDimension() == 2:
            for ie in range(cell.GetNumberOfEdges()):
                edge = cell.GetEdge(ie)
                edge_nodes = []
                for ien in range(edge.GetNumberOfPoints()):  # TODO, do less pedantic
                    edge_nodes.append(edge.GetPointId(ien))
                # print(tuple(sorted(edge_nodes)))
                # edge_nodes = tuple(sorted(edge_nodes))
                # if edge_nodes[0] == 402 and edge_nodes[1] == 403:
                #     print("HELLO")
                # all_fracture_edges.append(edge_nodes)
                fracture_edges[tuple(sorted(edge_nodes))] += 1
    # Edges that are mentioned once are boundary edges.
    external_fracture_nodes = []
    for kv in filter(lambda i: i[1] == 1, fracture_edges.items()):
        print(kv)
        for n in kv[0]:
            external_fracture_nodes.append(n)
        # external_fracture_nodes += kv[0]
    external_fracture_nodes = set(external_fracture_nodes)
    all_fracture_nodes = []
    for k in fracture_edges.keys():
        all_fracture_nodes += k
    all_fracture_nodes = set(all_fracture_nodes)

    # Warning, we're loosing multiplicities here of edges. Is it important?
    internal_fracture_nodes = all_fracture_nodes - external_fracture_nodes

    return cells_in_contact, internal_fracture_nodes

    # TODO Then split the nodes! (the nodes inside the fracture, not at the boundary) vtkMarkBoundaryFilter?
    # Done


# In [190]: for n in external_fracture_nodes:
#     ...:     pt = output_mesh.GetPoint(n)
#     ...:     print ( pt[0], pt[1], pt[2], n )


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
    cells_in_contact, internal_fracture_nodes = find_neighboring_cells(output_mesh)

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

