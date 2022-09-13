import argparse
import logging
from typing import List
import sys
import itertools

from collections import defaultdict

import numpy
import networkx as nx

import vtk
import vtk.util
import vtk.util.numpy_support


def id_list_to_tuple(i):
    o = []
    for v in range(i.GetNumberOfIds()):
        o.append(i.GetId(v))
    return tuple(o)

def duplicate_fracture_nodes(grid, connected_components, internal_fracture_nodes):
    pass


def color_fracture_sides(grid, internal_fracture_nodes, external_fracture_nodes):
    # Build an internal fracture node indicator for fast lookup.
    is_internal_fracture_node = numpy.zeros(grid.GetNumberOfPoints(), dtype=bool)
    for n in internal_fracture_nodes:
        is_internal_fracture_node[n] = True

    # One big elem_to_nodes: elems can be of all kinds (both tets and triangles, 2d and 3d).
    elem_to_nodes = [set()] * grid.GetNumberOfCells()
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        elem_to_nodes[i] = set(map(cell.GetPointId, range(cell.GetNumberOfPoints())))

    # We want to know which cell is touching each internal fracture node.
    int_frac_nodes_to_cells = defaultdict(set)
    int_frac_nodes_to_triangles = defaultdict(set)
    for cell_index, cell_nodes in enumerate(elem_to_nodes):
        cell = grid.GetCell(cell_index)
        for n in cell_nodes:
            if is_internal_fracture_node[n]:
                if cell.GetCellDimension() == 3:
                    int_frac_nodes_to_cells[n].add(cell_index)
                if cell.GetCellDimension() == 2:
                    int_frac_nodes_to_triangles[n].add(cell_index)

    int_frac_nodes_to_triangle_nodes = defaultdict(list)
    for node, triangles in int_frac_nodes_to_triangles.items():
        for triangle in triangles:
            tr = grid.GetCell(triangle)
            triangle_nodes = {tr.GetPointId(0), tr.GetPointId(1), tr.GetPointId(2)}
            int_frac_nodes_to_triangle_nodes[node].append(triangle_nodes)

    # Now we must find the neighboring cells which do not cross the fracture
    # To do this, we find all the connections which cross.
    # Any other connection will be a valid connection
    graph = nx.Graph()
    for node, cells in int_frac_nodes_to_cells.items():
        print(f"Working on node {node}")
        # n_cells = len(cells)
        for ic0, ic1 in itertools.combinations(cells, 2):
            assert ic0 != ic1
            c0 = grid.GetCell(ic0) # This is not thread safe, use GetCells instead
            faces0 = set()
            for nf in range(c0.GetNumberOfFaces()):
                f = c0.GetFace(nf)
                ff = tuple(sorted(map(f.GetPointId, range(f.GetNumberOfPoints()))))
                faces0.add(ff)
            faces1 = set()
            c1 = grid.GetCell(ic1)
            for nf in range(c1.GetNumberOfFaces()):
                f = c1.GetFace(nf)
                ff = tuple(sorted(map(f.GetPointId, range(f.GetNumberOfPoints()))))
                faces1.add(ff)
            intersection = faces1.intersection(faces0)
            # print (faces0, faces1, intersection)
            if len(intersection) == 1:
                common_face = set(intersection.pop())
                if common_face in int_frac_nodes_to_triangle_nodes[node]:
                    print(f"Crossing fracture between {ic0} and {ic1}.")
                else:
                    print(f"There is a non fracture-crossing between {ic0} and {ic1}.")
                    graph.add_edge(ic0, ic1)  # TODO do not rely on network to remove duplicates. let's do it ourselves.
            else:
                print(f"No connection between {ic0} and {ic1}.")
            # else:
            #     print("link")
            # Find common face.
    # Returning the tuple of connected components.
    # It's somehow equivalent to the sides of the fractures.
    # Each element of the tuple is a set containing all the cell indices of each component.
    return tuple(nx.connected_components(graph))


# # Displaying the sides of the fractures
# for i, cc in enumerate(nx.connected_components(G)):
#     for cc_cell in cc:
#         cell = grid.GetCell(cc_cell)
#         pt_ids = id_list_to_tuple(cell.GetPointIds())
#         cell_center = numpy.zeros(3)
#         for pt_id in pt_ids:
#             cell_center += numpy.array(grid.GetPoint(pt_id))
#         cell_center /= len(pt_ids)
#         print(cell_center[0], cell_center[1], cell_center[2], i)


def extract_fracture_nodes(grid):
    """
    :param grid:
    :return: internal_fracture_nodes, external_fracture_nodes
    """
    # It's necessary to find which nodes need to be split.
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
        # print(kv)
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

    return internal_fracture_nodes, external_fracture_nodes

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
    internal_fracture_nodes, external_fracture_nodes = extract_fracture_nodes(output_mesh)
    connected_components = color_fracture_sides(internal_fracture_nodes, extract_fracture_nodes)

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

