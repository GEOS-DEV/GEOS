import argparse
import logging
import sys
from typing import List, Tuple, Iterable
import itertools
from collections import defaultdict

import numpy
import networkx as nx

import vtk
import vtk.util
import vtk.util.numpy_support


def id_list_to_tuple(i):
    """
    Utility function to convert a vtkIdList into a tuple.
    This is mostly for convenient interactive debugging
    :param i: The vtk id list.
    :return: The tuple.
    """
    o = []
    for v in range(i.GetNumberOfIds()):
        o.append(i.GetId(v))
    return tuple(o)


def duplicate_fracture_nodes(grid: vtk.vtkUnstructuredGrid,
                             connected_components: Iterable[Iterable[int]],
                             internal_fracture_nodes: Iterable[int]) -> vtk.vtkUnstructuredGrid:
    """
    Splits the mesh on the internal_fracture_nodes.
    New nodes are then created and inserted in the new mesh.
    All the already existing cells are updated to take into account these new nodes.
    The 2d elements of the fracture point to duplicated nodes,
    but there is no specific action to be sure that their nodes belong to the same 3d element.
    Node positions will be correct anyway.
    :param grid: the vtk unstructured mesh
    :param connected_components:
        Each component contains an iterable of cells that touch the internal fracture nodes.
        All the cells in one same bucket will be based on the same duplicated node.
    :param internal_fracture_nodes: Iterable of all the internal fracture nodes.
    :return:
        A new vtk unstructured mesh with duplicated nodes.
        The support of the cells take the duplicated nodes into account.
    """
    is_internal_fracture_node = numpy.zeros(grid.GetNumberOfPoints(), dtype=bool)
    for n in internal_fracture_nodes:
        is_internal_fracture_node[n] = True

    num_points = grid.GetNumberOfPoints()

    # We are creating a new mesh.
    # The cells will be the same, except that their nodes may be duplicated or renumbered nodes.
    # The `new_cells` array will be modified in place, while the original grid remains unmodified.
    new_cells = vtk.vtkCellArray()
    new_cells.DeepCopy(grid.GetCells())

    # Building an indicator to find fracture cells.
    is_fracture_side_cells = numpy.zeros(new_cells.GetNumberOfCells(), dtype=bool)
    for fracture_side_cells in connected_components:
        for fracture_side_cell in fracture_side_cells:
            is_fracture_side_cells[fracture_side_cell] = True

    # Each non duplicated nodes is renumbered.
    # The idea is that I do not want to fill holes in the numbering,
    # I just want to increment the next number by one.
    # A negative number means that this is a duplicated number.
    non_duplicated_nodes = numpy.ones(num_points, dtype=int) * -1
    next_point_id = 0  # When a node needs to be duplicated, `next_point_id` is the next available node.
    for node_idx in range(num_points):
        if not is_internal_fracture_node[node_idx]:
            non_duplicated_nodes[node_idx] = next_point_id
            next_point_id += 1

    # Connected component index -> duplicated nodes.
    # It's not super useful for the moment to store the precise binding with the component.
    component_to_dup_nodes = {}

    # First, dealing with the renumbering of the cell touching the fracture.
    for component_idx, fracture_side_cells in enumerate(connected_components):
        duplicated_nodes = dict()
        for fracture_side_cell in fracture_side_cells:
            # `new_cell_point_ids` will contain the ids for the new cells.
            # Calling `new_cells.GetCellAtId` is a trick to get it at the correct size,
            # because all values should be overwritten.
            new_cell_point_ids = vtk.vtkIdList()
            new_cells.GetCellAtId(fracture_side_cell, new_cell_point_ids)
            for i in range(new_cell_point_ids.GetNumberOfIds()):
                new_cell_point_id = new_cell_point_ids.GetId(i)
                # Here we find or build the new number of the node.
                if is_internal_fracture_node[new_cell_point_id]:
                    assert non_duplicated_nodes[new_cell_point_id] == -1
                    if new_cell_point_id in duplicated_nodes:  # if already duplicated, take the same
                        new_cell_point_ids.SetId(i, duplicated_nodes[new_cell_point_id])
                    else:  # otherwise, duplicate
                        duplicated_nodes[new_cell_point_id] = next_point_id
                        new_cell_point_ids.SetId(i, next_point_id)
                        next_point_id += 1
                else:  # Here, it's not an internal fracture node which was not duplicated.
                    assert not non_duplicated_nodes[new_cell_point_id] == -1
                    new_cell_point_ids.SetId(i, non_duplicated_nodes[new_cell_point_id])
            new_cells.ReplaceCellAtId(fracture_side_cell, new_cell_point_ids)
        component_to_dup_nodes[component_idx] = duplicated_nodes

    # Now, we must not forget to modify the nodes of all the other cells,
    # because the nodes of all the points have been modified!
    for cell_idx in range(new_cells.GetNumberOfCells()):
        if is_fracture_side_cells[cell_idx]:
            continue  # This is for non fracture 3d elements!
        # Same trick to get the indices list at the correct size.
        new_cell_point_ids = vtk.vtkIdList()
        new_cells.GetCellAtId(cell_idx, new_cell_point_ids)
        if grid.GetCell(cell_idx).GetCellDimension() == 3:  # For
            for i in range(new_cell_point_ids.GetNumberOfIds()):
                new_cell_point_id = new_cell_point_ids.GetId(i)
                assert not non_duplicated_nodes[new_cell_point_id] == -1
                new_cell_point_ids.SetId(i, non_duplicated_nodes[new_cell_point_id])
            new_cells.ReplaceCellAtId(cell_idx, new_cell_point_ids)
        elif grid.GetCell(cell_idx).GetCellDimension() == 2:
            for i in range(new_cell_point_ids.GetNumberOfIds()):
                new_cell_point_id = new_cell_point_ids.GetId(i)
                # For fracture triangles, all the points have been duplicated.
                # We have to decide which node we pick.
                # We do in some way, but it's surely a bit bad.
                if non_duplicated_nodes[new_cell_point_id] != -1:
                    duplicated_node = non_duplicated_nodes[new_cell_point_id]
                else:
                    for dup_node in component_to_dup_nodes.values():  # TODO take something consistent
                        duplicated_node = dup_node.get(new_cell_point_id)
                        if duplicated_node: break  # TODO what if it fails!
                assert duplicated_node is not None
                new_cell_point_ids.SetId(i, duplicated_node)
            new_cells.ReplaceCellAtId(cell_idx, new_cell_point_ids)
        else:
            assert False

    # Now we duplicated nodes and we remove duplicated nodes
    num_new_points = num_points - len(internal_fracture_nodes) + sum(map(len, component_to_dup_nodes.values()))
    assert next_point_id == num_new_points

    old_points = grid.GetPoints()
    new_points = vtk.vtkPoints()
    new_points.SetNumberOfPoints(num_new_points)

    # Those two sets are for validation only.
    moved_nodes_from = set()
    moved_nodes_to = set()

    for from_, to in enumerate(non_duplicated_nodes):
        if not to == -1:
            new_points.SetPoint(to, old_points.GetPoint(from_))
            assert to not in moved_nodes_to
            moved_nodes_from.add(from_)
            moved_nodes_to.add(to)

    for dup_nodes in component_to_dup_nodes.values():
        for from_, to in dup_nodes.items():
            new_points.SetPoint(to, old_points.GetPoint(from_))
            assert to not in moved_nodes_to
            moved_nodes_from.add(from_)
            moved_nodes_to.add(to)

    assert moved_nodes_from == set(range(num_points))
    assert moved_nodes_to == set(range(num_new_points))

    # creating another mesh from some components of the input mesh.
    output = vtk.vtkUnstructuredGrid()
    output.SetPoints(new_points)
    output.SetCells(grid.GetCellTypesArray(), new_cells)  # The cell types are unchanged; we reuse the old cell types!

    # Defining the regions of the mesh. Volume will be `0`, fracture 2d elements will be `1`.
    attributes = []
    for cell_idx in range(output.GetNumberOfCells()):
        if grid.GetCell(cell_idx).GetCellDimension() == 3:
            attributes.append(0)
        else:
            attributes.append(1)
    att = vtk.util.numpy_support.numpy_to_vtk(attributes)
    att.SetName("attribute")
    output.GetCellData().AddArray(att)

    return output


def color_fracture_sides(grid: vtk.vtkUnstructuredGrid,
                         internal_fracture_nodes: Iterable[int]) -> Iterable[Iterable[int]]:
    """
    Takes all the cells which have at least one in common with the fracture nodes,
    and separate them into connected buckets of cells.
    Then returns an iterable of bucket of cells touching the fracture.
    It's guaranteed that no bucket of the
    :param grid: the vtk unstructured mesh
    :param internal_fracture_nodes: all the nodes that are inside the fracture (not on the boundary of the fracture).
    :return: All the buckets connected fracture cells.
    """
    # Build an internal fracture node indicator for fast lookup.
    is_internal_fracture_node = numpy.zeros(grid.GetNumberOfPoints(), dtype=bool)
    for n in internal_fracture_nodes:
        is_internal_fracture_node[n] = True

    # One big elem_to_nodes map: elems can be of all kinds (both tets and triangles, 2d and 3d).
    elem_to_nodes = [None] * grid.GetNumberOfCells()
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        elem_to_nodes[i] = set(map(cell.GetPointId, range(cell.GetNumberOfPoints())))

    # We want to know which cell is touching each internal fracture node.
    int_frac_nodes_to_cells = defaultdict(set)  # mapping from internal fracture nodes to touching (at least 1 node) 3d elements
    int_frac_nodes_to_triangles = defaultdict(set)  # mapping from internal fracture nodes to touching (at least 1 node) 2d elements
    for cell_index, cell_nodes in enumerate(elem_to_nodes):
        cell = grid.GetCell(cell_index)
        for n in cell_nodes:
            if is_internal_fracture_node[n]:
                if cell.GetCellDimension() == 3:
                    int_frac_nodes_to_cells[n].add(cell_index)
                if cell.GetCellDimension() == 2:
                    int_frac_nodes_to_triangles[n].add(cell_index)

    # Knowing which triangle touches which internal fracture node is a temporary,
    # since we need to extract the point ids of the touching triangles.
    int_frac_nodes_to_triangle_nodes = defaultdict(list)
    for node, triangles in int_frac_nodes_to_triangles.items():
        for triangle in triangles:
            tr = grid.GetCell(triangle)
            # triangle_nodes = {tr.GetPointId(0), tr.GetPointId(1), tr.GetPointId(2)}
            triangle_nodes = set()
            for i in range(tr.GetNumberOfPoints()):
                triangle_nodes.add(tr.GetPointId(i))
            int_frac_nodes_to_triangle_nodes[node].append(triangle_nodes)

    # Now we must find the neighboring cells which do not cross the fracture (so we know they are on the "same" side).
    # To do this, we first find all the connections which cross.
    # Any other connection (i.e. non-crossing) will be a valid connection then.
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
                    graph.add_edge(ic0, ic1)
            else:
                print(f"No connection between {ic0} and {ic1}.")
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


def extract_fracture_nodes(grid: vtk.vtkUnstructuredGrid) -> Tuple[Iterable[int], Iterable[int]]:
    """
    From the vtk unstructured grid, extract the nodes that are part of the fracture.
    All the 2d cells are assumed to define the fracture. Thus all the nodes of thos 2d cells are taken.
    Then the nodes are split in 2. Nodes that are inside the fracture and those at the boundary of the fracture.
    :param grid: the vtk unstructured mesh
    :return: internal_fracture_nodes, external_fracture_nodes
    """
    # It's necessary to find which nodes need to be split.
    # The nodes on the boundary do not need to be split.
    # To find the boundary edges, we find those which are only referenced by one unique 2d element.
    fracture_edges = defaultdict(int)
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        if cell.GetCellDimension() == 2:
            for ie in range(cell.GetNumberOfEdges()):
                edge = cell.GetEdge(ie)
                edge_nodes = []
                for ien in range(edge.GetNumberOfPoints()):  # TODO, do less pedantic
                    edge_nodes.append(edge.GetPointId(ien))
                fracture_edges[tuple(sorted(edge_nodes))] += 1
    # Edges that are mentioned once are boundary edges.
    external_fracture_nodes = []
    for kv in filter(lambda fe: fe[1] == 1, fracture_edges.items()):
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


def clean_input_mesh(vtk_input_file: str) -> vtk.vtkUnstructuredGrid:
    """
    Takes the hi24L_coarse mesh and remove all the vertex, lines... cells that are of no use for us.
    Then returns a clean vtk instance containing 2d and 3d cells.
    :param vtk_input_file: The initial vtk file name.
    :return: The vtk unstructured grid instance.
    """
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_input_file)
    reader.Update()
    input_ = reader.GetOutput()
    # We want to keep track of the cells types and the business zone (defined using the `CellEntityIds` field).
    # vtk type -> zone -> cell id
    keep_cells = defaultdict(lambda: defaultdict(list))
    # Getting the physical field.
    zones = input_.GetCellData().GetArray('CellEntityIds')
    assert zones.GetNumberOfValues() == input_.GetNumberOfCells()

    for i in range(input_.GetNumberOfCells()):
        zone, cell_type = zones.GetValue(i), input_.GetCellType(i)
        keep_cells[cell_type][zone].append(i)

    keep_cells = dict.copy(keep_cells)  # const version
    # We only want to keep triangles and for faults and the tets.
    faults = {1, 4, 7, 10, 11, 12}
    num_tets = len(keep_cells[vtk.VTK_TETRA][1])  # zone `1` is the volumic domain.
    num_triangles = 0
    for fault in faults:
        num_triangles += len(keep_cells[vtk.VTK_TRIANGLE][fault])
    # The new total of cells will be num_tets + num_triangles

    output_cells = vtk.vtkCellArray()
    output_cells.AllocateExact(num_tets + num_triangles, num_tets * 4 + num_triangles * 3)

    for tet_id in keep_cells[vtk.VTK_TETRA][1]:
        tet = input_.GetCell(tet_id)
        output_cells.InsertNextCell(tet)
    attributes = [0] * num_tets  # The volumic domain will be `0` in the new mesh.
    output_cell_types = [vtk.VTK_TETRA] * num_tets

    for region, triangle_ids in keep_cells[vtk.VTK_TRIANGLE].items():
        if region in faults:
            for triangle_id in triangle_ids:
                triangle = input_.GetCell(triangle_id)
                output_cells.InsertNextCell(triangle)
            attributes += [region] * len(triangle_ids)  # keep same id for fractures, don't increment (there may be holes in the numbering).
    output_cell_types += [vtk.VTK_TRIANGLE] * num_triangles

    assert len(attributes) == num_tets + num_triangles
    assert len(output_cell_types) == num_tets + num_triangles
    assert output_cells.GetNumberOfCells() == num_tets + num_triangles

    # Let's loop over all the cells and insert them
    # In the meantime, we should be creating the attribute array.
    att = vtk.util.numpy_support.numpy_to_vtk(attributes)
    att.SetName("attribute")

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


def main():
    vtk_input_file = "HI24L_coars.vtk"
    output_mesh = clean_input_mesh(vtk_input_file)
    internal_fracture_nodes, external_fracture_nodes = extract_fracture_nodes(output_mesh)
    connected_components = color_fracture_sides(output_mesh, internal_fracture_nodes)

    duplicated = duplicate_fracture_nodes(output_mesh, connected_components, internal_fracture_nodes)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName("dupli.vtu")
    writer.SetInputData(duplicated)
    writer.Write()


if __name__ == '__main__':
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s', level=logging.INFO)
    args = parse(sys.argv[1:])
    logging.info("Playing with hi24L coarse.")
    main()

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
