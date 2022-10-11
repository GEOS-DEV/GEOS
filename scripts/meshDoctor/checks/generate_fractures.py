from collections import defaultdict
from dataclasses import dataclass
import logging
from typing import (
    Tuple,
    Iterable,
    Dict,
    FrozenSet,
    Set,
    Iterator,
    Sequence,
    Collection,
)

import numpy

from vtkmodules.vtkCommonCore import (
    vtkIdList,
    vtkIntArray,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkUnstructuredGrid,
)
from vtkmodules.vtkFiltersGeometry import (
    vtkMarkBoundaryFilter,
)
from vtk.util.numpy_support import (
    vtk_to_numpy,
)

import networkx

from . import vtk_utils


@dataclass(frozen=True)
class Options:
    policy: str
    field: str
    field_values: FrozenSet[int]
    split_on_domain_boundary: bool
    output: str


@dataclass(frozen=True)
class Result:
    info: str


@dataclass(frozen=True)
class FractureNodesInfo:
    is_internal_fracture_node: Sequence[bool]
    is_boundary_fracture_node: Sequence[bool]


class FractureCellsInfo:
    def __init__(self, mesh: vtkUnstructuredGrid, cell_to_faces: Dict[int, Iterable[int]]):
        # For each cell (the key), gets the local (to the face, i.e. 0 to 3 for a tet) faces that are part of the fracture.
        # There can be multiple faces involved.
        self.cell_to_faces: Dict[int, Iterable[int]] = cell_to_faces
        tmp = defaultdict(list)
        for cell_id, faces_id in cell_to_faces.items():
            cell = mesh.GetCell(cell_id)
            for face_id in faces_id:
                face = cell.GetFace(face_id)
                face_point_ids = frozenset(_iter(face.GetPointIds()))
                tmp[face_point_ids] += cell_id, face_id
        # This field is really dedicated to the writing into the vtk field data.
        self.field_data: Collection[Tuple[int, int, int, int]] = tuple(tmp.values())
        for data in self.field_data:
            assert len(data) == 4


@dataclass(frozen=True)
class DuplicatedNodesInfo:
    # Spans the old number of nodes.
    # Each value is either the new number of the node or `-1` if it's a duplicated node.
    renumbered_nodes: Collection[int]
    # For each original node that has been duplicated (or more) as the key,
    # gets all the new numberings of the nodes.
    duplicated_nodes: Dict[int, Collection[int]]


def _iter(id_list: vtkIdList) -> Iterator[int]:
    """
    Utility function transforming a vtkIdList into an iterable to be used for building built-ins python containers.
    :param id_list: the vtkIdList.
    :return: The iterator.
    """
    for i in range(id_list.GetNumberOfIds()):
        yield id_list.GetId(i)


def __duplicate_fracture_nodes(mesh: vtkUnstructuredGrid,
                               connected_components: Iterable[Iterable[int]],
                               node_frac_info: FractureNodesInfo) -> Tuple[vtkUnstructuredGrid, DuplicatedNodesInfo]:
    """
    Splits the mesh on the fracture.
    New nodes are then created and inserted in the new mesh.
    All the already existing cells are updated to take into account these new nodes.
    The 2d elements of the fracture point to duplicated nodes,
    but there is no specific action to be sure that their nodes belong to the same 3d element.
    Node positions will be correct anyway.
    :param mesh: the vtk unstructured mesh
    :param connected_components:
        Each component contains an iterable of cells that touch the internal fracture nodes.
        All the cells in one same bucket will be based on the same duplicated node.
    :param internal_fracture_nodes: Iterable of all the internal fracture nodes.
    :return:
        A vtk unstructured mesh with duplicated nodes.
        The support of the cells take the duplicated nodes into account.
    """
    num_points: int = mesh.GetNumberOfPoints()
    num_cells: int = mesh.GetNumberOfCells()

    # We are creating a new mesh.
    # The cells will be the same, except that their nodes may be duplicated or renumbered nodes.
    # The `new_cells` array will be modified in place, while the original grid remains unmodified.
    new_cells = vtkCellArray()
    new_cells.DeepCopy(mesh.GetCells())

    # Building an indicator to find fracture cells.
    # A fracture cell is a cell which touches a fracture internal node.
    is_fracture_side_cells = numpy.zeros(num_cells, dtype=bool)  # Defaults to "not a fracture cell".
    for fracture_side_cells in connected_components:
        for fracture_side_cell in fracture_side_cells:
            is_fracture_side_cells[fracture_side_cell] = True

    # Each non duplicated nodes is renumbered.
    # The idea is that I do not want to fill holes in the numbering,
    # I just want to increment the next number by one.
    # A negative number means that this is a duplicated number.
    renumbered_nodes = numpy.ones(num_points, dtype=int) * -1  # Defaults to "moved".
    next_point_id = 0  # When a node needs to be duplicated, `next_point_id` is the next available node.
    for node_idx in range(num_points):
        if not node_frac_info.is_internal_fracture_node[node_idx]:
            renumbered_nodes[node_idx] = next_point_id
            next_point_id += 1

    # Connected component index -> duplicated nodes.
    # It's not super useful for the moment to store the precise binding with the component.
    component_to_dup_nodes: Dict[int, Dict[int, int]] = dict()

    # First, dealing with the renumbering of the cell touching the fracture.
    for component_idx, fracture_side_cells in enumerate(connected_components):
        dup_nodes: Dict[int, int] = dict()
        for fracture_side_cell in fracture_side_cells:
            # `new_cell_point_ids` will contain the ids for the new cells.
            # Calling `new_cells.GetCellAtId` is a trick to get it at the correct size,
            # because all values should be overwritten.
            new_cell_point_ids = vtkIdList()
            new_cells.GetCellAtId(fracture_side_cell, new_cell_point_ids)
            for i in range(new_cell_point_ids.GetNumberOfIds()):
                new_cell_point_id = new_cell_point_ids.GetId(i)
                # Here we find or build the new number of the node.
                if node_frac_info.is_internal_fracture_node[new_cell_point_id]:
                    assert renumbered_nodes[new_cell_point_id] == -1
                    if new_cell_point_id in dup_nodes:  # if already duplicated, take the same
                        new_cell_point_ids.SetId(i, dup_nodes[new_cell_point_id])
                    else:  # otherwise, duplicate
                        dup_nodes[new_cell_point_id] = next_point_id
                        new_cell_point_ids.SetId(i, next_point_id)
                        next_point_id += 1
                else:  # Here, it's not an internal fracture node which was not duplicated.
                    assert not renumbered_nodes[new_cell_point_id] == -1
                    new_cell_point_ids.SetId(i, renumbered_nodes[new_cell_point_id])
            new_cells.ReplaceCellAtId(fracture_side_cell, new_cell_point_ids)
        component_to_dup_nodes[component_idx] = dup_nodes

    # Now, we must not forget to modify the nodes of all the other cells,
    # because the nodes of all the points have been modified!
    # This is for non fracture 3d elements!
    for cell_idx in filter(lambda c: not is_fracture_side_cells[c], range(num_cells)):
        # Same trick to get the indices list at the correct size.
        new_cell_point_ids = vtkIdList()
        new_cells.GetCellAtId(cell_idx, new_cell_point_ids)
        # TODO Is there something to be done for 2D/3D cells?
        for i in range(new_cell_point_ids.GetNumberOfIds()):
            new_cell_point_id = new_cell_point_ids.GetId(i)
            assert not renumbered_nodes[new_cell_point_id] == -1
            new_cell_point_ids.SetId(i, renumbered_nodes[new_cell_point_id])
        new_cells.ReplaceCellAtId(cell_idx, new_cell_point_ids)

    # Now we finish the process for vtk points.
    num_new_points = num_points - sum(node_frac_info.is_internal_fracture_node) + sum(map(len, component_to_dup_nodes.values()))
    assert next_point_id == num_new_points

    old_points = mesh.GetPoints()
    new_points = vtkPoints()
    new_points.SetNumberOfPoints(num_new_points)

    for from_, to in enumerate(renumbered_nodes):
        if not to == -1:
            new_points.SetPoint(to, old_points.GetPoint(from_))

    duplicated_nodes: Dict[int, Set[int]] = defaultdict(set)
    for dup_nodes in component_to_dup_nodes.values():
        for from_, to in dup_nodes.items():
            duplicated_nodes[from_].add(to)

    for from_, tos in duplicated_nodes.items():
        for to in tos:
            new_points.SetPoint(to, old_points.GetPoint(from_))

    # Validation
    nodes_validation = numpy.zeros(num_new_points, dtype=bool)
    for n in range(num_points):
        if renumbered_nodes[n] == -1:
            assert n in duplicated_nodes
            for to in duplicated_nodes[n]:
                nodes_validation[to] = True
        else:
            assert n not in duplicated_nodes
            nodes_validation[renumbered_nodes[n]] = True
    assert all(nodes_validation)

    output = mesh.NewInstance()  # keeping the same instance type.
    output.SetPoints(new_points)
    output.SetCells(mesh.GetCellTypesArray(), new_cells)  # The cell types are unchanged; we reuse the old cell types!

    duplicated_nodes_info = DuplicatedNodesInfo(renumbered_nodes=renumbered_nodes,
                                                duplicated_nodes=duplicated_nodes)

    return output, duplicated_nodes_info


def __copy_fields(input_mesh: vtkUnstructuredGrid,
                  output_mesh: vtkUnstructuredGrid,
                  cell_frac_info: FractureCellsInfo,
                  duplicated_nodes_info: DuplicatedNodesInfo) -> None:
    """
    Given input and output meshes, copies the fields from the input into the output.
    Special care is taken because of the nodes renumbering.
    In case nodes are duplicated, then the value of the point field at this node is duplicated too.
    :param input_mesh: The input mesh from which we'll take the fields.
    :param output_mesh: The output mesh that will receive the fields. It should already have its points and cells defined consistently.
    :param cell_frac_info: This information is written into the field data of the vtk mesh.
    :param duplicated_nodes_info: This information is written into the field data of the vtk mesh and is used when copying the points fields.
    :return: None
    """
    # Copying the cell data.
    # The cells are the same, just their nodes support have changed.
    input_cell_data = input_mesh.GetCellData()
    for i in range(input_cell_data.GetNumberOfArrays()):
        input_array = input_cell_data.GetArray(i)
        logging.debug(f"Copying cell field {input_array.GetName()}")
        output_mesh.GetCellData().AddArray(input_array)
    # Copying the point data.
    input_point_data = input_mesh.GetPointData()
    for i in range(input_point_data.GetNumberOfArrays()):
        input_array = input_point_data.GetArray(i)
        logging.debug(f"Copying point field {input_array.GetName()}")
        tmp = input_array.NewInstance()
        tmp.SetName(input_array.GetName())
        tmp.SetNumberOfComponents(input_array.GetNumberOfComponents())
        tmp.SetNumberOfTuples(output_mesh.GetNumberOfPoints())
        # Now copy data, taking into account the duplicated nodes.
        for from_, to in duplicated_nodes_info.renumbered_nodes:
            if to != -1:
                tmp.SetTuple(to, from_, input_array)
        for from_, tos in duplicated_nodes_info.duplicated_nodes:
            # We duplicated the point information for duplicated nodes.
            for to in tos:
                tmp.SetTuple(to, from_, input_array)
        output_mesh.GetPointData().AddArray(tmp)

    # The "field data" will contain the fracture information
    field_data = input_mesh.GetFieldData()
    if field_data.GetNumberOfArrays() > 0:
        logging.warning(f"Copying field data that already has arrays. No modification was made w.r.t. nodes duplications.")

    # Building the field data for the fracture
    if field_data.HasArray("fracture_info"):
        logging.warning("Field data \"fracture_info\" already exists, nothing done.")
    else:
        frac_array = vtkIntArray()
        frac_array.SetName("fracture_info")  # TODO hard coded name.
        frac_array.SetNumberOfComponents(4)  # Warning the component has to be defined first...
        frac_array.SetNumberOfTuples(len(cell_frac_info.field_data))
        for i, data in enumerate(cell_frac_info.field_data):
            frac_array.SetTuple(i, data)
        field_data.AddArray(frac_array)

    # The nodes bindings field data
    if field_data.HasArray("duplicated_points_info"):
        logging.warning("Field data \"duplicated_points_info\" already exists, nothing done.")
    else:
        duplication_max_multiplicity = max(map(len, duplicated_nodes_info.duplicated_nodes.values()))
        nodes_array = vtkIntArray()
        nodes_array.SetName("duplicated_points_info")
        nodes_array.SetNumberOfComponents(duplication_max_multiplicity)  # Warning the component has to be defined first...
        nodes_array.SetNumberOfTuples(len(duplicated_nodes_info.duplicated_nodes))
        for i, data in enumerate(duplicated_nodes_info.duplicated_nodes.values()):
            tmp = list(data) + [-1] * (duplication_max_multiplicity - len(data))
            nodes_array.SetTuple(i, tmp)
        field_data.AddArray(nodes_array)

    output_mesh.SetFieldData(field_data)


def __color_fracture_sides(mesh: vtkUnstructuredGrid, cell_frac_info: FractureCellsInfo, node_frac_info: FractureNodesInfo) -> Iterable[Iterable[int]]:
    """
    Given all the cells that are in contact with the detected fracture,
    we separate them into bucket of connected cells touching the fractures.
    We do this because all the cells in one same bucket are connected and need to stay connected.
    So they need to share the same nodes as they did before the split.
    :return: All the buckets connected fracture cells.
    """
    def does_face_contain_boundary_node(_face_point_ids: Iterable[int]) -> bool:  # Small helper function.
        for face_point_id in _face_point_ids:
            if node_frac_info.is_boundary_fracture_node[face_point_id]:
                return True
        return False

    face_node_set_to_cell = defaultdict(list)  # For each face (defined by its node set), gives all the cells containing this face.
    for c, local_frac_f in cell_frac_info.cell_to_faces.items():
        cell = mesh.GetCell(c)
        for f in [i for i in range(cell.GetNumberOfFaces()) if i not in local_frac_f]:
            face = cell.GetFace(f)
            face_point_ids = frozenset(_iter(face.GetPointIds()))
            if not does_face_contain_boundary_node(face_point_ids):
                face_node_set_to_cell[face_point_ids].append(c)

    graph = networkx.Graph()
    # Let's add all the nodes first: all the cells are part of the graph.
    # It's important for isolated cells that have no connection with other cells.
    # If we rely on the edge creation to insert the nodes, then those cells would be forgotten.
    graph.add_nodes_from(cell_frac_info.cell_to_faces.keys())
    for cells in face_node_set_to_cell.values():
        assert 0 < len(cells) < 3
        if len(cells) == 2:
            graph.add_edge(cells[0], cells[1])

    return tuple(networkx.connected_components(graph))


def __find_boundary_nodes(mesh: vtkUnstructuredGrid) -> Sequence[int]:
    """
    Find all the points on the boundary of the mesh.
    :param mesh: the vtk unstructured mesh
    :return: A boundary point indicator array (value is 1 for a boundary point, 0 otherwise)
    """
    f = vtkMarkBoundaryFilter()
    f.SetBoundaryPointsName("boundary points")
    f.SetBoundaryCellsName("boundary cells")

    f.SetInputData(mesh)
    f.Update()

    output_mesh = f.GetOutputDataObject(0)
    point_data = output_mesh.GetPointData()
    assert point_data.GetNumberOfArrays() == 1
    assert point_data.GetArrayName(0) == "boundary points"
    is_boundary_point = point_data.GetArray(0)
    return vtk_to_numpy(is_boundary_point)


def __build_fracture_nodes(mesh: vtkUnstructuredGrid,
                           cell_frac_info: FractureCellsInfo,
                           split_on_domain_boundary: bool) -> FractureNodesInfo:
    """
    Given the description of the fracture in @p cell_frac_info, computes the underlying nodes.
    :param mesh: the vtk unstructured mesh
    :param cell_frac_info: The geometrical mesh description of the fracture.
    :param split_on_domain_boundary: If a fracture touches the boundary of the domain, what should we do?
    :return: The fracture nodes information.
    """
    fracture_edges: Dict[Tuple[int, int], int] = defaultdict(int)
    for c, fs in cell_frac_info.cell_to_faces.items():
        cell = mesh.GetCell(c)
        for f in fs:
            face = cell.GetFace(f)
            for e in range(face.GetNumberOfEdges()):
                edge = face.GetEdge(e)
                edge_nodes = []
                for i in range(edge.GetNumberOfPoints()):  # TODO, do less pedantic
                    edge_nodes.append(edge.GetPointId(i))
                fracture_edges[tuple(sorted(edge_nodes))] += 1  # TODO frozenset?

    boundary_fracture_edges = []
    # Boundary edges are seen twice because each 2d fracture element is seen twice too.
    # (Each side of the fracture).
    for kv in filter(lambda fe: fe[1] == 2, fracture_edges.items()):
        boundary_fracture_edges.append(kv[0])

    is_fracture_node = numpy.zeros(mesh.GetNumberOfPoints(), dtype=bool)
    for nodes in fracture_edges.keys():
        for n in nodes:
            is_fracture_node[n] = True

    is_boundary_fracture_node = numpy.zeros(mesh.GetNumberOfPoints(), dtype=bool)
    for edge in boundary_fracture_edges:
        for n in edge:
            is_boundary_fracture_node[n] = True

    if split_on_domain_boundary:
        # Here, we split on domain boundaries
        is_boundary_point = __find_boundary_nodes(mesh)
        # Here, `is_boundary_fracture_node` is filled properly.
        # In this `if` branch, we want all the nodes that only belong to boundary edges
        # to be in fact part of `is_internal_fracture_node`.
        # But for the nodes of the edges that have one node on the mesh boundary
        # and one node inside the domain should eventually be `internal fracture nodes`
        # and not `boundary fracture nodes`.
        move_to_boundary_fracture_node = numpy.zeros(mesh.GetNumberOfPoints(), dtype=bool)
        for n0, n1 in boundary_fracture_edges:
            if not (is_boundary_point[n0] and is_boundary_point[n1]):
                move_to_boundary_fracture_node[n0] = True
                move_to_boundary_fracture_node[n1] = True
        for n, move in enumerate(move_to_boundary_fracture_node):
            is_boundary_fracture_node[n] = move

    # Now compute the internal fracture nodes by "difference".
    is_internal_fracture_node = numpy.zeros(mesh.GetNumberOfPoints(), dtype=bool)
    for n in range(mesh.GetNumberOfPoints()):  # TODO duplicated, reorg the code.
        if is_fracture_node[n] and not is_boundary_fracture_node[n]:
            is_internal_fracture_node[n] = True

    return FractureNodesInfo(is_internal_fracture_node=is_internal_fracture_node,
                             is_boundary_fracture_node=is_boundary_fracture_node)


def __find_involved_cells(mesh: vtkUnstructuredGrid, options: Options) -> Tuple[FractureCellsInfo, FractureNodesInfo]:
    """
    To do the split, we need to find
        - all the cells that are touching the fracture,
        - the faces of those cells that are touching the fracture. TODO we may have no face (just an edge, for boundary purposes)
    For convenience, will also return all the nodes belonging to those faces. (Let's call them fracture nodes.)
    Also, there's something still not done about external fracture nodes and internal fracture nodes...
    But a choice needs to be made when we are at the boundary of the simulation domain.
    """
    cell_data = mesh.GetCellData()
    if cell_data.HasArray(options.field):
        f = vtk_to_numpy(cell_data.GetArray(options.field))
    else:
        raise ValueError(f"Cell field {options.field} does not exist in mesh, nothing done")

    cells_to_faces = defaultdict(list)
    is_fracture_node = numpy.zeros(mesh.GetNumberOfPoints(), dtype=bool)  # all False

    # For each face of each cell, we search for the unique neighbor cell (if it exists).
    # Then, if the 2 values of the two cells match the field requirements,
    # we store the cell and its local face index: this is indeed part of the surface that we'll need to be split.
    neighbor_cell_ids = vtkIdList()
    for cell_id in range(mesh.GetNumberOfCells()):
        if f[cell_id] not in options.field_values:  # No need to consider a cell if its field value is not in the target range.
            continue
        cell = mesh.GetCell(cell_id)
        for i in range(cell.GetNumberOfFaces()):
            face = cell.GetFace(i)
            mesh.GetCellNeighbors(cell_id, face.GetPointIds(), neighbor_cell_ids)
            assert neighbor_cell_ids.GetNumberOfIds() < 2
            for j in range(neighbor_cell_ids.GetNumberOfIds()):  # It's 0 or 1...
                neighbor_cell_id = neighbor_cell_ids.GetId(j)
                if f[neighbor_cell_id] != f[cell_id] and f[neighbor_cell_id] in options.field_values:
                    cells_to_faces[cell_id].append(i)
                    for k in range(face.GetNumberOfPoints()):
                        is_fracture_node[face.GetPointId(k)] = True

    logging.warning("The fracture split feature is still work in progress.")
    cell_frac_info: FractureCellsInfo = FractureCellsInfo(mesh, cells_to_faces)
    node_frac_info: FractureNodesInfo = __build_fracture_nodes(mesh, cell_frac_info, options.split_on_domain_boundary)
    # TODO what should happen if the fracture face of a cell is not to be split at all?

    logging.warning(cell_frac_info.cell_to_faces)
    logging.warning(node_frac_info)
    return cell_frac_info, node_frac_info


def __split_mesh_on_fracture(mesh, options: Options) -> vtkUnstructuredGrid:
    cell_frac_info, node_frac_info = __find_involved_cells(mesh, options)
    connected_cells = __color_fracture_sides(mesh, cell_frac_info, node_frac_info)
    output_mesh, duplicated_nodes_info = __duplicate_fracture_nodes(mesh, connected_cells, node_frac_info)
    __copy_fields(mesh, output_mesh, cell_frac_info, duplicated_nodes_info)
    return output_mesh


def __check(mesh, options: Options) -> Result:
    output_mesh = __split_mesh_on_fracture(mesh, options)
    if options.output:
        vtk_utils.write_mesh(output_mesh, options.output)
    # TODO provide statistics about what was actually performed (size of the fracture, number of split nodes...).
    return Result(info="OK")


def check(vtk_input_file: str, options: Options) -> Result:
    try:
        mesh = vtk_utils.read_mesh(vtk_input_file)
        return __check(mesh, options)
    except BaseException as e:
        logging.error(e)
        return Result(info="Something went wrong")
