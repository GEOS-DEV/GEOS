from collections import defaultdict
from dataclasses import dataclass
import logging
from typing import (
    Tuple,
    Iterable,
    Dict,
    Mapping,
    FrozenSet,
    List,
    Set,
    Sequence,
    Collection,
)

from tqdm import tqdm
import networkx
import numpy

from vtkmodules.vtkCommonCore import (
    vtkIdList,
    vtkIntArray,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    vtkCell,
    vtkCellArray,
    vtkPolygon,
    vtkUnstructuredGrid,
    VTK_POLYGON,
)
from vtkmodules.util.numpy_support import (
    vtk_to_numpy,
    numpy_to_vtk,
)
from vtkmodules.util.vtkConstants import VTK_ID_TYPE

from . import vtk_utils
from .vtk_utils import (
    vtk_iter,
    VtkOutput,
)


@dataclass(frozen=True)
class Options:
    policy: str
    field: str
    field_type: str
    field_values: FrozenSet[int]
    split_on_domain_boundary: bool
    vtk_output: VtkOutput
    vtk_fracture_output: VtkOutput


@dataclass(frozen=True)
class Result:
    info: str


@dataclass(frozen=True)
class FractureInfo:
    node_to_cells: Mapping[int, Iterable[int]] # For each _fracture_ node, gives all the cells that use this node.
    face_nodes: Iterable[Collection[int]]  # For each fracture face, returns the nodes of this face


def build_node_to_cells(mesh: vtkUnstructuredGrid,
                        face_nodes: Iterable[Iterable[int]]) -> Mapping[int, Iterable[int]]:
    node_to_cells: Dict[int, Set[int]] = defaultdict(set)  # TODO normally, just a list and not a set should be enough.

    fracture_nodes: Set[int] = set()
    for fns in face_nodes:
        for n in fns:
            fracture_nodes.add(n)

    for cell_id in tqdm(range(mesh.GetNumberOfCells()), desc="Computing the node to cells mapping"):
        cell_points: FrozenSet[int] = frozenset(vtk_iter(mesh.GetCell(cell_id).GetPointIds()))
        intersection: Iterable[int] = cell_points & fracture_nodes
        for node in intersection:
            node_to_cells[node].add(cell_id)

    return node_to_cells


def build_fracture_info(mesh: vtkUnstructuredGrid,
                        field: str,
                        field_values: FrozenSet[int]) -> FractureInfo:
    cell_data = mesh.GetCellData()
    if cell_data.HasArray(field):
        f = vtk_to_numpy(cell_data.GetArray(field))
    else:
        raise ValueError(f"Cell field {field} does not exist in mesh, nothing done")

    cells_to_faces: Dict[int, List[int]] = defaultdict(list)

    # For each face of each cell, we search for the unique neighbor cell (if it exists).
    # Then, if the 2 values of the two cells match the field requirements,
    # we store the cell and its local face index: this is indeed part of the surface that we'll need to be split.
    # neighbor_cell_ids = vtkIdList()
    for cell_id in tqdm(range(mesh.GetNumberOfCells()), desc="Computing the cell to faces mapping"):
        if f[cell_id] not in field_values:  # No need to consider a cell if its field value is not in the target range.
            continue
        cell = mesh.GetCell(cell_id)
        for i in range(cell.GetNumberOfFaces()):
            neighbor_cell_ids = vtkIdList()
            mesh.GetCellNeighbors(cell_id, cell.GetFace(i).GetPointIds(), neighbor_cell_ids)
            assert neighbor_cell_ids.GetNumberOfIds() < 2
            for j in range(neighbor_cell_ids.GetNumberOfIds()):    # It's 0 or 1...
                neighbor_cell_id = neighbor_cell_ids.GetId(j)
                if f[neighbor_cell_id] != f[cell_id] and f[neighbor_cell_id] in field_values:
                    cells_to_faces[cell_id].append(i)  # TODO add this (cell_is, face_id) information to the fracture_info?

    face_nodes: List[Collection[int]] = list()
    face_nodes_hashes: Set[FrozenSet[int]] = set()  # A temporary not to add multiple times the same face.
    for cell_id, faces_ids in tqdm(cells_to_faces.items(), desc="Extracting the faces of the fractures"):
        cell: vtkCell = mesh.GetCell(cell_id)
        for face_id in faces_ids:
            fn: Collection[int] = tuple(vtk_iter(cell.GetFace(face_id).GetPointIds()))
            fnh = frozenset(fn)
            if fnh not in face_nodes_hashes:
                face_nodes_hashes.add(fnh)
                face_nodes.append(fn)

    node_to_cells: Mapping[int, Iterable[int]] = build_node_to_cells(mesh, face_nodes)

    return FractureInfo(node_to_cells=node_to_cells, face_nodes=face_nodes)


def build_cell_to_cell_graph(mesh: vtkUnstructuredGrid,
                             fracture: FractureInfo) -> networkx.Graph:
    """
    Connects all the cells that touch the fracture by at least one node.
    Two cells are connected when they share at least a face which is not a face of the fracture.
    :param mesh: The input mesh.
    :param fracture: The fracture info.
    :return: The graph: each node of this graph is the index of the cell.
    There's an edge between two nodes of the graph if the cells share a face.
    """
    # Faces are identified by their nodes. But the order of those nodes may vary while referring to the same face.
    # Therefore we compute some kinds of hashes of those face to easily detect if a face is part of the fracture.
    tmp: List[FrozenSet[int]] = []
    for fn in fracture.face_nodes:
        tmp.append(frozenset(fn))
    face_hashes: FrozenSet[FrozenSet[int]] = frozenset(tmp)

    # We extract the list of the cells that touch the fracture by at least one node.
    cells: Set[int] = set()
    for cell_ids in fracture.node_to_cells.values():
        for cell_id in cell_ids:
            cells.add(cell_id)

    # Using the last precomputed containers, we're now building the dict which connects
    # every face (hash) of the fracture to the cells that touch the face...
    face_to_cells: Dict[FrozenSet[int], List[int]] = defaultdict(list)
    for cell_id in tqdm(cells, desc="Computing the cell to cell graph"):
        cell: vtkCell = mesh.GetCell(cell_id)
        for face_id in range(cell.GetNumberOfFaces()):
            face_hash: FrozenSet[int] = frozenset(vtk_iter(cell.GetFace(face_id).GetPointIds()))
            if face_hash not in face_hashes:
                face_to_cells[face_hash].append(cell_id)

    # ... eventually, when a face touches two cells, this means that those two cells share the same face
    # and should be connected in the final cell to cell graph.
    cell_to_cell = networkx.Graph()
    cell_to_cell.add_nodes_from(cells)
    cell_to_cell.add_edges_from(filter(lambda cs: len(cs) == 2, face_to_cells.values()))

    return cell_to_cell


def __identify_split(num_points: int,
                     cell_to_cell: networkx.Graph,
                     node_to_cells: Mapping[int, Iterable[int]]) -> Mapping[int, Mapping[int, int]]:
    """
    For each cell, compute the node indices replacements.
    :param num_points: Number of points in the whole mesh (not the fracture).
    :param cell_to_cell: The cell to cell graph (connection through common faces).
    :param node_to_cells: Maps the nodes of the fracture to the cells relying on this node.
    :return: For each cell (first key), returns a mapping from the current index
    and the new index that should replace the current index.
    Note that the current index and the new index can be identical: no replacement should be done then.
    """

    class NewIndex:
        """
        Returns the next available index.
        Note that the first time an index is met, the index itself is returned:
        we do not want to change an index if we do not have to.
        """
        def __init__(self, num_nodes: int):
            self.__current_last_index = num_nodes - 1
            self.__seen = set()

        def __call__(self, index) -> int:
            if index in self.__seen:
                self.__current_last_index += 1
                return self.__current_last_index
            else:
                self.__seen.add(index)
                return index

    build_new_index = NewIndex(num_points)
    result: Dict[int, Dict[int, int]] = defaultdict(dict)
    for node, cells in tqdm(sorted(node_to_cells.items()),  # Iteration over `sorted` nodes to have a predictable result for tests.
                            desc="Identifying the node duplications"):
        for connected_cells in networkx.connected_components(cell_to_cell.subgraph(cells)):
            # Each group of connect cells need around `node` must consider the same `node`.
            # Separate groups must have different (duplicated) nodes.
            new_index: int = build_new_index(node)
            for cell in connected_cells:
                result[cell][node] = new_index
    return result


def __perform_split(mesh: vtkUnstructuredGrid,
                    cell_to_node_mapping: Mapping[int, Mapping[int, int]]) -> vtkUnstructuredGrid:
    """
    Split the main 3d mesh based on the node duplication information contained in @p cell_to_node_mapping
    :param mesh: The main 3d mesh.
    :param cell_to_node_mapping: For each cell, gives the nodes that must be duplicated and their new index.
    :return: The main 3d mesh split at the fracture location.
    """
    added_points: Set[int] = set()
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            if i != o:
                added_points.add(o)
    num_new_points: int = mesh.GetNumberOfPoints() + len(added_points)

    # Creating the new points for the new mesh.
    old_points: vtkPoints = mesh.GetPoints()
    new_points = vtkPoints()
    new_points.SetNumberOfPoints(num_new_points)
    # Copying old points into the new container.
    for p in range(old_points.GetNumberOfPoints()):
        new_points.SetPoint(p, old_points.GetPoint(p))
    # Creating the new collocated/duplicated points based on the old points positions.
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            if i != o:
                new_points.SetPoint(o, old_points.GetPoint(i))

    # We are creating a new mesh.
    # The cells will be the same, except that their nodes may be duplicated or renumbered nodes.
    # The `new_cells` array will be modified in place.
    new_cells = vtkCellArray()
    new_cells.DeepCopy(mesh.GetCells())

    for c in tqdm(range(mesh.GetNumberOfCells()), desc="Performing the mesh split"):
        node_mapping: Mapping[int, int] = cell_to_node_mapping.get(c, {})
        # Extracting the point ids of the cell.
        # The values will be (potentially) overwritten in place, before being sent back into the cell.
        cell_point_ids = vtkIdList()
        new_cells.GetCellAtId(c, cell_point_ids)
        for i in range(cell_point_ids.GetNumberOfIds()):
            current_point_id: int = cell_point_ids.GetId(i)
            new_point_id: int = node_mapping.get(current_point_id, current_point_id)
            cell_point_ids.SetId(i, new_point_id)
        new_cells.ReplaceCellAtId(c, cell_point_ids)

    output = mesh.NewInstance()
    output.SetPoints(new_points)
    output.SetCells(mesh.GetCellTypesArray(), new_cells)  # The cell types are unchanged; we reuse the old cell types!
    return output


# def __generate_fracture_mesh(mesh: vtkUnstructuredGrid,
def __generate_fracture_mesh(mesh_points: vtkPoints,
                             fracture_info: FractureInfo,
                             cell_to_node_mapping: Mapping[int, Mapping[int, int]]) -> vtkUnstructuredGrid:
    """
    Generates the mesh of the fracture.
    :param mesh_points: The points of the main 3d mesh.
    :param fracture_info: The fracture description.
    :param cell_to_node_mapping: For each cell, gives the nodes that must be duplicated and their new index.
    :return: The fracture mesh.
    """
    logging.info("Generating the meshes")
    fracture_nodes: Collection[int] = tuple(fracture_info.node_to_cells.keys())
    num_points: int = len(fracture_nodes)

    points = vtkPoints()
    points.SetNumberOfPoints(num_points)
    tmp: Dict[int, int] = {}  # Building the node mapping, from 3d mesh nodes to 2d fracture nodes.
    for i, n in enumerate(fracture_nodes):
        coords: Tuple[float, float, float] = mesh_points.GetPoint(n)
        points.SetPoint(i, coords)
        tmp[n] = i

    polygons = vtkCellArray()
    for ns in fracture_info.face_nodes:
        polygon = vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(len(ns))
        for i, n in enumerate(ns):
            polygon.GetPointIds().SetId(i, tmp[n])
        polygons.InsertNextCell(polygon)

    # # TODO use an extended version of the `tmp` up there? Or vice versa?
    buckets: Dict[int, Set[int]] = defaultdict(set)
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            k: int = tmp[min(i, o)]
            buckets[k].update((i, o))

    assert set(buckets.keys()) == set(range(num_points))
    max_duplicated_nodes: int = max(map(len, buckets.values()))
    collocated_nodes = numpy.ones((num_points, max_duplicated_nodes), dtype=int) * -1
    for i, bucket in buckets.items():
        for j, val in enumerate(bucket):
            collocated_nodes[i, j] = val
    array = numpy_to_vtk(collocated_nodes, array_type=VTK_ID_TYPE)
    array.SetName("collocated_nodes")

    fracture_mesh = vtkUnstructuredGrid()  # We could be using vtkPolyData, but it's not supported by GEOSX for now.
    fracture_mesh.SetPoints(points)
    fracture_mesh.SetCells([VTK_POLYGON] * polygons.GetNumberOfCells(), polygons)
    fracture_mesh.GetPointData().AddArray(array)
    return fracture_mesh


def __split_mesh_on_fracture(mesh: vtkUnstructuredGrid,
                             options: Options) -> Tuple[vtkUnstructuredGrid, vtkUnstructuredGrid]:
    fracture: FractureInfo = build_fracture_info(mesh, options.field, options.field_values)
    cell_to_cell: networkx.Graph = build_cell_to_cell_graph(mesh, fracture)
    cell_to_node_mapping: Mapping[int, Mapping[int, int]] = __identify_split(mesh.GetNumberOfPoints(),
                                                                             cell_to_cell,
                                                                             fracture.node_to_cells)
    output_mesh: vtkUnstructuredGrid = __perform_split(mesh, cell_to_node_mapping)
    fractured_mesh: vtkUnstructuredGrid = __generate_fracture_mesh(mesh.GetPoints(), fracture, cell_to_node_mapping)
    return output_mesh, fractured_mesh


def __check(mesh, options: Options) -> Result:
    output_mesh, fracture_mesh = __split_mesh_on_fracture(mesh, options)
    vtk_utils.write_mesh(output_mesh, options.vtk_output)
    vtk_utils.write_mesh(fracture_mesh, options.vtk_fracture_output)
    # TODO provide statistics about what was actually performed (size of the fracture, number of split nodes...).
    return Result(info="OK")


def check(vtk_input_file: str, options: Options) -> Result:
    try:
        mesh = vtk_utils.read_mesh(vtk_input_file)
        return __check(mesh, options)
    except BaseException as e:
        logging.error(e)
        return Result(info="Something went wrong")
