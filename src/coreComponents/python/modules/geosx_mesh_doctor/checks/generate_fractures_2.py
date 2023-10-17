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
    face_nodes: Iterable[Iterable[int]]  # For each fracture face, returns the nodes of this face


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
    for cell_id in range(mesh.GetNumberOfCells()):
        if f[cell_id] not in field_values:  # No need to consider a cell if its field value is not in the target range.
            continue
        cell = mesh.GetCell(cell_id)
        for i in range(cell.GetNumberOfFaces()):
            neighbor_cell_ids = vtkIdList()
            mesh.GetCellNeighbors(cell_id, cell.GetFace(i).GetPointIds(), neighbor_cell_ids)
            # face = cell.GetFace(i)
            # mesh.GetCellNeighbors(cell_id, face.GetPointIds(), neighbor_cell_ids)
            assert neighbor_cell_ids.GetNumberOfIds() < 2
            for j in range(neighbor_cell_ids.GetNumberOfIds()):    # It's 0 or 1...
                neighbor_cell_id = neighbor_cell_ids.GetId(j)
                if f[neighbor_cell_id] != f[cell_id] and f[neighbor_cell_id] in field_values:
                    cells_to_faces[cell_id].append(i)  # TODO add this (cell_is, face_id) information to the fracture_info?

    face_nodes: List[Iterable[int]] = list()
    face_nodes_hashes: Set[FrozenSet[int]] = set()
    for cell_id, faces_ids in cells_to_faces.items():
        cell: vtkCell = mesh.GetCell(cell_id)
        for face_id in faces_ids:
            fn: Tuple[int, ...] = tuple(vtk_iter(cell.GetFace(face_id).GetPointIds()))
            # face: vtkCell = cell.GetFace(face_id)
            # fn: Tuple[int, ...] = tuple(vtk_iter(face.GetPointIds()))
            fnh = frozenset(fn)
            if fnh not in face_nodes_hashes:
                face_nodes_hashes.add(fnh)
                face_nodes.append(fn)

    # computing the node_to_cells mapping
    node_to_cells: Dict[int, Set[int]] = defaultdict(set)  # TODO normally, just a list and not a set should be enough.

    fracture_nodes: Set[int] = set()
    for fns in face_nodes:
        for n in fns:
            fracture_nodes.add(n)

    for cell_id in range(mesh.GetNumberOfCells()):
        # cell: vtkCell = mesh.GetCell(cell_id)
        # cell_points: Set[int] = set(vtk_iter(cell.GetPointIds()))
        cell_points: FrozenSet[int] = frozenset(vtk_iter(mesh.GetCell(cell_id).GetPointIds()))
        intersection: Iterable[int] = cell_points & fracture_nodes
        for node in intersection:
            node_to_cells[node].add(cell_id)

    return FractureInfo(node_to_cells=node_to_cells, face_nodes=face_nodes)


def build_cell_to_cell_graph(mesh: vtkUnstructuredGrid,
                             fracture: FractureInfo) -> networkx.Graph:
    # TODO connect through boundary edges if needed
    tmp: List[FrozenSet[int]] = []
    for fn in fracture.face_nodes:
        tmp.append(frozenset(fn))
    face_hashes: FrozenSet[FrozenSet[int]] = frozenset(tmp)

    cells: Set[int] = set()
    for cell_ids in fracture.node_to_cells.values():
        for cell_id in cell_ids:
            cells.add(cell_id)

    face_to_cells: Dict[FrozenSet[int], List[int]] = defaultdict(list)
    for cell_id in cells:
        cell: vtkCell = mesh.GetCell(cell_id)
        for face_id in range(cell.GetNumberOfFaces()):
            face_hash: FrozenSet[int] = frozenset(vtk_iter(cell.GetFace(face_id).GetPointIds()))
            if face_hash not in face_hashes:
                face_to_cells[face_hash].append(cell_id)

    cell_to_cell = networkx.Graph()
    cell_to_cell.add_nodes_from(cells)
    cell_to_cell.add_edges_from(filter(lambda cs: len(cs) == 2, face_to_cells.values()))

    return cell_to_cell


class NewIndex:
    def __init__(self, num_nodes: int):
        self.__next_index = num_nodes
        self.__seen = set()

    def __call__(self, index) -> int:
        if index not in self.__seen:
            self.__seen.add(index)
            return index
        else:
            result = self.__next_index
            self.__next_index += 1
            return result


def __identify_split(num_points: int,
                     cell_to_cell: networkx.Graph,
                     node_to_cells: Mapping[int, Iterable[int]]) -> Mapping[int, Mapping[int, int]]:
    """
    :param mesh:
    :param cell_to_cell:
    :param node_to_cells:
    :return: For each cell, the node mapping to perform when splitting.
    """
    build_new_index = NewIndex(num_points)
    result: Dict[int, Dict[int, int]] = defaultdict(dict)
    for node, cells in sorted(node_to_cells.items()):  # TODO can I iterate of _sorted_ nodes?
        for connected_cells in networkx.connected_components(cell_to_cell.subgraph(cells)):
            new_index: int = build_new_index(node)
            for cell in connected_cells:  # Same new index for all the cells in the same connected component.
                result[cell][node] = new_index
    return result


def __perform_split(mesh: vtkUnstructuredGrid,
                    cell_to_node_mapping: Mapping[int, Mapping[int, int]]) -> vtkUnstructuredGrid:
    added_points: Set[int] = set()
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            if i != o:
                added_points.add(o)
    num_new_points: int = mesh.GetNumberOfPoints() + len(added_points)

    old_points: vtkPoints = mesh.GetPoints()
    new_points = vtkPoints()
    new_points.SetNumberOfPoints(num_new_points)
    # Copying old points into the new container.
    for p in range(old_points.GetNumberOfPoints()):
        new_points.SetPoint(p, old_points.GetPoint(p))
    # Creating the new collocated/duplicated points based on the old points.
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            if i != o:
                new_points.SetPoint(o, old_points.GetPoint(i))

    # We are creating a new mesh.
    # The cells will be the same, except that their nodes may be duplicated or renumbered nodes.
    # The `new_cells` array will be modified in place, while the original grid remains unmodified.
    new_cells = vtkCellArray()
    new_cells.DeepCopy(mesh.GetCells())

    for c in range(mesh.GetNumberOfCells()):
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
    # for c in range(mesh.GetNumberOfCells()):
    #     node_mapping: Dict[int, int] = cell_node_split.get(c)
    #     if not node_mapping:
    #         continue
    #     new_cell_point_ids = vtkIdList()
    #     new_cells.GetCellAtId(c, new_cell_point_ids)
    #     for i in range(new_cell_point_ids.GetNumberOfIds()):
    #         new_cell_point_id = new_cell_point_ids.GetId(i)
    #         new_i: int = node_mapping.get(new_cell_point_id)
    #         if new_i is not None:
    #             new_cell_point_ids.SetId(i, new_i)
    #     new_cells.ReplaceCellAtId(c, new_cell_point_ids)

    output = mesh.NewInstance()  # keeping the same instance type.
    output.SetPoints(new_points)
    output.SetCells(mesh.GetCellTypesArray(), new_cells)  # The cell types are unchanged; we reuse the old cell types!
    return output


def __generate_fracture_mesh(mesh: vtkUnstructuredGrid,
                             fracture_info: FractureInfo,
                             cell_to_node_mapping: Mapping[int, Mapping[int, int]]) -> vtkUnstructuredGrid:
    fracture_nodes: Collection[int] = tuple(fracture_info.node_to_cells.keys())
    num_points: int = len(fracture_nodes)

    points = vtkPoints()
    points.SetNumberOfPoints(num_points)
    tmp: Dict[int, int] = {}  # Building the node mapping, from 3d mesh nodes to 2d fracture nodes.
    for i, n in enumerate(fracture_nodes):
        coords: Tuple[float, float, float] = mesh.GetPoint(n)
        points.SetPoint(i, coords)
        tmp[n] = i

    polygons = vtkCellArray()
    for ns in fracture_info.face_nodes:
        polygon = vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(len(ns))
        for i, n in enumerate(ns):
            polygon.GetPointIds().SetId(i, tmp[n])
        polygons.InsertNextCell(polygon)

    buckets: Dict[int, Set[int]] = defaultdict(set)
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            k: int = tmp[min(i, o)]
            buckets[k].update((i, o))
    # # TODO use an extended version of the `tmp` up there? Or vice versa?
    # buckets: Dict[int, List[int]] = defaultdict(list)  # Bucket of collocated nodes (maps 2d fracture nodes to the 3d mesh nodes).
    # for node_mapping in cell_to_node_mapping.values():
    #     for i, o in node_mapping.items():
    #         if i != o:
    #             buckets[tmp[i]].append(o)

    assert set(buckets.keys()) == set(range(num_points))
    max_duplicated_nodes: int = max(map(len, buckets.values()))
    collocated_nodes = numpy.ones((num_points, max_duplicated_nodes), dtype=int) * -1
    for i, bucket in buckets.items():
        for j, val in enumerate(bucket):
            collocated_nodes[i, j] = val
    # for i, n in enumerate(fracture_nodes):
    #     for j, m in enumerate(buckets.get(n, ())):  # TODO confirm usage of `.get(n, ())`?
    #         collocated_nodes[i, j] = m
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
    cell_to_node_mapping: Mapping[int, Mapping[int, int]] = __identify_split(mesh.GetNumberOfPoints(), cell_to_cell, fracture.node_to_cells)
    output_mesh: vtkUnstructuredGrid = __perform_split(mesh, cell_to_node_mapping)
    return output_mesh, __generate_fracture_mesh(mesh, fracture, cell_to_node_mapping)


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
