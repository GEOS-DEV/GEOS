from collections import defaultdict
from dataclasses import dataclass
import logging
from typing import (
    Collection,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Mapping,
    Optional,
    Set,
    Sequence,
    Tuple,
)

from tqdm import tqdm
import networkx
import numpy

from vtkmodules.vtkCommonCore import (
    vtkIdList,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    vtkCell,
    vtkCellArray,
    vtkPolygon,
    vtkUnstructuredGrid,
    VTK_POLYGON,
    VTK_POLYHEDRON,
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
    to_vtk_id_list,
)
from .vtk_polyhedron import (
    FaceStream,
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
                        options: Options) -> FractureInfo:
    field = options.field
    field_values = options.field_values
    if options.policy == "field":
        cell_data = mesh.GetCellData()
        if cell_data.HasArray(field):
            f = vtk_to_numpy(cell_data.GetArray(field))
        else:
            raise ValueError(f"Cell field {field} does not exist in mesh, nothing done")
        cells_to_faces: Dict[int, List[int]] = defaultdict(list)
        # For each face of each cell, we search for the unique neighbor cell (if it exists).
        # Then, if the 2 values of the two cells match the field requirements,
        # we store the cell and its local face index: this is indeed part of the surface that we'll need to be split.
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
    
    if options.policy == "internal_surfaces":
        node_to_cells = {}
        face_nodes = []
        attribute = mesh.GetCellData().GetArray(field)
        for c in range( mesh.GetNumberOfCells() ):
            cell = mesh.GetCell(c)
            if cell.GetCellDimension() == 2:
                if attribute.GetValue(c) in field_values:
                    nodes = []
                    for v in range(cell.GetNumberOfPoints()):
                        node_to_cells[cell.GetPointId(v)] = []
                        nodes.append(cell.GetPointId(v))
                    face_nodes.append(tuple(nodes))

        for c in range( mesh.GetNumberOfCells() ):
            cell = mesh.GetCell(c)
            if cell.GetCellDimension() == 3:
                for v in range(cell.GetNumberOfPoints()):
                    if cell.GetPointId(v) in node_to_cells:
                        node_to_cells[cell.GetPointId(v)].append(c)


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
            self.__seen: Set[int] = set()

        def __call__(self, index: int) -> int:
            if index in self.__seen:
                self.__current_last_index += 1
                return self.__current_last_index
            else:
                self.__seen.add(index)
                return index

    build_new_index = NewIndex(num_points)
    result: Dict[int, Dict[int, int]] = defaultdict(dict)
    for node, cells in tqdm(sorted(node_to_cells.items()),  # Iteration over `sorted` nodes to have a predictable result for tests.
                            desc="Identifying the node splits"):
        for connected_cells in networkx.connected_components(cell_to_cell.subgraph(cells)):
            # Each group of connect cells need around `node` must consider the same `node`.
            # Separate groups must have different (duplicated) nodes.
            new_index: int = build_new_index(node)
            for cell in connected_cells:
                result[cell][node] = new_index
    return result


def __copy_fields(old_mesh: vtkUnstructuredGrid,
                  new_mesh: vtkUnstructuredGrid,
                  collocated_nodes: Sequence[int]) -> None:
    """
    Copies the fields from the old mesh to the new one.
    Point data will be duplicated for collocated nodes.
    :param old_mesh: The mesh before the split.
    :param new_mesh: The mesh after the split. Will receive the fields in place.
    :param collocated_nodes: New index to old index.
    :return: None
    """
    # Copying the cell data.
    # The cells are the same, just their nodes support have changed.
    input_cell_data = old_mesh.GetCellData()
    for i in range(input_cell_data.GetNumberOfArrays()):
        input_array = input_cell_data.GetArray(i)
        logging.info(f"Copying cell field \"{input_array.GetName()}\".")
        new_mesh.GetCellData().AddArray(input_array)

    # Copying field data. This data is a priori not related to geometry.
    input_field_data = old_mesh.GetFieldData()
    for i in range(input_field_data.GetNumberOfArrays()):
        input_array = input_field_data.GetArray(i)
        logging.info(f"Copying field data \"{input_array.GetName()}\".")
        new_mesh.GetFieldData().AddArray(input_array)

    # Copying the point data.
    input_point_data = old_mesh.GetPointData()
    for i in range(input_point_data.GetNumberOfArrays()):
        input_array = input_point_data.GetArray(i)
        logging.info(f"Copying point field \"{input_array.GetName()}\"")
        tmp = input_array.NewInstance()
        tmp.SetName(input_array.GetName())
        tmp.SetNumberOfComponents(input_array.GetNumberOfComponents())
        tmp.SetNumberOfTuples(new_mesh.GetNumberOfPoints())
        for p in range(tmp.GetNumberOfTuples()):
            tmp.SetTuple(p, input_array.GetTuple(collocated_nodes[p]))
        new_mesh.GetPointData().AddArray(tmp)


def __perform_split(old_mesh: vtkUnstructuredGrid,
                    cell_to_node_mapping: Mapping[int, Mapping[int, int]]) -> vtkUnstructuredGrid:
    """
    Split the main 3d mesh based on the node duplication information contained in @p cell_to_node_mapping
    :param old_mesh: The main 3d mesh.
    :param cell_to_node_mapping: For each cell, gives the nodes that must be duplicated and their new index.
    :return: The main 3d mesh split at the fracture location.
    """
    added_points: Set[int] = set()
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            if i != o:
                added_points.add(o)
    num_new_points: int = old_mesh.GetNumberOfPoints() + len(added_points)

    # Creating the new points for the new mesh.
    old_points: vtkPoints = old_mesh.GetPoints()
    new_points = vtkPoints()
    new_points.SetNumberOfPoints(num_new_points)
    collocated_nodes = numpy.ones(num_new_points, dtype=int) * -1
    # Copying old points into the new container.
    for p in range(old_points.GetNumberOfPoints()):
        new_points.SetPoint(p, old_points.GetPoint(p))
        collocated_nodes[p] = p
    # Creating the new collocated/duplicated points based on the old points positions.
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            if i != o:
                new_points.SetPoint(o, old_points.GetPoint(i))
                collocated_nodes[o] = i
    collocated_nodes.flags.writeable = False

    # We are creating a new mesh.
    # The cells will be the same, except that their nodes may be duplicated or renumbered nodes.
    # In vtk, the polyhedron and the standard cells are managed differently.
    # Also, it looks like the internal representation is being modified
    # (see https://gitlab.kitware.com/vtk/vtk/-/merge_requests/9812)
    # so we'll try nothing fancy for the moment.
    # Maybe in the future using a `DeepCopy` of the vtkCellArray can be considered?
    # The cell point ids could be modified in place then.
    new_mesh = old_mesh.NewInstance()
    new_mesh.SetPoints(new_points)
    new_mesh.Allocate(old_mesh.GetNumberOfCells())

    for c in tqdm(range(old_mesh.GetNumberOfCells()), desc="Performing the mesh split"):
        node_mapping: Mapping[int, int] = cell_to_node_mapping.get(c, {})
        cell: vtkCell = old_mesh.GetCell(c)
        cell_type: int = cell.GetCellType()
        # For polyhedron, we'll manipulate the face stream directly.
        if cell_type == VTK_POLYHEDRON:
            face_stream = vtkIdList()
            old_mesh.GetFaceStream(c, face_stream)
            new_face_nodes: List[List[int]] = []
            for face_nodes in FaceStream.build_from_vtk_id_list(face_stream).face_nodes:
                new_point_ids = []
                for current_point_id in face_nodes:
                    new_point_id: int = node_mapping.get(current_point_id, current_point_id)
                    new_point_ids.append(new_point_id)
                new_face_nodes.append(new_point_ids)
            new_mesh.InsertNextCell(cell_type, to_vtk_id_list(FaceStream(new_face_nodes).dump()))
        else:
            # For the standard cells, we extract the point ids of the cell directly.
            # Then the values will be (potentially) overwritten in place, before being sent back into the cell.
            cell_point_ids: vtkIdList = cell.GetPointIds()
            for i in range(cell_point_ids.GetNumberOfIds()):
                current_point_id: int = cell_point_ids.GetId(i)
                new_point_id: int = node_mapping.get(current_point_id, current_point_id)
                cell_point_ids.SetId(i, new_point_id)
            new_mesh.InsertNextCell(cell_type, cell_point_ids)

    __copy_fields(old_mesh, new_mesh, collocated_nodes)

    return new_mesh


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

    is_node_duplicated = numpy.zeros(mesh_points.GetNumberOfPoints(), dtype=bool)  # defaults to False
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            if not is_node_duplicated[i]:
                is_node_duplicated[i] = i != o

    # Some elements can have all their nodes not duplicated.
    # In this case, it's mandatory not get rid of this element
    # because the neighboring 3d elements won't follow.
    face_nodes: List[Collection[int]] = []
    discarded_face_nodes: Set[Iterable[int]] = set()
    for ns in fracture_info.face_nodes:
        if any(map(is_node_duplicated.__getitem__, ns)):
            face_nodes.append(ns)
        else:
            discarded_face_nodes.add(ns)

    if discarded_face_nodes:
        # tmp = []
        # for dfns in discarded_face_nodes:
        #     tmp.append(", ".join(map(str, dfns)))
        msg: str = "(" + '), ('.join(map(lambda dfns: ", ".join(map(str, dfns)), discarded_face_nodes)) + ")"
        # logging.info(f"The {len(tmp)} faces made of nodes ({'), ('.join(tmp)}) were/was discarded from the fracture mesh because none of their/its nodes were duplicated.")
        # print(f"The {len(tmp)} faces made of nodes ({'), ('.join(tmp)}) were/was discarded from the fracture mesh because none of their/its nodes were duplicated.")
        print(f"The faces made of nodes [{msg}] were/was discarded from the fracture mesh because none of their/its nodes were duplicated.")

    fracture_nodes_tmp = numpy.ones(mesh_points.GetNumberOfPoints(), dtype=int) * -1
    for ns in face_nodes:
        for n in ns:
            fracture_nodes_tmp[n] = n
    fracture_nodes: Collection[int] = tuple(filter(lambda n: n > -1, fracture_nodes_tmp))
    num_points: int = len(fracture_nodes)
    points = vtkPoints()
    points.SetNumberOfPoints(num_points)
    node_3d_to_node_2d: Dict[int, int] = {}  # Building the node mapping, from 3d mesh nodes to 2d fracture nodes.
    for i, n in enumerate(fracture_nodes):
        coords: Tuple[float, float, float] = mesh_points.GetPoint(n)
        points.SetPoint(i, coords)
        node_3d_to_node_2d[n] = i

    polygons = vtkCellArray()
    for ns in face_nodes:
        polygon = vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(len(ns))
        for i, n in enumerate(ns):
            polygon.GetPointIds().SetId(i, node_3d_to_node_2d[n])
        polygons.InsertNextCell(polygon)

    buckets: Dict[int, Set[int]] = defaultdict(set)
    for node_mapping in cell_to_node_mapping.values():
        for i, o in node_mapping.items():
            k: Optional[int] = node_3d_to_node_2d.get(min(i, o))
            if k is not None:
                buckets[k].update((i, o))

    assert set(buckets.keys()) == set(range(num_points))
    max_duplicated_nodes: int = max(map(len, buckets.values())) if buckets.values() else 0
    collocated_nodes = numpy.ones((num_points, max_duplicated_nodes), dtype=int) * -1
    for i, bucket in buckets.items():
        for j, val in enumerate(bucket):
            collocated_nodes[i, j] = val
    array = numpy_to_vtk(collocated_nodes, array_type=VTK_ID_TYPE)
    array.SetName("collocated_nodes")

    fracture_mesh = vtkUnstructuredGrid()  # We could be using vtkPolyData, but it's not supported by GEOS for now.
    fracture_mesh.SetPoints(points)
    if polygons.GetNumberOfCells() > 0:
        fracture_mesh.SetCells([VTK_POLYGON] * polygons.GetNumberOfCells(), polygons)
    fracture_mesh.GetPointData().AddArray(array)
    return fracture_mesh


def __split_mesh_on_fracture(mesh: vtkUnstructuredGrid,
                             options: Options) -> Tuple[vtkUnstructuredGrid, vtkUnstructuredGrid]:
    fracture: FractureInfo = build_fracture_info(mesh, options)
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
