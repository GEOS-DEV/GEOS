from collections import deque, defaultdict
from dataclasses import dataclass
import logging
from typing import List, Tuple, Iterator, Sequence, Iterable

import numpy

from vtkmodules.vtkCommonCore import (
    vtkIdList,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_POLYHEDRON,
    vtkCellArray,
    vtkUnstructuredGrid,
)

import vtk
import networkx


from . import vtk_utils


@dataclass(frozen=True)
class Options:
    dummy: float


@dataclass(frozen=True)
class Result:
    dummy: float


def _iter(id_list: vtkIdList) -> Iterator[int]:  # TODO duplicated
    """
    Utility function transforming a vtkIdList into an iterable to be used for building built-ins python containers.
    :param id_list: the vtkIdList.
    :return: The iterator.
    """
    for i in range(id_list.GetNumberOfIds()):
        yield id_list.GetId(i)


def to_vtk_id_list(data):  # TODO move to utility
    result = vtkIdList()
    result.Allocate(len(data))
    for d in data:
        result.InsertNextId(d)
    return result


def parse_face_stream(ids: vtkIdList):  # TODO move to FaceStream.build_from_vtk_id_list
    """
    Parses the face stream raw information and converts it into a tuple of tuple of integers,
    each tuple of integer being the nodes of a face.
    :param ids:
    :return:
    """
    result = []
    it = _iter(ids)
    num_faces = next(it)
    try:
        while True:
            num_nodes = next(it)
            tmp = []
            for i in range(num_nodes):
                tmp.append(next(it))
            result.append(tuple(tmp))
    except StopIteration:
        pass
    assert len(result) == num_faces
    return tuple(result)


class FaceStream:
    """
    Helper class to manipulate the vtk face streams.
    """
    def __init__(self, data):
        self.__data = data

    @staticmethod
    def build_from_vtk_id_list(ids: vtkIdList):
        return FaceStream(parse_face_stream(ids))

    @property
    def face_nodes(self) -> Iterable[Sequence[int]]:
        return iter(self.__data)

    def __getitem__(self, face_index):
        return self.__data[face_index]

    def flip_faces(self, face_indices):
        result = []
        for face_index, face_nodes in enumerate(self.__data):
            result.append(tuple(reversed(face_nodes)) if face_index in face_indices else face_nodes)
        return FaceStream(tuple(result))

    def dump(self) -> Sequence[int]:
        result = [len(self.__data)]
        for face_nodes in self.__data:
            result.append(len(face_nodes))
            result += face_nodes
        return tuple(result)

    def __repr__(self):
        result = [str(len(self.__data))]
        for face_nodes in self.__data:
            result.append(str(len(face_nodes)))
            result.append(", ".join(map(str, face_nodes)))
        return ",\n".join(result)


def __analyze_face_stream(face_stream: FaceStream):
    edges_to_face_indices = defaultdict(list)
    for face_index, face_nodes in enumerate(face_stream.face_nodes):
        shifted_face_nodes = deque(face_nodes)
        shifted_face_nodes.rotate(-1)
        for face_indices in zip(face_nodes, shifted_face_nodes):
            edges_to_face_indices[frozenset(face_indices)].append(face_index)
    for face_indices in edges_to_face_indices.values():
        assert len(face_indices) == 2
    # TODO check that there is only one link (edge) between two faces.
    # Computing the graph degree for validation
    degrees = defaultdict(int)
    for face_indices in edges_to_face_indices.values():
        for face_index in face_indices:
            degrees[face_index] += 1
    for face_index, degree in degrees.items():
        assert len(face_stream[face_index]) == degree
    # Connecting the faces. Neighbor faces with consistent normals (i.e. facing both inward or outward)
    # will be connected. This will allow us to extract connected components with consistent orientations.
    # Another step will then reconcile all the components such that all the normals of the cell
    # will consistently point outward.
    graph = networkx.Graph()
    for edge, face_indices in edges_to_face_indices.items():
        face_index_0, face_index_1 = face_indices
        face_nodes_0 = face_stream[face_index_0] + (face_stream[face_index_0][0],)
        face_nodes_1 = face_stream[face_index_1] + (face_stream[face_index_1][0],)
        node_0, node_1 = edge
        order_0 = 1 if face_nodes_0[face_nodes_0.index(node_0) + 1] == node_1 else -1
        order_1 = 1 if face_nodes_1[face_nodes_1.index(node_0) + 1] == node_1 else -1
        # Same order of nodes means that the normals of the faces or not both (in|out)ward.
        if order_0 * order_1 == 1:
            logging.warning(f"Inconsistent order between faces {face_index_0} and {face_index_1}.")
        else:
            logging.warning(f"Consistent order between faces {face_index_0} and {face_index_1}.")
            graph.add_edge(face_index_0, face_index_1)
    # Consider adding the wrong and bad connections between the faces in thegraph,
    # then return the graph and let the algo continue.
    return tuple(networkx.connected_components(graph))
    # TODO I should be able to colour the graphs of the connected components with two colors.
    #      Then I have two choices.


def __check(mesh, options: Options) -> Result:
    points = mesh.GetPoints()
    output_mesh = vtkUnstructuredGrid()
    output_mesh.SetPoints(points)

    for ic in range(mesh.GetNumberOfCells()):
        c = mesh.GetCell(ic)
        if not c.GetCellType() == VTK_POLYHEDRON:
            continue
        pt_ids = vtk.vtkIdList()
        mesh.GetFaceStream(ic, pt_ids)
        face_stream = FaceStream.build_from_vtk_id_list(pt_ids)
        connected_components = __analyze_face_stream(face_stream)
        assert len(connected_components) in (1, 2)  # TODO Some tricky cells may have more than 2?
        print(connected_components)
        # TODO flip the faces, check if faces are now directed outwards.
        # TODO One trick, for prisms, could be to have the "big faces" not facing each other.
        cc = connected_components  # alias
        if len(cc) > 1:  # TODO this is totally random (and wrong)...
            fs = face_stream.flip_faces(cc[0] if len(cc[0]) < len(cc[1]) else cc[1])
        else:
            fs = face_stream
        new_face_ids = to_vtk_id_list(fs.dump())
        output_mesh.InsertNextCell(VTK_POLYHEDRON, new_face_ids)
    return output_mesh


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
