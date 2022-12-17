from collections import deque, defaultdict
from dataclasses import dataclass
import logging
from typing import List, Tuple, Iterator, Sequence, Iterable

import numpy

from vtkmodules.vtkCommonCore import (
    vtkIdList,
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


def parse_face_stream(ids: vtkIdList):
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
    # Computing the graph degree
    degrees = defaultdict(int)
    for face_indices in edges_to_face_indices.values():
        for face_index in face_indices:
            degrees[face_index] += 1
    for face_index, degree in degrees.items():
        assert len(face_stream[face_index]) == degree
    # for pair in edges_to_face_indices.values():
    #     assert len(pair) == 2
    #     graph.add_edge(*pair)
    # for face_index, face_nodes in enumerate(face_stream):  # networkx is overkill for this
    #     assert graph.degree(face_index) == len(face_nodes)
    graph = networkx.Graph()
    for edge, face_indices in edges_to_face_indices.items():
        face_index_0, face_index_1 = face_indices
        face_nodes_0 = face_stream[face_index_0] + (face_stream[face_index_0][0],)
        face_nodes_1 = face_stream[face_index_1] + (face_stream[face_index_1][0],)
        node_0, node_1 = edge
        order_0 = 1 if face_nodes_0[face_nodes_0.index(node_0) + 1] == node_1 else -1
        order_1 = 1 if face_nodes_1[face_nodes_1.index(node_0) + 1] == node_1 else -1
        if order_0 * order_1 == 1:
            logging.warning(f"Inconsistent order between faces {face_index_0} and {face_index_1}.")
            # graph.add_edge(face_index_0, face_index_1, connection=-1)  # TODO do not add this connection if we want it removed!
        else:
            logging.warning(f"Consistent order between faces {face_index_0} and {face_index_1}.")
            # graph.add_edge(face_index_0, face_index_1, connection=1)
            graph.add_edge(face_index_0, face_index_1)

    # for e in graph.edges:
    #     if graph.edges[e]["connection"] == -1:  # TODO Use dict and constant for connection
    #         graph.remove_edge(*e)
    return tuple(networkx.connected_components(graph))
    # return graph


def __check(mesh, options: Options) -> Result:
    num_cells = mesh.GetNumberOfCells()
    for ic in range(num_cells):
        c = mesh.GetCell(ic)
        if not c.GetCellType() == vtk.VTK_POLYHEDRON:
            continue
        pt_ids = vtk.vtkIdList()
        mesh.GetFaceStream(ic, pt_ids)
        face_stream = FaceStream.build_from_vtk_id_list(pt_ids)
        connected_components = __analyze_face_stream(face_stream)
        assert len(connected_components) in (1, 2)
        # TODO flip the faces, check if faces are now directed outwards


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
