from collections import defaultdict
from dataclasses import dataclass
from typing import (
    Collection,
    Dict,
    Iterable,
    List,
    Sequence,
)

from vtkmodules.vtkCommonCore import (
    vtkIdList,
)

import networkx

from .vtk_utils import (
    vtk_iter,
)


@dataclass(frozen=True)
class Options:
    dummy: float


@dataclass(frozen=True)
class Result:
    dummy: float


def parse_face_stream(ids: vtkIdList) -> Sequence[Sequence[int]]:
    """
    Parses the face stream raw information and converts it into a tuple of tuple of integers,
    each tuple of integer being the nodes of a face.
    :param ids: The raw vtk face stream.
    :return: The tuple of tuple of integers.
    """
    result = []
    it = vtk_iter(ids)
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
    assert sum(map(len, result)) + len(result) + 1 == ids.GetNumberOfIds()

    return tuple(result)


class FaceStream:
    """
    Helper class to manipulate the vtk face streams.
    """
    def __init__(self, data: Sequence[Sequence[int]]):
        # self.__data contains the list of faces nodes, like it appears in vtk face streams.
        # Except that the additional size information is removed
        # in favor of the __len__ of the containers.
        self.__data: Sequence[Sequence[int]] = data

    @staticmethod
    def build_from_vtk_id_list(ids: vtkIdList):
        """
        Builds a FaceStream from the raw vtk face stream.
        :param ids: The vtk face stream.
        :return: A new FaceStream instance.
        """
        return FaceStream(parse_face_stream(ids))

    @property
    def face_nodes(self) -> Iterable[Sequence[int]]:
        """
        Iterate on the nodes of all the faces.
        :return: An iterator.
        """
        return iter(self.__data)

    @property
    def num_faces(self) -> int:
        """
        Number of faces in the face stream
        :return: An integer
        """
        return len(self.__data)

    @property
    def support_point_ids(self) -> Collection[int]:
        """
        The list of all (unique) support points of the face stream, in no specific order.
        :return: The set of all the point ids.
        """
        tmp: List[int] = []
        for nodes in self.face_nodes:
            tmp += nodes
        return frozenset(tmp)

    @property
    def num_support_points(self) -> int:
        """
        The number of unique support nodes of the polyhedron.
        :return: An integer.
        """
        return len(self.support_point_ids)

    def __getitem__(self, face_index) -> Sequence[int]:
        """
        The support point ids for the `face_index` face.
        :param face_index: The face index (within the face stream).
        :return: A tuple containing all the point ids.
        """
        return self.__data[face_index]

    def flip_faces(self, face_indices):
        """
        Returns a new FaceStream instance with the face indices defined in face_indices flipped.,
        :param face_indices: The faces (local) indices to flip.
        :return: A newly created instance.
        """
        result = []
        for face_index, face_nodes in enumerate(self.__data):
            result.append(tuple(reversed(face_nodes)) if face_index in face_indices else face_nodes)
        return FaceStream(tuple(result))

    def dump(self) -> Sequence[int]:
        """
        Returns the face stream awaited by vtk, but in a python container.
        The content can be used, once converted to a vtkIdList, to define another polyhedron in vtk.
        :return: The face stream in a python container.
        """
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


def build_face_to_face_connectivity_through_edges(face_stream: FaceStream, add_compatibility=False) -> networkx.Graph:
    """
    Given a face stream/polyhedron, builds the connections between the faces.
    Those connections happen when two faces share an edge.
    :param face_stream: The face stream description of the polyhedron.
    :param add_compatibility: Two faces are considered compatible if their normals point in the same direction (inwards or outwards).
        If `add_compatibility=True`, we add a `compatible={"-", "+"}` flag on the edges
        to indicate that the two connected faces are compatible or not.
        If `add_compatibility=False`, non-compatible faces are simply not connected by any edge.
    :return: A graph which nodes are actually the faces of the polyhedron.
        Two nodes of the graph are connected if they share an edge.
    """
    edges_to_face_indices = defaultdict(list)
    for face_index, face_nodes in enumerate(face_stream.face_nodes):
        # Each edge is defined by two nodes. We do a small trick to loop on consecutive points.
        for face_indices in zip(face_nodes, face_nodes[1:] + (face_nodes[0], )):
            edges_to_face_indices[frozenset(face_indices)].append(face_index)
    # We are doing here some small validations w.r.t. the connections of the faces
    # which may only make sense in the context of numerical simulations.
    # As such, an error will be thrown in case the polyhedron is not closed.
    # So there may be a lack of absolute genericity, and the code may evolve if needed.
    for face_indices in edges_to_face_indices.values():
        assert len(face_indices) == 2
    # Computing the graph degree for validation
    degrees: Dict[int, int] = defaultdict(int)
    for face_indices in edges_to_face_indices.values():
        for face_index in face_indices:
            degrees[face_index] += 1
    for face_index, degree in degrees.items():
        assert len(face_stream[face_index]) == degree
    # Validation that there is one unique edge connecting two faces.
    face_indices_to_edge_index = defaultdict(list)
    for edge_index, face_indices in edges_to_face_indices.items():
        face_indices_to_edge_index[frozenset(face_indices)].append(edge_index)
    for edge_indices in face_indices_to_edge_index.values():
        assert len(edge_indices) == 1
    # Connecting the faces. Neighbor faces with consistent normals (i.e. facing both inward or outward)
    # will be connected. This will allow us to extract connected components with consistent orientations.
    # Another step will then reconcile all the components such that all the normals of the cell
    # will consistently point outward.
    graph = networkx.Graph()
    graph.add_nodes_from(range(face_stream.num_faces))
    for edge, face_indices in edges_to_face_indices.items():
        face_index_0, face_index_1 = face_indices
        face_nodes_0 = face_stream[face_index_0] + (face_stream[face_index_0][0],)
        face_nodes_1 = face_stream[face_index_1] + (face_stream[face_index_1][0],)
        node_0, node_1 = edge
        order_0 = 1 if face_nodes_0[face_nodes_0.index(node_0) + 1] == node_1 else -1
        order_1 = 1 if face_nodes_1[face_nodes_1.index(node_0) + 1] == node_1 else -1
        # Same order of nodes means that the normals of the faces are not both in the same "direction" (inward or outward).
        if order_0 * order_1 == 1:
            if add_compatibility:
                graph.add_edge(face_index_0, face_index_1, compatible="-")
        else:
            if add_compatibility:
                graph.add_edge(face_index_0, face_index_1, compatible="+")
            else:
                graph.add_edge(face_index_0, face_index_1)
    return graph
