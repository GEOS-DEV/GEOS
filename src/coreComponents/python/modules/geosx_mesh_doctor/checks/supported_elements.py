from dataclasses import dataclass
import logging
import multiprocessing
from typing import Dict, Sequence, FrozenSet, Set

from tqdm import tqdm

import networkx
import numpy

# from vtkmodules.vtkFiltersGeneral import (
#     vtkCellValidator
# )
# from vtkmodules.vtkCommonCore import (
#     vtkOutputWindow,
#     vtkFileOutputWindow
# )
# from vtk.util.numpy_support import (
#     vtk_to_numpy,
# )
from vtkmodules.vtkCommonCore import (
    vtkIdList,
)
from vtkmodules.vtkCommonDataModel import (
    vtkCellTypes,
    VTK_HEXAGONAL_PRISM,
    VTK_HEXAHEDRON,
    VTK_PENTAGONAL_PRISM,
    VTK_POLYHEDRON,
    VTK_PYRAMID,
    VTK_TETRA,
    VTK_VOXEL,
    VTK_WEDGE,
)

from . import vtk_utils
from .vtk_polyhedron import build_cell_graph, FaceStream

@dataclass(frozen=True)
class Options:
    num_proc: int
    chunk_size: int


@dataclass(frozen=True)
class Result:
    unsupported_std_elements_types: Set[int]  # list of unsupported types
    unsupported_polyhedron_elements: Sequence[int]  # list of polyhedron elements that could not be converted to supported std elements


# def __check(mesh) -> Result:
#     # Dealing with std elements (i.e. non polyhedron)
#     cell_types = vtkCellTypes()
#     mesh.GetCellTypes(cell_types)
#     cell_types = _convert(cell_types)
#     supported_cell_types = {
#         VTK_HEXAGONAL_PRISM,
#         VTK_HEXAHEDRON,
#         VTK_PENTAGONAL_PRISM,
#         VTK_POLYHEDRON,
#         VTK_PYRAMID,
#         VTK_TETRA,
#         VTK_VOXEL,
#         VTK_WEDGE
#     }
#     unsupported_std_elements_types = cell_types - supported_cell_types
#     # Dealing with polyhedron elements.
#     unsupported_polyhedron_elements = []
#     is_polyhedron_supported = IsPolyhedronConvertible()
#     for ic in tqdm(range(mesh.GetNumberOfCells())):
#         if mesh.GetCellType(ic) == VTK_POLYHEDRON:
#             pt_ids = vtkIdList()
#             mesh.GetFaceStream(ic, pt_ids)
#             face_stream = FaceStream.build_from_vtk_id_list(pt_ids)
#             converted_type_name = is_polyhedron_supported(mesh, face_stream)
#             if converted_type_name:
#                 # logging.debug(f"Polyhedron cell {ic} can be converted into \"{converted_type_name}\"")
#                 pass
#             else:
#                 unsupported_polyhedron_elements.append(ic)
#     return Result(unsupported_std_elements_types=unsupported_std_elements_types,
#                   unsupported_polyhedron_elements=unsupported_polyhedron_elements)



# def _convert(cell_types) -> FrozenSet[int]:  # TODO move to utilities
#     tmp = []
#     for i in range(cell_types.GetNumberOfValues()):
#         tmp.append(cell_types.GetValue(i))
#     return frozenset(tmp)
def _convert(cell_types) -> FrozenSet[int]:  # TODO move to utilities
    tmp = []
    for i in range(cell_types.GetNumberOfTypes()):
        tmp.append(cell_types.GetCellType(i))
    return frozenset(tmp)


# class IsPolyhedronConvertible:
#     def __init__(self):
#         # self.__equivalent_graphs = dict()
#         # # Tetrahedron
#         # self.__equivalent_graphs["Tetrahedron"] = networkx.complete_graph(4)
#         def build_prism_graph(n, name):
#             tmp = networkx.cycle_graph(n)
#             for node in range(n):
#                 tmp.add_edge(node, n)
#                 tmp.add_edge(node, n + 1)
#             tmp.name = name
#             return tmp
#         # # Pyramid
#         # self.__equivalent_graphs["Pyramid"] = build_prism_graph(4)
#         # self.__equivalent_graphs["Pyramid"].remove_node(5)  # also removes the associated edges.
#         # # Prisms
#         # self.__equivalent_graphs["Wedge"] = build_prism_graph(3)
#         # self.__equivalent_graphs["Hexahedron"] = build_prism_graph(4)
#         # self.__equivalent_graphs["Prism5"] = build_prism_graph(5)
#         # self.__equivalent_graphs["Prism6"] = build_prism_graph(6)
#         # self.__equivalent_graphs["Prism7"] = build_prism_graph(7)
#         # self.__equivalent_graphs["Prism8"] = build_prism_graph(8)
#         # self.__equivalent_graphs["Prism9"] = build_prism_graph(9)
#         # self.__equivalent_graphs["Prism10"] = build_prism_graph(10)
#         # self.__equivalent_graphs["Prism11"] = build_prism_graph(11)
#
#         tet_graph = networkx.complete_graph(4)
#         tet_graph.name = "Tetrahedron"
#         pyr_graph = build_prism_graph(4, "Pyramid")
#         pyr_graph.remove_node(5)  # also removes the associated edges.
#         self.__test = {
#             4: (tet_graph,),
#             5: (pyr_graph, build_prism_graph(3, "Wedge")),
#             6: (build_prism_graph(4, "Hexahedron"),),
#             7: (build_prism_graph(5, "Prism5"),),
#             8: (build_prism_graph(6, "Prism6"),),
#             9: (build_prism_graph(7, "Prism7"),),
#             10: (build_prism_graph(8, "Prism8"),),
#             11: (build_prism_graph(9, "Prism9"),),
#             12: (build_prism_graph(10, "Prism10"),),
#             13: (build_prism_graph(11, "Prism11"),),
#         }
#         # In order not to test all the combinations all the time.
#         # self.__face_stream_sig = {
#         #     frozenset({3: 4}): VTK_TETRA,
#         #     frozenset({3: 2, 4: 3}): VTK_WEDGE,
#         #     frozenset({4: 6}): VTK_HEXAHEDRON,
#         #     frozenset({4: 5, 5: 2}): VTK_PENTAGONAL_PRISM,
#         #     frozenset({4: 6, 6: 2}): VTK_PENTAGONAL_PRISM,
#         # } # TODO Prims 7+!
#
#     def __call__(self, mesh, face_stream) -> str:
#         cell_graph = build_cell_graph(face_stream, add_compatibility=True)
#         for reference_graph in self.__test[cell_graph.order()]:
#             if networkx.is_isomorphic(reference_graph, cell_graph):
#                 return str(reference_graph.name)
#         return ""
#     # def __call__(self, mesh, face_stream) -> str:
#     #     g = build_cell_graph(face_stream, add_compatibility=True)
#     #     for name, reference_graph in self.__equivalent_graphs.items():  # TODO improve if too slow
#     #         if networkx.is_isomorphic(reference_graph, g):
#     #             return name
#     #     return ""


MESH = None  # for multiprocessing, vtkUnstructuredGrid cannot be pickled. Let's use a global variable instead.

class IsPolyhedronConvertible:
    def __init__(self):
        def build_prism_graph(n, name):
            tmp = networkx.cycle_graph(n)
            for node in range(n):
                tmp.add_edge(node, n)
                tmp.add_edge(node, n + 1)
            tmp.name = name
            return tmp
        # building the reference graphs
        tet_graph = networkx.complete_graph(4)
        tet_graph.name = "Tetrahedron"
        pyr_graph = build_prism_graph(4, "Pyramid")
        pyr_graph.remove_node(5)  # also removes the associated edges.
        self.__reference_graphs = {
            4: (tet_graph,),
            5: (pyr_graph, build_prism_graph(3, "Wedge")),
            6: (build_prism_graph(4, "Hexahedron"),),
            7: (build_prism_graph(5, "Prism5"),),
            8: (build_prism_graph(6, "Prism6"),),
            9: (build_prism_graph(7, "Prism7"),),
            10: (build_prism_graph(8, "Prism8"),),
            11: (build_prism_graph(9, "Prism9"),),
            12: (build_prism_graph(10, "Prism10"),),
            13: (build_prism_graph(11, "Prism11"),),
        }

    def __is_polyhedron_supported(self, face_stream) -> str:
        cell_graph = build_cell_graph(face_stream, add_compatibility=True)
        for reference_graph in self.__reference_graphs[cell_graph.order()]:
            if networkx.is_isomorphic(reference_graph, cell_graph):
                return str(reference_graph.name)
        return ""

    def __call__(self, ic) -> int:
        """
        Checks if a vtk polyhedron cell can be converted into a supported GEOSX element.
        :param ic: The index of the vtk element.
        :return: -1 if the polyhedron vtk element can be converted into a supported element type. -1 otherwise.
        """
        if MESH.GetCellType(ic) != VTK_POLYHEDRON:
            return -1
        pt_ids = vtkIdList()
        MESH.GetFaceStream(ic, pt_ids)
        face_stream = FaceStream.build_from_vtk_id_list(pt_ids)
        converted_type_name = self.__is_polyhedron_supported(face_stream)
        if converted_type_name:
            logging.warning(f"Polyhedron cell {ic} can be converted into \"{converted_type_name}\"")
            return -1
        else:
            logging.warning(f"Polyhedron cell {ic} cannot be converted into any supported element.")
            return ic


def __check(mesh, options: Options) -> Result:
    # Dealing with std elements (i.e. non polyhedron)
    cell_types = vtkCellTypes()
    mesh.GetCellTypes(cell_types)
    cell_types = _convert(cell_types)
    # TODO use cell_types = mesh.GetDistinctCellTypesArray()
    # cell_types = _convert(mesh.GetDistinctCellTypesArray())
    supported_cell_types = {
        VTK_HEXAGONAL_PRISM,
        VTK_HEXAHEDRON,
        VTK_PENTAGONAL_PRISM,
        VTK_POLYHEDRON,
        VTK_PYRAMID,
        VTK_TETRA,
        VTK_VOXEL,
        VTK_WEDGE
    }
    unsupported_std_elements_types = cell_types - supported_cell_types

    # Dealing with polyhedron elements.
    global MESH  # for multiprocessing, vtkUnstructuredGrid cannot be pickled. Let's use a global variable instead.
    MESH = mesh
    num_cells = mesh.GetNumberOfCells()
    result = numpy.ones(num_cells, dtype=int) * -1
    f = IsPolyhedronConvertible()
    with multiprocessing.Pool(processes=options.num_proc) as pool:
        generator = pool.imap_unordered(f, range(num_cells), chunksize=options.chunk_size)
        for i, val in enumerate(tqdm(generator, total=num_cells)):
            result[i] = val
    unsupported_polyhedron_elements = [i for i in result if i > -1]
    return Result(unsupported_std_elements_types=unsupported_std_elements_types,
                  unsupported_polyhedron_elements=unsupported_polyhedron_elements)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
