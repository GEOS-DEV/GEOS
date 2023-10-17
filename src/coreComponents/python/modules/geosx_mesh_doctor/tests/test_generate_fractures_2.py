from dataclasses import dataclass

from typing import (
    Tuple,
    Iterable,
    Iterator,
    Dict,
    Mapping,
    FrozenSet,
    List,
    Set,
    Sequence,
    Collection,
)

import numpy

import pytest

from vtkmodules.vtkCommonDataModel import (
    vtkUnstructuredGrid,
)


from vtkmodules.util.numpy_support import (
    numpy_to_vtk,
    vtk_to_numpy,
)
import sys
sys.path.append("..")

# import sys
# sys.path.insert(0, "/Users/j0436735/CLionProjects/GEOS/src/coreComponents/python/modules/geosx_mesh_doctor")
# sys.path.insert(0, "/Users/j0436735/CLionProjects/GEOS/src/coreComponents/python/modules/geosx_mesh_doctor/tests")

from checks.generate_cube import build_rectilinear_blocks_mesh, XYZ
from checks.generate_fractures_2 import __split_mesh_on_fracture, Options
# from ..checks.generate_cube import build_rectilinear_blocks_mesh, XYZ
# from ..checks.generate_fractures_2 import __split_mesh_on_fracture, Options


@dataclass(frozen=True)
class TestCase:
    input_mesh: vtkUnstructuredGrid
    options: Options
    collocated_nodes: Tuple[Tuple[int, ...], ...]
    result: Tuple[int, ...]


def __build_test_case(xs: Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray],
                      attribute: Iterable[int]):
    xyz = XYZ(*xs)

    mesh: vtkUnstructuredGrid = build_rectilinear_blocks_mesh((xyz, ))

    ref = numpy.array(attribute, dtype=int)
    assert len(ref) == mesh.GetNumberOfCells()
    attr = numpy_to_vtk(ref)
    attr.SetName("attribute")
    mesh.GetCellData().AddArray(attr)

    options = Options(policy="field",
                      field="attribute",
                      field_type="cells",
                      field_values=frozenset(attribute),
                      vtk_output=None,
                      vtk_fracture_output=None,
                      split_on_domain_boundary=True)
    return mesh, options


class Incrementor:
    def __init__(self, start):
        self.__val = start

    def next(self, num: int) -> Iterable[int]:
        self.__val += num
        return range(self.__val - num, self.__val)


def __generate_test_data() -> Iterator[Tuple[vtkUnstructuredGrid, Options]]:
    tmp0 = numpy.arange(3, dtype=float)

    # Split in 2
    mesh, options = __build_test_case((tmp0, tmp0, tmp0), (0, 1, 0, 1, 0, 1, 0, 1))
    yield TestCase(input_mesh=mesh, options=options,
                   collocated_nodes=tuple(map(lambda i: (1 + 3 * i, 27 + i), range(9))),
                   result=(9 * 4, 8, 9, 4))

    # Split in 3
    inc = Incrementor(27)
    collocated_nodes = (
        (1, *inc.next(1)),
        (3, *inc.next(1)),
        (4, *inc.next(2)),
        (7, *inc.next(1)),
        (1 + 9, *inc.next(1)),
        (3 + 9, *inc.next(1)),
        (4 + 9, *inc.next(2)),
        (7 + 9, *inc.next(1)),
        (1 + 18, *inc.next(1)),
        (3 + 18, *inc.next(1)),
        (4 + 18, *inc.next(2)),
        (7 + 18, *inc.next(1)),
    )
    mesh, options = __build_test_case((tmp0, tmp0, tmp0), (0, 1, 2, 1, 0, 1, 2, 1))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(9 * 4 + 6, 8, 12, 6))

    # Split in 8
    inc = Incrementor(27)
    collocated_nodes = (
        (1, *inc.next(1)),
        (3, *inc.next(1)),
        (4, *inc.next(3)),
        (5, *inc.next(1)),
        (7, *inc.next(1)),
        (0 + 9, *inc.next(1)),
        (1 + 9, *inc.next(3)),
        (2 + 9, *inc.next(1)),
        (3 + 9, *inc.next(3)),
        (4 + 9, *inc.next(7)),
        (5 + 9, *inc.next(3)),
        (6 + 9, *inc.next(1)),
        (7 + 9, *inc.next(3)),
        (8 + 9, *inc.next(1)),
        (1 + 18, *inc.next(1)),
        (3 + 18, *inc.next(1)),
        (4 + 18, *inc.next(3)),
        (5 + 18, *inc.next(1)),
        (7 + 18, *inc.next(1)),
    )
    mesh, options = __build_test_case((tmp0, tmp0, tmp0), range(8))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(8 * 8, 8, 3 * 3 * 3 - 8, 12))


@pytest.mark.parametrize("expected", __generate_test_data())
def test_generate_fracture(expected: TestCase):
    main_mesh, fracture_mesh = __split_mesh_on_fracture(expected.input_mesh, expected.options)
    assert main_mesh.GetNumberOfPoints() == expected.result[0]
    assert main_mesh.GetNumberOfCells() == expected.result[1]
    assert fracture_mesh.GetNumberOfPoints() == expected.result[2]
    assert fracture_mesh.GetNumberOfCells() == expected.result[3]

    collocated_nodes: numpy.ndarray = vtk_to_numpy(fracture_mesh.GetPointData().GetArray("collocated_nodes"))
    res = tuple(sorted(map(lambda bucket: tuple(sorted(filter(lambda i: i != -1, bucket))), collocated_nodes)))
    assert res == expected.collocated_nodes
    assert len(collocated_nodes) == expected.result[2]
