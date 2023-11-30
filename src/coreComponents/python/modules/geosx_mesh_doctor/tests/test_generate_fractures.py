from dataclasses import dataclass

from typing import (
    Tuple,
    Iterable,
    Iterator,
)

import numpy

import pytest

from vtkmodules.vtkCommonDataModel import (
    vtkUnstructuredGrid,
    VTK_POLYHEDRON,
    VTK_HEXAHEDRON,
)
from vtkmodules.util.numpy_support import (
    numpy_to_vtk,
)

from checks.vtk_utils import (
    to_vtk_id_list,
)

from checks.check_fractures import format_collocated_nodes
from checks.generate_cube import build_rectilinear_blocks_mesh, XYZ
from checks.generate_fractures import __split_mesh_on_fracture, Options


@dataclass(frozen=True)
class TestCase:
    __test__ = False
    input_mesh: vtkUnstructuredGrid
    options: Options
    collocated_nodes: Tuple[Tuple[int, ...], ...]
    result: Tuple[int, ...]


def __build_test_case(xs: Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray],
                      attribute: Iterable[int],
                      field_values: Iterable[int] = None):
    xyz = XYZ(*xs)

    mesh: vtkUnstructuredGrid = build_rectilinear_blocks_mesh((xyz, ))

    ref = numpy.array(attribute, dtype=int)
    assert len(ref) == mesh.GetNumberOfCells()
    attr = numpy_to_vtk(ref)
    attr.SetName("attribute")
    mesh.GetCellData().AddArray(attr)

    if field_values is None:
        fv = frozenset(attribute)
    else:
        fv = field_values

    options = Options(policy="field",
                      field="attribute",
                      field_type="cells",
                      field_values=fv,
                      vtk_output=None,
                      vtk_fracture_output=None,
                      split_on_domain_boundary=True)
    return mesh, options


# Utility class to generate the new indices of the newly created collocated nodes.
class Incrementor:
    def __init__(self, start):
        self.__val = start

    def next(self, num: int) -> Iterable[int]:
        self.__val += num
        return range(self.__val - num, self.__val)


def __generate_test_data() -> Iterator[Tuple[vtkUnstructuredGrid, Options]]:
    two_nodes = numpy.arange(2, dtype=float)
    three_nodes = numpy.arange(3, dtype=float)
    four_nodes = numpy.arange(4, dtype=float)

    # Split in 2
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), (0, 1, 0, 1, 0, 1, 0, 1))
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
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), (0, 1, 2, 1, 0, 1, 2, 1))
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
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), range(8))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(8 * 8, 8, 3 * 3 * 3 - 8, 12))

    # Straight notch
    inc = Incrementor(27)
    collocated_nodes = (
        (1, *inc.next(1)),
        (4,),
        (1 + 9, *inc.next(1)),
        (4 + 9,),
        (1 + 18, *inc.next(1)),
        (4 + 18,),
    )
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), (0, 1, 2, 2, 0, 1, 2, 2), field_values=(0, 1))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(3 * 3 * 3 + 3, 8, 6, 2))

    # L-shaped notch
    inc = Incrementor(27)
    collocated_nodes = (
        (1, *inc.next(1)),
        (4, *inc.next(1)),
        (7, *inc.next(1)),
        (1 + 9, *inc.next(1)),
        (4 + 9,),
        (7 + 9,),
        (1 + 18, *inc.next(1)),
        (4 + 18,),
    )
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), (0, 1, 0, 1, 0, 1, 2, 2), field_values=(0, 1))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(3 * 3 * 3 + 5, 8, 8, 3))

    # 3x1x1 split
    inc = Incrementor(2 * 2 * 4)
    collocated_nodes = (
        (1, *inc.next(1)),
        (2, *inc.next(1)),
        (5, *inc.next(1)),
        (6, *inc.next(1)),
        (1 + 8, *inc.next(1)),
        (2 + 8, *inc.next(1)),
        (5 + 8, *inc.next(1)),
        (6 + 8, *inc.next(1)),
    )
    mesh, options = __build_test_case((four_nodes, two_nodes, two_nodes), (0, 1, 2))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(6 * 4, 3, 2 * 4, 2))

    # Discarded fracture element if no node duplication.
    collocated_nodes = ()
    mesh, options = __build_test_case((three_nodes, four_nodes, four_nodes), [0, ] * 8 + [1, 2] + [0, ] * 8, field_values=(1, 2))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(3 * 4 * 4, 2 * 3 * 3, 0, 0))

    # Fracture on a corner
    inc = Incrementor(3 * 4 * 4)
    collocated_nodes = (
        (1 + 12,),
        (4 + 12,),
        (7 + 12,),
        (1 + 12 * 2, *inc.next(1)),
        (4 + 12 * 2, *inc.next(1)),
        (7 + 12 * 2,),
        (1 + 12 * 3, *inc.next(1)),
        (4 + 12 * 3, *inc.next(1)),
        (7 + 12 * 3,),
    )
    mesh, options = __build_test_case((three_nodes, four_nodes, four_nodes), [0, ] * 6 + [1, 2, 1, 2, 0, 0, 1, 2, 1, 2, 0, 0], field_values=(1, 2))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(3 * 4 * 4 + 4, 2 * 3 * 3, 9, 4))

    # Generate mesh with 2 hexs, one being a standard hex, the other a 42 hex.
    inc = Incrementor(3 * 2 * 2)
    collocated_nodes = (
        (1, *inc.next(1)),
        (1 + 3, *inc.next(1)),
        (1 + 6, *inc.next(1)),
        (1 + 9, *inc.next(1)),
    )
    mesh, options = __build_test_case((three_nodes, two_nodes, two_nodes), (0, 1))
    polyhedron_mesh = vtkUnstructuredGrid()
    polyhedron_mesh.SetPoints(mesh.GetPoints())
    polyhedron_mesh.Allocate(2)
    polyhedron_mesh.InsertNextCell(VTK_HEXAHEDRON, to_vtk_id_list((1, 2, 5, 4, 7, 8, 10, 11)))
    poly = to_vtk_id_list([6] + [4, 0, 1, 7, 6] + [4, 1, 4, 10, 7] + [4, 4, 3, 9, 10] + [4, 3, 0, 6, 9] + [4, 6, 7, 10, 9] + [4, 1, 0, 3, 4])
    polyhedron_mesh.InsertNextCell(VTK_POLYHEDRON, poly)
    polyhedron_mesh.GetCellData().AddArray(mesh.GetCellData().GetArray("attribute"))

    yield TestCase(input_mesh=polyhedron_mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(4 * 4, 2, 4, 1))


@pytest.mark.parametrize("expected", __generate_test_data())
def test_generate_fracture(expected: TestCase):
    main_mesh, fracture_mesh = __split_mesh_on_fracture(expected.input_mesh, expected.options)
    assert main_mesh.GetNumberOfPoints() == expected.result[0]
    assert main_mesh.GetNumberOfCells() == expected.result[1]
    assert fracture_mesh.GetNumberOfPoints() == expected.result[2]
    assert fracture_mesh.GetNumberOfCells() == expected.result[3]

    res = format_collocated_nodes(fracture_mesh)
    assert res == expected.collocated_nodes
    assert len(res) == expected.result[2]
