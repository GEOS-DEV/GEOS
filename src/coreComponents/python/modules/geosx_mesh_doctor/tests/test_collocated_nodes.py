from typing import Iterator, Tuple

import pytest

from vtkmodules.vtkCommonCore import (
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_TETRA,
    vtkCellArray,
    vtkTetra,
    vtkUnstructuredGrid,
)

from checks.collocated_nodes import Options, __check


def get_points() -> Iterator[Tuple[vtkPoints, int]]:
    """
    Generates the data for the cases.
    One case has two nodes at the exact same position.
    The other has two differente nodes
    :return: Generator to (vtk points, number of expected duplicated locations)
    """
    for p0, p1 in ((0, 0, 0), (1, 1, 1)), ((0, 0, 0), (0, 0, 0)):
        points = vtkPoints()
        points.SetNumberOfPoints(2)
        points.SetPoint(0, p0)
        points.SetPoint(1, p1)
        num_nodes_bucket = 1 if p0 == p1 else 0
        yield points, num_nodes_bucket


@pytest.mark.parametrize("data", get_points())
def test_simple_collocated_points(data: Tuple[vtkPoints, int]):
    points, num_nodes_bucket = data

    mesh = vtkUnstructuredGrid()
    mesh.SetPoints(points)

    result = __check(mesh, Options(tolerance=1.e-12))

    assert len(result.wrong_support_elements) == 0
    assert len(result.nodes_buckets) == num_nodes_bucket
    if num_nodes_bucket == 1:
        assert len(result.nodes_buckets[0]) == points.GetNumberOfPoints()


def test_wrong_support_elements():
    points = vtkPoints()
    points.SetNumberOfPoints(4)
    points.SetPoint(0, (0, 0, 0))
    points.SetPoint(1, (1, 0, 0))
    points.SetPoint(2, (0, 1, 0))
    points.SetPoint(3, (0, 0, 1))

    cell_types = [VTK_TETRA]
    cells = vtkCellArray()
    cells.AllocateExact(1, 4)

    tet = vtkTetra()
    tet.GetPointIds().SetId(0, 0)
    tet.GetPointIds().SetId(1, 1)
    tet.GetPointIds().SetId(2, 2)
    tet.GetPointIds().SetId(3, 0)  # Intentionally wrong
    cells.InsertNextCell(tet)

    mesh = vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells(cell_types, cells)

    result = __check(mesh, Options(tolerance=1.e-12))

    assert len(result.nodes_buckets) == 0
    assert len(result.wrong_support_elements) == 1
    assert result.wrong_support_elements[0] == 0
