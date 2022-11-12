import numpy

from vtkmodules.vtkCommonCore import (
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_TETRA,
    vtkCellArray,
    vtkTetra,
    vtkUnstructuredGrid,
)


from checks.elements_volumes import Options, __check


def test_simple_collocated_points():
    # creating a simple tetrahedron
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
    tet.GetPointIds().SetId(3, 3)
    cells.InsertNextCell(tet)

    mesh = vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells(cell_types, cells)

    result = __check(mesh, Options(min_volume=1.))

    assert len(result.element_volumes) == 1
    assert result.element_volumes[0][0] == 0
    assert abs(result.element_volumes[0][1] - 1./6.) < 10 * numpy.finfo(float).eps

    result = __check(mesh, Options(min_volume=0.))

    assert len(result.element_volumes) == 0
