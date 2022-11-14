from vtkmodules.vtkCommonCore import (
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_HEXAHEDRON,
    vtkCellArray,
    vtkHexahedron,
    vtkUnstructuredGrid,
)


from checks.self_intersecting_elements import Options, __check


def test_jumbled_hex():
    # creating a simple hexahedron
    points = vtkPoints()
    points.SetNumberOfPoints(8)
    points.SetPoint(0, (0, 0, 0))
    points.SetPoint(1, (1, 0, 0))
    points.SetPoint(2, (1, 1, 0))
    points.SetPoint(3, (0, 1, 0))
    points.SetPoint(4, (0, 0, 1))
    points.SetPoint(5, (1, 0, 1))
    points.SetPoint(6, (1, 1, 1))
    points.SetPoint(7, (0, 1, 1))

    cell_types = [VTK_HEXAHEDRON]
    cells = vtkCellArray()
    cells.AllocateExact(1, 8)

    tet = vtkHexahedron()
    tet.GetPointIds().SetId(0, 0)
    tet.GetPointIds().SetId(1, 1)
    tet.GetPointIds().SetId(2, 3)  # Intentionally wrong
    tet.GetPointIds().SetId(3, 2)  # Intentionally wrong
    tet.GetPointIds().SetId(4, 4)
    tet.GetPointIds().SetId(5, 5)
    tet.GetPointIds().SetId(6, 6)
    tet.GetPointIds().SetId(7, 7)
    cells.InsertNextCell(tet)

    mesh = vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells(cell_types, cells)

    result = __check(mesh, Options(tolerance=0.))

    assert len(result.jumbled_elements) == 1
    assert result.jumbled_elements[0] == 0
