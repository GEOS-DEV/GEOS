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

    hex = vtkHexahedron()
    hex.GetPointIds().SetId(0, 0)
    hex.GetPointIds().SetId(1, 1)
    hex.GetPointIds().SetId(2, 3)  # Intentionally wrong
    hex.GetPointIds().SetId(3, 2)  # Intentionally wrong
    hex.GetPointIds().SetId(4, 4)
    hex.GetPointIds().SetId(5, 5)
    hex.GetPointIds().SetId(6, 6)
    hex.GetPointIds().SetId(7, 7)
    cells.InsertNextCell(hex)

    mesh = vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells(cell_types, cells)

    result = __check(mesh, Options(tolerance=0.))

    assert len(result.intersecting_faces_elements) == 1
    assert result.intersecting_faces_elements[0] == 0
