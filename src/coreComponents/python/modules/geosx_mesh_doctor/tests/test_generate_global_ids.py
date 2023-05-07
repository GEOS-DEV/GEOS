from vtkmodules.vtkCommonCore import (
    vtkPoints, )
from vtkmodules.vtkCommonDataModel import (
    VTK_VERTEX,
    vtkCellArray,
    vtkUnstructuredGrid,
    vtkVertex,
)

from checks.generate_global_ids import __build_global_ids


def test_generate_global_ids():
    points = vtkPoints()
    points.InsertNextPoint(0, 0, 0)

    vertex = vtkVertex()
    vertex.GetPointIds().SetId(0, 0)

    vertices = vtkCellArray()
    vertices.InsertNextCell(vertex)

    mesh = vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells([VTK_VERTEX], vertices)

    __build_global_ids(mesh, True, True)

    global_cell_ids = mesh.GetCellData().GetGlobalIds()
    global_point_ids = mesh.GetPointData().GetGlobalIds()
    assert global_cell_ids.GetNumberOfValues() == 1
    assert global_point_ids.GetNumberOfValues() == 1
