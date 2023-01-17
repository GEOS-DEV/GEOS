import os
# import pathlib

import pytest

from vtkmodules.vtkCommonCore import (
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_POLYHEDRON,
    vtkPolyhedron,
    vtkUnstructuredGrid,
)

from checks.supported_elements import Options, check, __check
from checks.vtk_polyhedron import to_vtk_id_list


@pytest.mark.parametrize("base_name",
                         ("supportedElements.vtk", "supportedElementsAsVTKPolyhedra.vtk"))
def test_supported_elements(base_name):
    # p = pathlib.Path(".")
    # supported_elements_file_name = os.path.join(p.absolute(), "../../../../unitTests/meshTests", base_name)
    supported_elements_file_name = os.path.join("/docker-exchange", base_name)
    options = Options(chunk_size=1, num_proc=4)
    result = check(supported_elements_file_name, options)
    assert not result.unsupported_std_elements_types
    assert not result.unsupported_polyhedron_elements


# def make_dodecahedron():
#     dodecahedron = vtkPolyhedron()
#
#     for i in range(0, 20):
#         dodecahedron.GetPointIds().InsertNextId(i)
#
#     dodecahedron.GetPoints().InsertNextPoint(1.21412, 0, 1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(0.375185, 1.1547, 1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(-0.982247, 0.713644, 1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(-0.982247, -0.713644, 1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(0.375185, -1.1547, 1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(1.96449, 0, 0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(0.607062, 1.86835, 0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(-1.58931, 1.1547, 0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(-1.58931, -1.1547, 0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(0.607062, -1.86835, 0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(1.58931, 1.1547, -0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(-0.607062, 1.86835, -0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(-1.96449, 0, -0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(-0.607062, -1.86835, -0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(1.58931, -1.1547, -0.375185)
#     dodecahedron.GetPoints().InsertNextPoint(0.982247, 0.713644, -1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(-0.375185, 1.1547, -1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(-1.21412, 0, -1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(-0.375185, -1.1547, -1.58931)
#     dodecahedron.GetPoints().InsertNextPoint(0.982247, -0.713644, -1.58931)
#
#     faces = [12,  # number of faces
#              5, 0, 1, 2, 3, 4,  # number of ids on face, ids
#              5, 0, 5, 10, 6, 1,
#              5, 1, 6, 11, 7, 2,
#              5, 2, 7, 12, 8, 3,
#              5, 3, 8, 13, 9, 4,
#              5, 4, 9, 14, 5, 0,
#              5, 15, 10, 5, 14, 19,
#              5, 16, 11, 6, 10, 15,
#              5, 17, 12, 7, 11, 16,
#              5, 18, 13, 8, 12, 17,
#              5, 19, 14, 9, 13, 18,
#              5, 19, 18, 17, 16, 15]
#
#     dodecahedron.SetFaces(faces)
#     dodecahedron.Initialize()
#
#     return dodecahedron, faces


def make_dodecahedron():
    points = (
        (1.21412, 0, 1.58931),
        (0.375185, 1.1547, 1.58931),
        (-0.982247, 0.713644, 1.58931),
        (-0.982247, -0.713644, 1.58931),
        (0.375185, -1.1547, 1.58931),
        (1.96449, 0, 0.375185),
        (0.607062, 1.86835, 0.375185),
        (-1.58931, 1.1547, 0.375185),
        (-1.58931, -1.1547, 0.375185),
        (0.607062, -1.86835, 0.375185),
        (1.58931, 1.1547, -0.375185),
        (-0.607062, 1.86835, -0.375185),
        (-1.96449, 0, -0.375185),
        (-0.607062, -1.86835, -0.375185),
        (1.58931, -1.1547, -0.375185),
        (0.982247, 0.713644, -1.58931),
        (-0.375185, 1.1547, -1.58931),
        (-1.21412, 0, -1.58931),
        (-0.375185, -1.1547, -1.58931),
        (0.982247, -0.713644, -1.58931)
    )

    faces = (12,  # number of faces
             5, 0, 1, 2, 3, 4,  # number of ids on face, ids
             5, 0, 5, 10, 6, 1,
             5, 1, 6, 11, 7, 2,
             5, 2, 7, 12, 8, 3,
             5, 3, 8, 13, 9, 4,
             5, 4, 9, 14, 5, 0,
             5, 15, 10, 5, 14, 19,
             5, 16, 11, 6, 10, 15,
             5, 17, 12, 7, 11, 16,
             5, 18, 13, 8, 12, 17,
             5, 19, 14, 9, 13, 18,
             5, 19, 18, 17, 16, 15)

    p = vtkPoints()
    p.Allocate(len(points))
    for coords in points:
        p.InsertNextPoint(coords)

    f = to_vtk_id_list(faces)

    return p, f


def test_dodecahedron():
    # d, faces = make_dodecahedron()
    points, faces = make_dodecahedron()
    # cell_types = [VTK_POLYHEDRON]
    # cells = vtkCellArray()
    # cells.AllocateExact(1, d.GetPointIds().GetNumberOfIds())
    # cells.AllocateEstimate(1, 300)
    # cells.InsertNextCell(d)
    mesh = vtkUnstructuredGrid()
    mesh.Allocate(1)
    # mesh.SetPoints(d.GetPoints())
    mesh.SetPoints(points)
    # assert len(faces) == 73
    # mesh.InsertNextCell(VTK_POLYHEDRON, _mk_vtk_id_list(faces))
    # mesh.InsertNextCell(VTK_POLYHEDRON, to_vtk_id_list(faces))
    mesh.InsertNextCell(VTK_POLYHEDRON, faces)
    # mesh.SetCells(cell_types, cells)

    result = __check(mesh, Options(num_proc=1, chunk_size=1))
    assert set(result.unsupported_polyhedron_elements) == {0}
    assert set(result.unsupported_std_elements_types) == set()
