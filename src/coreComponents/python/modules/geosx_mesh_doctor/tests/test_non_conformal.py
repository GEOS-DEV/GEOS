from dataclasses import dataclass
from typing import Iterable

import numpy

from vtkmodules.vtkCommonCore import (
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_HEXAHEDRON,
    vtkCellArray,
    vtkHexahedron,
    vtkRectilinearGrid,
    vtkUnstructuredGrid,
)
from vtk.util.numpy_support import (
    numpy_to_vtk,
    vtk_to_numpy,
)

from checks.non_conformal import Options, __check


@dataclass(frozen=True)
class XYZ:
    x: numpy.array
    y: numpy.array
    z: numpy.array


def __build_mesh(xyzs: Iterable[XYZ]):  # TODO move to testing utilities
    rgs = []
    for xyz in xyzs:
        rg = vtkRectilinearGrid()
        rg.SetDimensions(len(xyz.x), len(xyz.y), len(xyz.z))
        rg.SetXCoordinates(numpy_to_vtk(xyz.x))
        rg.SetYCoordinates(numpy_to_vtk(xyz.y))
        rg.SetZCoordinates(numpy_to_vtk(xyz.z))
        rgs.append(rg)

    num_points = sum(map(lambda r: r.GetNumberOfPoints(), rgs))
    num_cells = sum(map(lambda r: r.GetNumberOfCells(), rgs))

    points = vtkPoints()
    points.Allocate(num_points)
    for rg in rgs:
        for i in range(rg.GetNumberOfPoints()):
            points.InsertNextPoint(rg.GetPoint(i))

    cell_types = [VTK_HEXAHEDRON] * num_cells
    cells = vtkCellArray()
    cells.AllocateExact(num_cells, num_cells * 8)

    m = (0, 1, 3, 2, 4, 5, 7, 6)  # VTK_VOXEL and VTK_HEXAHEDRON do not share the same ordering.
    offset = 0
    for rg in rgs:
        for i in range(rg.GetNumberOfCells()):
            c = rg.GetCell(i)
            new_cell = vtkHexahedron()
            for j in range(8):
                new_cell.GetPointIds().SetId(j, offset + c.GetPointId(m[j]))
            cells.InsertNextCell(new_cell)
        offset += 8 * rg.GetNumberOfCells()

    mesh = vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells(cell_types, cells)

    return mesh


def test_two_close_hexs():
    delta = 1.e-6
    tmp = numpy.arange(2, dtype=float)
    xyz0 = XYZ(tmp, tmp, tmp)
    xyz1 = XYZ(tmp + 1 + delta, tmp, tmp)
    mesh = __build_mesh((xyz0, xyz1))

    # Close enough, but points tolerance is too strict to consider the faces matching.
    options = Options(angle_tolerance=1., point_tolerance=delta / 2, face_tolerance=delta * 2)
    results = __check(mesh, options)
    assert len(results.non_conformal_cells) == 1
    assert set(results.non_conformal_cells[0]) == {0, 1}

    # Close enough, and points tolerance is loose enough to consider the faces matching.
    options = Options(angle_tolerance=1., point_tolerance=delta * 2, face_tolerance=delta * 2)
    results = __check(mesh, options)
    assert len(results.non_conformal_cells) == 0


def test_two_distant_hexs():
    delta = 1
    tmp = numpy.arange(2, dtype=float)
    xyz0 = XYZ(tmp, tmp, tmp)
    xyz1 = XYZ(tmp + 1 + delta, tmp, tmp)
    mesh = __build_mesh((xyz0, xyz1))

    options = Options(angle_tolerance=1., point_tolerance=delta / 2., face_tolerance=delta / 2.)

    results = __check(mesh, options)
    assert len(results.non_conformal_cells) == 0


def test_two_close_shifted_hexs():
    delta_x, delta_y = 1.e-6, 0.5
    tmp = numpy.arange(2, dtype=float)
    xyz0 = XYZ(tmp, tmp, tmp)
    xyz1 = XYZ(tmp + 1 + delta_x, tmp + delta_y, tmp + delta_y)
    mesh = __build_mesh((xyz0, xyz1))

    options = Options(angle_tolerance=1., point_tolerance=delta_x * 2, face_tolerance=delta_x * 2)

    results = __check(mesh, options)
    assert len(results.non_conformal_cells) == 1
    assert set(results.non_conformal_cells[0]) == {0, 1}


def test_big_elem_next_to_small_elem():
    delta = 1.e-6
    tmp = numpy.arange(2, dtype=float)
    xyz0 = XYZ(tmp, tmp + 1, tmp + 1)
    xyz1 = XYZ(3 * tmp + 1 + delta, 3 * tmp, 3 * tmp)
    mesh = __build_mesh((xyz0, xyz1))

    options = Options(angle_tolerance=1., point_tolerance=delta * 2, face_tolerance=delta * 2)

    results = __check(mesh, options)
    assert len(results.non_conformal_cells) == 1
    assert set(results.non_conformal_cells[0]) == {0, 1}
