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
)


@dataclass(frozen=True)
class XYZ:
    x: numpy.array
    y: numpy.array
    z: numpy.array


def build_rectilinear_blocks_mesh(xyzs: Iterable[XYZ]):
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
