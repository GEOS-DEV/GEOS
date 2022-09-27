# import logging
import numpy

import vtk

# from collections import defaultdict
from vtk.util.numpy_support import numpy_to_vtk


def main():
    r0 = vtk.vtkRectilinearGrid()
    r0.SetDimensions(5 + 1, 10 + 1, 1 + 1)
    r0.SetXCoordinates(numpy_to_vtk(numpy.ones(-5, 0 + 1, 1, dtype=float)))
    r0.SetYCoordinates(numpy_to_vtk(numpy.arange( 0, 10 + 1, 1, dtype=float)))
    r0.SetZCoordinates(numpy_to_vtk(numpy.arange( 0, 1 + 1, 1, dtype=float)))

    r1 = vtk.vtkRectilinearGrid()
    r1.SetDimensions(5 + 1, 5 + 1, 1 + 1)
    r1.SetXCoordinates(numpy_to_vtk(numpy.arange( 0, 5 + 1, 1, dtype=float)))
    r1.SetYCoordinates(numpy_to_vtk(numpy.arange( 0, 5 + 1, 1, dtype=float)))
    r1.SetZCoordinates(numpy_to_vtk(numpy.arange( 0, 1 + 1, 1, dtype=float)))

    r2 = vtk.vtkRectilinearGrid()
    r2.SetDimensions(5 + 1, 5 + 1, 1 + 1)
    r2.SetXCoordinates(numpy_to_vtk(numpy.arange( 0, 5 + 1, 1, dtype=float)))
    r2.SetYCoordinates(numpy_to_vtk(numpy.arange( 5, 10 + 1, 1, dtype=float)))
    r2.SetZCoordinates(numpy_to_vtk(numpy.arange( 0, 1 + 1, 1, dtype=float)))

    rs = r0, r1, r2

    # num_points_0 = r0.GetNumberOfPoints()
    # num_points_1 = r1.GetNumberOfPoints()
    # num_points_2 = r2.GetNumberOfPoints()
    # num_cells_0 = r0.GetNumberOfCells()
    # num_cells_1 = r1.GetNumberOfCells()
    # num_cells_2 = r2.GetNumberOfCells()

    num_points, num_cells = 0, 0
    for r in rs:
        num_points += r.GetNumberOfPoints()
        num_cells += r.GetNumberOfCells()

    points = vtk.vtkPoints()
    points.Allocate(num_points)
    # points.Allocate(num_points_0 + num_points_1 + num_points_2)

    # for i in range(num_points_0):
    #     points.InsertNextPoint(r0.GetPoint(i))
    # for i in range(num_points_1):
    #     points.InsertNextPoint(r1.GetPoint(i))
    # for i in range(num_points_2):
    #     points.InsertNextPoint(r2.GetPoint(i))
    for r in rs:
        for i in range(r.GetNumberOfPoints()):
            points.InsertNextPoint(r.GetPoint(i))


    # num_cells = num_cells_0 + num_cells_1 + num_cells_2
    cell_types = [vtk.VTK_HEXAHEDRON] * num_cells
    cells = vtk.vtkCellArray()
    cells.AllocateExact(num_cells, num_cells * 8)

    m = (0, 1, 3, 2, 4, 5, 7, 6)  # VTK_VOXEL and VTK_HEXAHEDRON do not share the same ordering.

    points_offset = 0
    for r in rs:
        for i in range(r.GetNumberOfCells()):
            c = r.GetCell(i)
            new_cell = vtk.vtkHexahedron()
            for j in range(8):
                new_cell.GetPointIds().SetId(j, c.GetPointId(m[j]) + points_offset)
            cells.InsertNextCell(new_cell)
        points_offset += r.GetNumberOfPoints()

    # for i in range(num_cells_0):
    #     c = r0.GetCell(i)
    #     new_cell = vtk.vtkHexahedron()
    #     for j in range(8):
    #         new_cell.GetPointIds().SetId(j, c.GetPointId(m[j]))
    #     cells.InsertNextCell(new_cell)
    # for i in range(num_cells_1):
    #     c = r1.GetCell(i)
    #     new_cell = vtk.vtkHexahedron()
    #     for j in range(8):
    #         new_cell.GetPointIds().SetId(j, c.GetPointId(m[j]) + num_points_0)
    #     cells.InsertNextCell(new_cell)
    # for i in range(num_cells_2):
    #     c = r2.GetCell(i)
    #     new_cell = vtk.vtkHexahedron()
    #     for j in range(8):
    #         new_cell.GetPointIds().SetId(j, c.GetPointId(m[j]) + num_points_0 + num_points_1)
    #     cells.InsertNextCell(new_cell)

    mesh = vtk.vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells(cell_types, cells)

    # For points...
    point_global_ids = vtk.vtkIdTypeArray()
    point_global_ids.SetName("GLOBAL_IDS_POINTS")
    point_global_ids.Allocate(mesh.GetNumberOfPoints())
    for i in range(mesh.GetNumberOfPoints()):
        point_global_ids.InsertNextValue(i)
    mesh.GetPointData().SetGlobalIds(point_global_ids)
    # ... then for cells.
    cells_global_ids = vtk.vtkIdTypeArray()
    cells_global_ids.SetName("GLOBAL_IDS_CELLS")
    cells_global_ids.Allocate(mesh.GetNumberOfCells())
    for i in range(mesh.GetNumberOfCells()):
        cells_global_ids.InsertNextValue(i)
    mesh.GetCellData().SetGlobalIds(cells_global_ids)

    data = (
        (4, 1, 50, 0),
        (9, 1, 55, 0),
        (14, 1, 60, 0),
        (19, 1, 65, 0),
        (24, 1, 70, 0),
        (29, 1, 75, 0),
        (34, 1, 80, 0),
        (39, 1, 85, 0),
        (44, 1, 90, 0),
        (49, 1, 95, 0),
        (70, 3, 75, 2),
        (71, 3, 76, 2),
        (72, 3, 77, 2),
        (73, 3, 78, 2),
        (74, 3, 79, 2)
    )
    assert len(data) == 15

    da = vtk.vtkIntArray()
    da.SetName("fracture_info")
    da.SetNumberOfComponents(4)
    da.SetNumberOfTuples(15)
    for i in range(da.GetNumberOfTuples()):
        # e0, f0, e1, f1 = data[i]
        da.SetTuple(i, data[i])
    fracture_field_data = vtk.vtkFieldData()
    fracture_field_data.AddArray(da)

    mesh.SetFieldData(fracture_field_data)

    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileTypeToASCII()
    writer.SetFileName("/docker-exchange/frac-t.vtk")
    writer.SetInputData(mesh)
    writer.Write()


if __name__ == '__main__':
    main()
