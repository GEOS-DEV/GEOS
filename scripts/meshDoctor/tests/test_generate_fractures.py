import numpy

import vtk  # TODO use new pyvtk style
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from checks.generate_fractures import Options, __split_mesh_on_fracture, __find_involved_cells, __color_fracture_sides


def __build_mesh() -> vtk.vtkUnstructuredGrid:
    # creation of a 3 x 3 x 1 grid
    rg = vtk.vtkRectilinearGrid()
    rg.SetDimensions(4, 4, 2)
    rg.SetXCoordinates(numpy_to_vtk(numpy.arange(4, dtype=float)))
    rg.SetYCoordinates(numpy_to_vtk(numpy.arange(4, dtype=float)))
    rg.SetZCoordinates(numpy_to_vtk(numpy.arange(2, dtype=float)))

    num_points = rg.GetNumberOfPoints()
    num_cells = rg.GetNumberOfCells()

    points = vtk.vtkPoints()
    points.Allocate(num_points)
    for i in range(num_points):
        points.InsertNextPoint(rg.GetPoint(i))

    cell_types = [vtk.VTK_HEXAHEDRON] * num_cells
    cells = vtk.vtkCellArray()
    cells.AllocateExact(num_cells, num_cells * 8)

    m = (0, 1, 3, 2, 4, 5, 7, 6)  # VTK_VOXEL and VTK_HEXAHEDRON do not share the same ordering.

    for i in range(rg.GetNumberOfCells()):
        c = rg.GetCell(i)
        new_cell = vtk.vtkHexahedron()
        for j in range(8):
            new_cell.GetPointIds().SetId(j, c.GetPointId(m[j]))
        cells.InsertNextCell(new_cell)

    mesh = vtk.vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells(cell_types, cells)

    return mesh


def test_generate_fracture():
    mesh = __build_mesh()

    ref = numpy.array((0, 1, 2, 1, 1, 2, 1, 1, 2), dtype=int)
    attribute = numpy_to_vtk(ref)
    attribute.SetName("attribute")
    mesh.GetCellData().AddArray(attribute)

    options = Options(policy="field", field="attribute", field_values={0, 1, 2}, output="")
    cell_frac_info, node_frac_info = __find_involved_cells(mesh, options)
    assert set(cell_frac_info.cell_to_faces.keys()) == {0, 1, 2, 3, 4, 5, 7, 8}
    assert len(cell_frac_info.field_data) == 5

    connected_cells = __color_fracture_sides(mesh, cell_frac_info, node_frac_info)
    cc = set(map(frozenset, connected_cells))
    assert cc == {frozenset({0}), frozenset({1, 3, 4, 7}), frozenset({2, 5, 8})}

    output_mesh = __split_mesh_on_fracture(mesh, options)
    assert mesh.GetNumberOfCells() == output_mesh.GetNumberOfCells()
    assert mesh.GetNumberOfPoints() + 6 + 8 == output_mesh.GetNumberOfPoints()

    cell_data = output_mesh.GetCellData()
    assert cell_data.HasArray("attribute")
    attribute = vtk_to_numpy(cell_data.GetArray("attribute"))
    assert not any(ref - attribute)
