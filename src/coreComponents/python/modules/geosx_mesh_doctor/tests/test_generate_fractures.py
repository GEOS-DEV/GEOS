import numpy

from vtkmodules.util.numpy_support import (
    numpy_to_vtk,
    vtk_to_numpy,
)

import sys
sys.path.insert(0, "/Users/j0436735/CLionProjects/GEOS/src/coreComponents/python/modules/geosx_mesh_doctor")

from checks.generate_fractures import Options, __split_mesh_on_fracture, __color_fracture_sides, vtk_utils, FractureInfo
from checks.generate_cube import (
    build_rectilinear_blocks_mesh,
    XYZ,
)


def test_generate_fracture():
    tmp0 = numpy.arange(4, dtype=float)
    tmp1 = numpy.arange(2, dtype=float)
    xyz = XYZ(tmp0, tmp0, tmp1)
    mesh = build_rectilinear_blocks_mesh((xyz, ))

    ref = numpy.array((0, 1, 2, 1, 1, 2, 1, 1, 2), dtype=int)
    attribute = numpy_to_vtk(ref)
    attribute.SetName("attribute")
    mesh.GetCellData().AddArray(attribute)

    options = Options(policy="field",
                      field="attribute",
                      field_type="cells",
                      field_values={0, 1, 2},
                      vtk_output=None,
                      vtk_fracture_output=None,
                      split_on_domain_boundary=True)
    frac = FractureInfo(mesh, options)
    assert set(frac.cell_to_faces.keys()) == {0, 1, 2, 3, 4, 5, 7, 8}
    assert numpy.array_equal(frac.is_internal_fracture_node,
                             numpy.array(([0, 1, 1, 0, 1, 1, 1, 0] + [0, 0, 1, 0] * 2) * 2, dtype=bool))

    connected_cells = __color_fracture_sides(mesh,
                                             frac.cell_to_faces,
                                             frac.is_fracture_node)
    cc = set(map(frozenset, connected_cells))
    assert cc == {frozenset({0}), frozenset({1, 3, 4, 7}), frozenset({2, 5, 8})}

    output_mesh, fracture_mesh = __split_mesh_on_fracture(mesh, options)
    assert mesh.GetNumberOfCells() == output_mesh.GetNumberOfCells()
    assert mesh.GetNumberOfPoints() + 6 + 8 == output_mesh.GetNumberOfPoints()
    assert fracture_mesh.GetNumberOfCells() == 5
    assert fracture_mesh.GetNumberOfPoints() == 6 + 8
    dup = fracture_mesh.GetPointData().GetArray("duplicated_nodes")
    assert dup.GetNumberOfComponents() == 2
    assert dup.GetNumberOfTuples() == fracture_mesh.GetNumberOfPoints()

    cell_data = output_mesh.GetCellData()
    assert cell_data.HasArray("attribute")
    attribute = vtk_to_numpy(cell_data.GetArray("attribute"))
    assert not any(ref - attribute)


def test_generate_fracture2():
    tmp0 = numpy.arange(3, dtype=float)
    tmp1 = numpy.arange(2, dtype=float)
    xyz = XYZ(tmp0, tmp0, tmp1)
    mesh = build_rectilinear_blocks_mesh((xyz, ))

    ref = numpy.array((0, 1, 1, 1), dtype=int)
    attribute = numpy_to_vtk(ref)
    attribute.SetName("attribute")
    mesh.GetCellData().AddArray(attribute)

    options = Options(policy="field",
                      field="attribute",
                      field_type="cells",
                      field_values={0, 1},
                      vtk_output=None,
                      vtk_fracture_output=None,
                      split_on_domain_boundary=True)
    frac = FractureInfo(mesh, options)
    # assert set(frac.cell_to_faces.keys()) == {0, 1, 2, 3, 4, 5, 7, 8}
    # assert numpy.array_equal(frac.is_internal_fracture_node,
    #                          numpy.array(([0, 1, 1, 0, 1, 1, 1, 0] + [0, 0, 1, 0] * 2) * 2, dtype=bool))

    connected_cells = __color_fracture_sides(mesh,
                                             frac.cell_to_faces,
                                             frac.is_fracture_node)
    cc = set(map(frozenset, connected_cells))
    assert cc == {frozenset({0}), frozenset({1, 2, 3})}

    # output_mesh, fracture_mesh = __split_mesh_on_fracture(mesh, options)
    # assert mesh.GetNumberOfCells() == output_mesh.GetNumberOfCells()
    # assert mesh.GetNumberOfPoints() + 6 + 8 == output_mesh.GetNumberOfPoints()
    # assert fracture_mesh.GetNumberOfCells() == 5
    # assert fracture_mesh.GetNumberOfPoints() == 6 + 8
    # dup = fracture_mesh.GetPointData().GetArray("duplicated_nodes")
    # assert dup.GetNumberOfComponents() == 2
    # assert dup.GetNumberOfTuples() == fracture_mesh.GetNumberOfPoints()
    #
    # cell_data = output_mesh.GetCellData()
    # assert cell_data.HasArray("attribute")
    # attribute = vtk_to_numpy(cell_data.GetArray("attribute"))
    # assert not any(ref - attribute)


if __name__ == '__main__':
    test_generate_fracture()
    test_generate_fracture2()
