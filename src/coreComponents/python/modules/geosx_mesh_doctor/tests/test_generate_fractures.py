import numpy

from vtk.util.numpy_support import (
    numpy_to_vtk,
    vtk_to_numpy,
)

from checks.generate_fractures import Options, __split_mesh_on_fracture, __find_involved_cells, __color_fracture_sides
from .test_utils import (
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
                      field_values={0, 1, 2},
                      output="",
                      split_on_domain_boundary=True)
    cell_frac_info, node_frac_info = __find_involved_cells(mesh, options)
    assert set(cell_frac_info.cell_to_faces.keys()) == {0, 1, 2, 3, 4, 5, 7, 8}
    assert len(cell_frac_info.field_data) == 5
    assert numpy.array_equal(node_frac_info.is_internal_fracture_node,
                             numpy.array(([0, 1, 1, 0, 1, 1, 1, 0] + [0, 0, 1, 0] * 2) * 2, dtype=bool))

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
