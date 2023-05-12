from checks.generate_cube import __build, Options, FieldInfo


def test_generate_cube():
    options = Options(
        vtk_output=None,
        generate_cells_global_ids=True,
        generate_points_global_ids=False,
        xs=(0, 5, 10),
        ys=(0, 4, 8),
        zs=(0, 1),
        nxs=(5, 2),
        nys=(1, 1),
        nzs=(1,),
        fields=(
            FieldInfo(name="test", dimension=2, support="CELLS"),
        )
    )
    output = __build(options)
    assert output.GetNumberOfCells() == 14
    assert output.GetNumberOfPoints() == 48
    assert output.GetCellData().GetArray("test").GetNumberOfComponents() == 2
    assert output.GetCellData().GetGlobalIds()
    assert not output.GetPointData().GetGlobalIds()
