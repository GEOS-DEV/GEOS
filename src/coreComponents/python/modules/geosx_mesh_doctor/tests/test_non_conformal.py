import numpy

from checks.non_conformal import Options, __check
from checks.generate_cube import (
    build_rectilinear_blocks_mesh,
    XYZ,
)


def test_two_close_hexs():
    delta = 1.e-6
    tmp = numpy.arange(2, dtype=float)
    xyz0 = XYZ(tmp, tmp, tmp)
    xyz1 = XYZ(tmp + 1 + delta, tmp, tmp)
    mesh = build_rectilinear_blocks_mesh((xyz0, xyz1))

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
    mesh = build_rectilinear_blocks_mesh((xyz0, xyz1))

    options = Options(angle_tolerance=1., point_tolerance=delta / 2., face_tolerance=delta / 2.)

    results = __check(mesh, options)
    assert len(results.non_conformal_cells) == 0


def test_two_close_shifted_hexs():
    delta_x, delta_y = 1.e-6, 0.5
    tmp = numpy.arange(2, dtype=float)
    xyz0 = XYZ(tmp, tmp, tmp)
    xyz1 = XYZ(tmp + 1 + delta_x, tmp + delta_y, tmp + delta_y)
    mesh = build_rectilinear_blocks_mesh((xyz0, xyz1))

    options = Options(angle_tolerance=1., point_tolerance=delta_x * 2, face_tolerance=delta_x * 2)

    results = __check(mesh, options)
    assert len(results.non_conformal_cells) == 1
    assert set(results.non_conformal_cells[0]) == {0, 1}


def test_big_elem_next_to_small_elem():
    delta = 1.e-6
    tmp = numpy.arange(2, dtype=float)
    xyz0 = XYZ(tmp, tmp + 1, tmp + 1)
    xyz1 = XYZ(3 * tmp + 1 + delta, 3 * tmp, 3 * tmp)
    mesh = build_rectilinear_blocks_mesh((xyz0, xyz1))

    options = Options(angle_tolerance=1., point_tolerance=delta * 2, face_tolerance=delta * 2)

    results = __check(mesh, options)
    assert len(results.non_conformal_cells) == 1
    assert set(results.non_conformal_cells[0]) == {0, 1}
