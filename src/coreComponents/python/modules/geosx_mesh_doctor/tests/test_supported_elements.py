import os
from typing import Tuple

import pytest

from vtkmodules.vtkCommonCore import (
    vtkIdList,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_POLYHEDRON,
    vtkUnstructuredGrid,
)

from checks.supported_elements import Options, check, __check
from checks.vtk_polyhedron import parse_face_stream, build_face_to_face_connectivity_through_edges, FaceStream
from checks.vtk_utils import (
    to_vtk_id_list,
)


@pytest.mark.parametrize("base_name",
                         ("supportedElements.vtk", "supportedElementsAsVTKPolyhedra.vtk"))
def test_supported_elements(base_name) -> None:
    """
    Testing that the supported elements are properly detected as supported!
    :param base_name: Supported elements are provided as standard elements or polyhedron elements.
    """
    directory = os.path.dirname(os.path.realpath(__file__))
    supported_elements_file_name = os.path.join(directory, "../../../../unitTests/meshTests", base_name)
    options = Options(chunk_size=1, num_proc=4)
    result = check(supported_elements_file_name, options)
    assert not result.unsupported_std_elements_types
    assert not result.unsupported_polyhedron_elements


def make_dodecahedron() -> Tuple[vtkPoints, vtkIdList]:
    """
    Returns the points and faces for a dodecahedron.
    This code was adapted from an official vtk example.
    :return: The tuple of points and faces (as vtk instances).
    """
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


def test_dodecahedron() -> None:
    """
    Tests that a dodecahedron is not supported by GEOSX.
    """
    points, faces = make_dodecahedron()
    mesh = vtkUnstructuredGrid()
    mesh.Allocate(1)
    mesh.SetPoints(points)
    mesh.InsertNextCell(VTK_POLYHEDRON, faces)

    result = __check(mesh, Options(num_proc=1, chunk_size=1))
    assert set(result.unsupported_polyhedron_elements) == {0}
    assert not result.unsupported_std_elements_types


def test_parse_face_stream() -> None:
    _, faces = make_dodecahedron()
    result = parse_face_stream(faces)
    expected = (
        (0, 1, 2, 3, 4),
        (0, 5, 10, 6, 1),
        (1, 6, 11, 7, 2),
        (2, 7, 12, 8, 3),
        (3, 8, 13, 9, 4),
        (4, 9, 14, 5, 0),
        (15, 10, 5, 14, 19),
        (16, 11, 6, 10, 15),
        (17, 12, 7, 11, 16),
        (18, 13, 8, 12, 17),
        (19, 14, 9, 13, 18),
        (19, 18, 17, 16, 15)
    )
    assert result == expected
    face_stream = FaceStream.build_from_vtk_id_list(faces)
    assert face_stream.num_faces == 12
    assert face_stream.num_support_points == 20
