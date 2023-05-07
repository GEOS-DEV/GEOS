from dataclasses import dataclass
from typing import Generator

import pytest

from vtkmodules.vtkCommonCore import (
    vtkIdList,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_POLYHEDRON,
    vtkUnstructuredGrid,
)

import numpy

from checks.reorient_mesh import reorient_mesh
from checks.vtk_polyhedron import FaceStream
from checks.vtk_utils import (
    to_vtk_id_list,
    vtk_iter,
)


@dataclass(frozen=True)
class Expected:
    mesh: vtkUnstructuredGrid
    face_stream: FaceStream


def __build_test_meshes() -> Generator[Expected, None, None]:
    # Creating the support nodes for the polyhedron.
    # It has a C shape and is actually non-convex, non star-shaped.
    front_nodes = numpy.array((
        (0, 0, 0),
        (3, 0, 0),
        (3, 1, 0),
        (1, 1, 0),
        (1, 2, 0),
        (3, 2, 0),
        (3, 3, 0),
        (0, 3, 0),
    ), dtype=float)
    front_nodes = numpy.array(front_nodes, dtype=float)
    back_nodes = front_nodes - (0., 0., 1.)

    n = len(front_nodes)

    points = vtkPoints()
    points.Allocate(2 * n)
    for coords in front_nodes:
        points.InsertNextPoint(coords)
    for coords in back_nodes:
        points.InsertNextPoint(coords)

    # Creating the polyhedron with faces all directed outward.
    faces = []
    # Creating the side faces
    for i in range(n):
        faces.append(
            (i % n + n, (i + 1) % n + n, (i + 1) % n, i % n)
        )
    # Creating the front faces
    faces.append(tuple(range(n)))
    faces.append(tuple(reversed(range(n, 2 * n))))
    face_stream = FaceStream(faces)

    # Creating multiple meshes, each time with one unique polyhedron,
    # but with different "face flip status".
    # First case, no face is flipped.
    mesh = vtkUnstructuredGrid()
    mesh.Allocate(1)
    mesh.SetPoints(points)
    mesh.InsertNextCell(VTK_POLYHEDRON, to_vtk_id_list(
        face_stream.dump()
    ))
    yield Expected(mesh=mesh, face_stream=face_stream)

    # Here, two faces are flipped.
    mesh = vtkUnstructuredGrid()
    mesh.Allocate(1)
    mesh.SetPoints(points)
    mesh.InsertNextCell(VTK_POLYHEDRON, to_vtk_id_list(
        face_stream.flip_faces((1, 2)).dump()
    ))
    yield Expected(mesh=mesh, face_stream=face_stream)

    # Last, all faces are flipped.
    mesh = vtkUnstructuredGrid()
    mesh.Allocate(1)
    mesh.SetPoints(points)
    mesh.InsertNextCell(VTK_POLYHEDRON, to_vtk_id_list(
        face_stream.flip_faces(range(len(faces))).dump()
    ))
    yield Expected(mesh=mesh, face_stream=face_stream)


@pytest.mark.parametrize("expected", __build_test_meshes())
def test_reorient_polyhedron(expected: Expected):
    output_mesh = reorient_mesh(expected.mesh, range(expected.mesh.GetNumberOfCells()))
    assert output_mesh.GetNumberOfCells() == 1
    assert output_mesh.GetCell(0).GetCellType() == VTK_POLYHEDRON
    face_stream_ids = vtkIdList()
    output_mesh.GetFaceStream(0, face_stream_ids)
    # Note that the following makes a raw (but simple) check.
    # But one may need to be more precise some day,
    # since triangular faces (0, 1, 2) and (1, 2, 0) should be considered as equivalent.
    # And the current simpler check does not consider this case.
    assert tuple(vtk_iter(face_stream_ids)) == expected.face_stream.dump()
