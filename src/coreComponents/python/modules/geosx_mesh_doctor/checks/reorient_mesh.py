import logging
from typing import (
    Dict,
    FrozenSet,
    Iterator,
    List,
    Tuple,
)

import numpy

from tqdm import tqdm

import networkx

from vtkmodules.vtkCommonCore import (
    vtkIdList,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_POLYHEDRON,
    VTK_TRIANGLE,
    vtkCellArray,
    vtkPolyData,
    vtkPolygon,
    vtkUnstructuredGrid,
    vtkTetra,
)
from vtkmodules.vtkFiltersCore import (
    vtkTriangleFilter,
)
from .vtk_utils import (
    to_vtk_id_list,
)

from .vtk_polyhedron import (
    FaceStream,
    build_face_to_face_connectivity_through_edges,
)


def __compute_volume(mesh_points: vtkPoints, face_stream: FaceStream) -> float:
    """
    Computes the volume of a polyhedron element (defined by its face_stream).
    :param mesh_points: The mesh points, needed to compute the volume.
    :param face_stream: The vtk face stream.
    :return: The volume of the element.
    :note: The faces of the polyhedron are triangulated and the volumes of the tetrahedra
    from the barycenter to the triangular bases are summed.
    The normal of each face plays critical role,
    since the volume of each tetrahedron can be positive or negative.
    """
    # Triangulating the envelope of the polyhedron for further volume computation.
    polygons = vtkCellArray()
    for face_nodes in face_stream.face_nodes:
        polygon = vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(len(face_nodes))
        # We use the same global points numbering for the polygons than for the input mesh.
        # There will be a lot of points in the poly data that won't be used as a support for the polygons.
        # But the algorithm deals with it, and it's actually faster (and easier) to do this
        # than to renumber and allocate a new fit-for-purpose set of points just for the polygons.
        for i, point_id in enumerate(face_nodes):
            polygon.GetPointIds().SetId(i, point_id)
        polygons.InsertNextCell(polygon)
    polygon_poly_data = vtkPolyData()
    polygon_poly_data.SetPoints(mesh_points)
    polygon_poly_data.SetPolys(polygons)

    f = vtkTriangleFilter()
    f.SetInputData(polygon_poly_data)
    f.Update()
    triangles = f.GetOutput()
    # Computing the barycenter that will be used as the tip of all the tetra which mesh the polyhedron.
    # (The basis of all the tetra being the triangles of the envelope).
    # We could take any point, not only the barycenter.
    # But in order to work with figure of the same magnitude, let's compute the barycenter.
    tmp_barycenter = numpy.empty((face_stream.num_support_points, 3), dtype=float)
    for i, point_id in enumerate(face_stream.support_point_ids):
        tmp_barycenter[i, :] = mesh_points.GetPoint(point_id)
    barycenter = tmp_barycenter[:, 0].mean(), tmp_barycenter[:, 1].mean(), tmp_barycenter[:, 2].mean()
    # Looping on all the triangles of the envelope of the polyhedron, creating the matching tetra.
    # Then the volume of all the tetra are added to get the final polyhedron volume.
    cell_volume = 0.
    for i in range(triangles.GetNumberOfCells()):
        triangle = triangles.GetCell(i)
        assert triangle.GetCellType() == VTK_TRIANGLE
        p = triangle.GetPoints()
        cell_volume += vtkTetra.ComputeVolume(barycenter, p.GetPoint(0), p.GetPoint(1), p.GetPoint(2))
    return cell_volume


def __select_and_flip_faces(mesh_points: vtkPoints,
                            colors: Dict[FrozenSet[int], int],
                            face_stream: FaceStream) -> FaceStream:
    """
    Given a polyhedra, given that we were able to paint the faces in two colors,
    we now need to select which faces/color to flip such that the volume of the element is positive.
    :param mesh_points: The mesh points, needed to compute the volume.
    :param colors: Maps the nodes of each connected component (defined as a frozenset) to its color.
    :param face_stream: the polyhedron.
    :return: The face stream that leads to a positive volume.
    """
    # Flipping either color 0 or 1.
    color_to_nodes: Dict[int, List[int]] = {0: [], 1: []}
    for connected_components_indices, color in colors.items():
        color_to_nodes[color] += connected_components_indices
    # This implementation works even if there is one unique color.
    # Admittedly, there will be one face stream that won't be flipped.
    fs: Tuple[FaceStream, FaceStream] = face_stream.flip_faces(color_to_nodes[0]), face_stream.flip_faces(color_to_nodes[1])
    volumes = __compute_volume(mesh_points, fs[0]), __compute_volume(mesh_points, fs[1])
    # We keep the flipped element for which the volume is largest
    # (i.e. positive, since they should be the opposite of each other).
    return fs[numpy.argmax(volumes)]


def __reorient_element(mesh_points: vtkPoints, face_stream_ids: vtkIdList) -> vtkIdList:
    """
    Considers a vtk face stream and flips the appropriate faces to get an element with normals directed outwards.
    :param mesh_points: The mesh points, needed to compute the volume.
    :param face_stream_ids: The raw vtk face stream, not converted into a more practical python class.
    :return: The raw vtk face stream with faces properly flipped.
    """
    face_stream = FaceStream.build_from_vtk_id_list(face_stream_ids)
    face_graph = build_face_to_face_connectivity_through_edges(face_stream, add_compatibility=True)
    # Removing the non-compatible connections to build the non-connected components.
    g = networkx.Graph()
    g.add_nodes_from(face_graph.nodes)
    g.add_edges_from(filter(lambda uvd: uvd[2]["compatible"] == "+", face_graph.edges(data=True)))
    connected_components = tuple(networkx.connected_components(g))
    # Squashing all the connected nodes that need to receive the normal direction flip (or not) together.
    quotient_graph = networkx.algorithms.quotient_graph(face_graph, connected_components)
    # Coloring the new graph lets us know how which cluster of faces need to eventually receive the same flip.
    # W.r.t. the nature of our problem (a normal can be directed inwards or outwards),
    # two colors should be enough to color the face graph.
    # `colors` maps the nodes of each connected component to its color.
    colors: Dict[FrozenSet[int], int] = networkx.algorithms.greedy_color(quotient_graph)
    assert len(colors) in (1, 2)
    # We now compute the face stream which generates outwards normal vectors.
    flipped_face_stream = __select_and_flip_faces(mesh_points, colors, face_stream)
    return to_vtk_id_list(flipped_face_stream.dump())


def reorient_mesh(mesh, cell_indices: Iterator[int]) -> vtkUnstructuredGrid:
    """
    Reorient the polyhedron elements such that they all have their normals directed outwards.
    :param mesh: The input vtk mesh.
    :param cell_indices: We may need to only flip a limited number of polyhedron cells (only on the boundary for example).
    :return: The vtk mesh with the desired polyhedron cells directed outwards.
    """
    num_cells = mesh.GetNumberOfCells()
    # Building an indicator/predicate from the list
    needs_to_be_reoriented = numpy.zeros(num_cells, dtype=bool)
    for ic in cell_indices:
        needs_to_be_reoriented[ic] = True

    output_mesh = mesh.NewInstance()
    # I did not manage to call `output_mesh.CopyStructure(mesh)` because I could not modify the polyhedron in place.
    # Therefore, I insert the cells one by one...
    output_mesh.SetPoints(mesh.GetPoints())
    logging.info("Reorienting the polyhedron cells to enforce normals directed outward.")
    with tqdm(total=needs_to_be_reoriented.sum(), desc="Reorienting polyhedra") as progress_bar:  # For smoother progress, we only update on reoriented elements.
        for ic in range(num_cells):
            cell = mesh.GetCell(ic)
            cell_type = cell.GetCellType()
            if cell_type == VTK_POLYHEDRON:
                face_stream_ids = vtkIdList()
                mesh.GetFaceStream(ic, face_stream_ids)
                if needs_to_be_reoriented[ic]:
                    new_face_stream_ids = __reorient_element(mesh.GetPoints(), face_stream_ids)
                else:
                    new_face_stream_ids = face_stream_ids
                output_mesh.InsertNextCell(VTK_POLYHEDRON, new_face_stream_ids)
            else:
                output_mesh.InsertNextCell(cell_type, cell.GetPointIds())
            if needs_to_be_reoriented[ic]:
                progress_bar.update(1)
    assert output_mesh.GetNumberOfCells() == mesh.GetNumberOfCells()
    return output_mesh
