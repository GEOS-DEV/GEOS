from dataclasses import dataclass
import itertools
import logging
import math
from typing import List, Tuple, Dict, FrozenSet
import numpy
import scipy.optimize

from tqdm import tqdm

import networkx

from vtkmodules.vtkCommonCore import (
    vtkIdList,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_POLYHEDRON,
    VTK_TRIANGLE,
    vtkBoundingBox,
    vtkCellArray,
    vtkPolyData,
    vtkPolygon,
    vtkPolyhedron,
    vtkUnstructuredGrid,
    vtkStaticPointLocator,
    vtkTetra,
)
from vtkmodules.vtkFiltersCore import (
    vtkCellCenters,
    vtkPolyDataNormals,
)
from vtkmodules.vtkFiltersGeometry import (
    vtkMarkBoundaryFilter,
    vtkDataSetSurfaceFilter,
)
# from vtk.util.numpy_support import (
#     vtk_to_numpy, )
import vtk  # TODO Fix for reference

# from . import vtk_utils
import vtk_utils

# from .vtk_polyhedron import (
from vtk_polyhedron import (  # TODO rename to vtk_polyhedron_utils
    FaceStream,
    build_cell_graph,
    to_vtk_id_list,
)


@dataclass(frozen=True)
class Options:
    angle_tolerance: float
    point_tolerance: float
    face_tolerance: float


@dataclass(frozen=True)
class Result:
    non_conformal_cells: List[Tuple[int, int]]


def get_cell_field_by_name(mesh, field_name):
    cd = mesh.GetCellData()
    for i in range(cd.GetNumberOfArrays()):
        if cd.GetArrayName(i) == field_name:
            return cd.GetArray(i)


def distance_between_cells(x: Tuple[float, float, float, float], i: int, j: int, boundary_mesh: vtkPolyData) -> float:
    assert len(x) == 4
    xi, yi, xj, yj = x
    sub_id = vtk.reference(0)
    # First cell
    ci = boundary_mesh.GetCell(i)
    assert ci.IsPrimaryCell()
    pi = numpy.empty(3)
    wi = numpy.empty(ci.GetNumberOfPoints())
    ci.EvaluateLocation(sub_id, (xi, yi, 0), pi, wi)
    # Second cell
    cj = boundary_mesh.GetCell(j)
    assert cj.IsPrimaryCell()
    pj = numpy.empty(3)
    wj = numpy.empty(cj.GetNumberOfPoints())
    cj.EvaluateLocation(sub_id, (xj, yj, 0), pj, wj)
    return numpy.linalg.norm(pj - pi)


def __build_boundary_mesh(mesh, options: Options) -> Result:
    ORIGINAL_POINTS = "ORIGINAL_POINTS"
    ORIGINAL_CELLS = "ORIGINAL_CELLS"

    f = vtkDataSetSurfaceFilter()
    f.PassThroughCellIdsOn()
    f.PassThroughPointIdsOn()
    f.FastModeOff()

    f.SetOriginalCellIdsName(ORIGINAL_CELLS)
    f.SetOriginalPointIdsName(ORIGINAL_POINTS)

    boundary_mesh = vtkPolyData()
    f.UnstructuredGridExecute(mesh, boundary_mesh)

    cc = vtkCellCenters()
    cc.VertexCellsOn()
    cc.CopyArraysOff()
    cc.SetInputData(boundary_mesh)
    cc.Update()
    cell_centers_mesh = cc.GetOutput()
    cell_centers: List[Tuple[float, float, float]] = []
    for ic in range(cell_centers_mesh.GetNumberOfCells()):
        cell_centers.append(numpy.array(cell_centers_mesh.GetPoint(cell_centers_mesh.GetCell(ic).GetPointId(0))))

    n = vtkPolyDataNormals()
    n.ConsistencyOff()  # TODO review the options!
    n.FlipNormalsOn()
    n.AutoOrientNormalsOff()
    n.ComputeCellNormalsOn()
    n.SetInputData(boundary_mesh)
    n.Update()
    normals = get_cell_field_by_name(n.GetOutput(), "Normals")
    assert normals
    assert normals.GetName() == "Normals"
    assert normals.GetNumberOfComponents() == 3
    assert normals.GetNumberOfTuples() == boundary_mesh.GetNumberOfCells()

    return boundary_mesh, normals, cell_centers


def __check(mesh, options: Options) -> Result:
    boundary_mesh, normals, cell_centers = __build_boundary_mesh(mesh, options)
    cos_theta = abs(math.cos(numpy.deg2rad(options.angle_tolerance)))

    non_conformal_cells = []

    num_cells = boundary_mesh.GetNumberOfCells()

    # Precomputing the bounding boxes of all the boundary cells.
    bounding_boxes = []
    for i in range(num_cells):
        bb = vtkBoundingBox(boundary_mesh.GetCell(i).GetBounds())
        bb.Inflate(options.face_tolerance)
        bounding_boxes.append(bb)

    # Looping on all the pairs of boundary cells. We'll hopefully discard most of the pairs.
    for i, j in tqdm(itertools.combinations(range(num_cells), r=2),
                     desc="Non conformal elements",
                     total=(num_cells - 1) * num_cells // 2):
        # Discarding pairs that are obviously too far from each other.
        if not bounding_boxes[i].Intersects(bounding_boxes[j]):
            continue
        # Discarding pairs that are not facing each others (with a threshold).
        normal_i, normal_j = normals.GetTuple3(i), normals.GetTuple3(j)
        if numpy.dot(normal_i, normal_j) > -cos_theta:  # opposite directions only (can be facing or not)
            continue
        # cci, ccj = cell_centers[i], cell_centers[j]
        # direction_i = numpy.dot(ccj - cci, normal_i)  # TODO use closest point to polygon instead (for ccj)
        # direction_j = numpy.dot(cci - ccj, normal_j)  # TODO use
        # # print(f"{i}, {j}, {scalar_product}, {direction_i}, {direction_j}")
        # if direction_i < 0 or direction_j < 0:  # checking direction now.
        #     continue
        logging.debug(f"Computing the distance for faces {i} and {j}.")
        result = scipy.optimize.minimize(lambda x: distance_between_cells(x, i, j, boundary_mesh),  # TODO that somehow sucks
                                         (0.5, 0.5, 0.5, 0.5),  # TODO do a better screening!
                                         method="Nelder-Mead",
                                         bounds=((0, 1), (0, 1), (0, 1), (0, 1))
                                         )
        # TODO deal with error code.
        logging.debug(f"Result between {i} and {j} is {result.fun}")
        if abs(result.fun) < options.face_tolerance:
            non_conformal_cells.append((i, j))

    original_cells = get_cell_field_by_name(boundary_mesh, "ORIGINAL_CELLS")  # TODO name is copied
    assert original_cells
    tmp = []
    for i, j in non_conformal_cells:
        tmp.append((original_cells.GetValue(i), original_cells.GetValue(j)))

    return Result(non_conformal_cells=tmp)


def __compute_volume(points: vtkPoints, face_stream: FaceStream) -> float:
    """
    Computes the volume of a polyhedron element, based on its face_stream
    :param face_stream:
    :return: The volume of the element.
    """
    polygons = vtkCellArray()
    for face_nodes in face_stream.face_nodes:
        polygon = vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(len(face_nodes))
        for i, point_id in enumerate(face_nodes):
            polygon.GetPointIds().SetId(i, point_id)
        polygons.InsertNextCell(polygon)
    polygon_poly_data = vtkPolyData()
    polygon_poly_data.SetPoints(points)  # TODO Is this properly working? Is the number of points OK?
    polygon_poly_data.SetPolys(polygons)
    f = vtk.vtkTriangleFilter()
    f.SetInputData(polygon_poly_data)
    f.Update()
    triangles = f.GetOutput()
    # Compute the barycenter TODO extract as a function?
    tmp_barycenter = numpy.empty((face_stream.num_support_points, 3), dtype=float)
    for i, point_id in enumerate(face_stream.support_point_ids):
        tmp_barycenter[i, :] = points.GetPoint(point_id)
    # x, y, z = tmp_barycenter[:, 0].mean(), tmp_barycenter[:, 1].mean(), tmp_barycenter[:, 2].mean()  # TODO dirty
    barycenter = tmp_barycenter[:, 0].mean(), tmp_barycenter[:, 1].mean(), tmp_barycenter[:, 2].mean()  # TODO dirty
    cell_volume = 0.
    for i in range(triangles.GetNumberOfCells()):
        triangle = triangles.GetCell(i)
        assert triangle.GetCellType() == VTK_TRIANGLE
        p = triangle.GetPoints()
        cell_volume += vtkTetra.ComputeVolume(p.GetPoint(0), p.GetPoint(1), p.GetPoint(2), barycenter)  # TODO invert the order to reach positive volume (more natural)
    return cell_volume


def __select_flip(mesh_points: vtkPoints,
                  colors: Dict[FrozenSet[int], int],
                  face_stream: FaceStream) -> FaceStream:
    """

    :param points:
    :param connected_components: List of connected faces.
    :param colors: Maps  the connected component index to its color.
    :param face_stream: the face.
    :return:
    """
    # Flipping either color 0 or 1.
    color_to_nodes: Dict[int, List[int]] = {0: [], 1: []}
    for connected_components_indices, color in colors.items():
        color_to_nodes[color] += connected_components_indices
    # TODO this should work even if there's only one unique color.
    fs = face_stream.flip_faces(color_to_nodes[0]), face_stream.flip_faces(color_to_nodes[1])
    volumes: Tuple[float, float] = __compute_volume(mesh_points, fs[0]), __compute_volume(mesh_points, fs[1])
    # print(volumes)
    return fs[numpy.argmin(volumes)]  # TODO we take the min because the barycenter is supposed to be "in the center" while the normals are outward.


def __reorient_element(mesh_points: vtkPoints, face_stream_ids: vtkIdList) -> vtkIdList:
    face_stream = FaceStream.build_from_vtk_id_list(face_stream_ids)
    face_graph = build_cell_graph(face_stream, add_compatibility=True)
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
    colors: Dict[FrozenSet[int], int] = networkx.algorithms.greedy_color(quotient_graph)
    assert len(colors) in (1, 2)
    # We now compute the face stream which generates outwards normal vectors.
    flipped_face_stream = __select_flip(mesh_points, colors, face_stream)
    return to_vtk_id_list(flipped_face_stream.dump())


def __reorient_mesh(mesh, cell_indices) -> vtkUnstructuredGrid:  # TODO do it only on boundary cells?
    output = vtkUnstructuredGrid()
    output.SetPoints(mesh.GetPoints())
    # output.AllocateExact(mesh.GetNumberOfCells())  # TODO what's the size here?
    logging.info("Reorienting the polyhedron cells to enforce normals directed outward.")
    for ic in tqdm(cell_indices, desc="Reorienting polyhedra"):
        cell = mesh.GetCell(ic)
        cell_type = cell.GetCellType()
        if cell_type != VTK_POLYHEDRON:
            output.InsertNextCell(cell_type, cell.GetPointIds())
        else:
            face_stream_ids = vtkIdList()
            mesh.GetFaceStream(ic, face_stream_ids)
            reoriented_faces = __reorient_element(mesh.GetPoints(), face_stream_ids)
            output.InsertNextCell(VTK_POLYHEDRON, reoriented_faces)
    assert output.GetNumberOfCells() == mesh.GetNumberOfCells()
    return output


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)

# mesh = vtk_utils.read_mesh('/docker-exchange/torn_wedges.vtu')
# # mesh = vtk_utils.read_mesh('/docker-exchange/torn_pebi.vtu')
# output = __reorient_mesh(mesh, range(mesh.GetNumberOfCells()))
# vtk_utils.write_mesh(output, "/docker-exchange/outward.vtk")
# result = __check(output, Options(angle_tolerance=10, point_tolerance=0.1, face_tolerance=1.))

# face_stream = FS
# points = PTS
