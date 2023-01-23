from dataclasses import dataclass
import itertools
import logging
import math
from typing import List, Tuple, Dict, FrozenSet, Iterator
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


def get_cell_field_by_name(mesh, field_name):  # TODO move to utilities
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


def build_boundary_mesh(mesh, consistency=True, flip_normals=True) -> Result:  # TODO Move to the BoundaryMesh class # TODO check the flip normals...
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
    n.SetConsistency(consistency)  # TODO review the options!
    n.SetFlipNormals(flip_normals)  # TODO review the options!
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


class BoundaryMesh:
    def __init__(self, mesh):
        self.mesh = mesh
        self.boundary_mesh, self.__normals, self.cell_centers = build_boundary_mesh(mesh)
        self.original_cells = get_cell_field_by_name(self.boundary_mesh, "ORIGINAL_CELLS")
        assert self.original_cells
        cells_to_reorient = map(self.original_cells.GetValue, range(self.original_cells.GetNumberOfValues()))
        reoriented_mesh = reorient_mesh(mesh, cells_to_reorient)
        # reoriented_mesh = reorient_mesh(mesh, range(mesh.GetNumberOfCells()))  # TODO only some cells...
        self.re_boundary_mesh, self.re_normals, self.re_cell_centers = build_boundary_mesh(reoriented_mesh, consistency=False)
    def GetNumberOfCells(self):
        return self.boundary_mesh.GetNumberOfCells()
    def bounds(self, i):
        return self.boundary_mesh.GetCell(i).GetBounds()
    def __underlying_original_cell_type(self, i):
        return self.mesh.GetCell(self.original_cells.GetValue(i)).GetCellType()
    def normals(self, i):
        if self.__underlying_original_cell_type(i) == VTK_POLYHEDRON:
            return self.re_normals.GetTuple3(i)
        else:
            return self.__normals.GetTuple3(i)
    def GetCell(self, i):
        if self.__underlying_original_cell_type(i) == VTK_POLYHEDRON:
            return self.re_boundary_mesh.GetCell(i)
        else:
            return self.boundary_mesh.GetCell(i)
    def GetPoint(self, i):
        return self.boundary_mesh.GetPoint(i)


def build_poly_data_for_extrusion(i, boundary_mesh):
    cell_i = boundary_mesh.GetCell(i)
    cp_i = cell_i.NewInstance()
    cp_i.DeepCopy(cell_i)
    points_ids_mapping_i = []
    for ii in range(cell_i.GetNumberOfPoints()):
        cp_i.GetPointIds().SetId(ii, ii)
        points_ids_mapping_i.append(cell_i.GetPointId(ii))
    polygons_i = vtkCellArray()
    polygons_i.InsertNextCell(cp_i)
    points_i = vtkPoints()
    points_i.SetNumberOfPoints(len(points_ids_mapping_i))
    for ii, v in enumerate(points_ids_mapping_i):
        points_i.SetPoint(ii, boundary_mesh.GetPoint(v))
    polygon_poly_data_i = vtkPolyData()
    polygon_poly_data_i.SetPoints(points_i)
    polygon_poly_data_i.SetPolys(polygons_i)
    return polygon_poly_data_i


def test(i, j, normal_i, normal_j, boundary_mesh, face_tolerance):
    polygon_poly_data_i = build_poly_data_for_extrusion(i, boundary_mesh)
    polygon_poly_data_j = build_poly_data_for_extrusion(j, boundary_mesh)

    # polygons_i = vtkCellArray()
    # polygons_i.InsertNextCell(boundary_mesh.GetCell(i))
    # polygons_j = vtkCellArray()
    # polygons_j.InsertNextCell(boundary_mesh.GetCell(j))
    # polygon_poly_data_i = vtkPolyData()
    # polygon_poly_data_i.SetPoints(boundary_mesh.GetPoints())  # TODO Is this properly working? Is the number of points OK?
    # polygon_poly_data_i.SetPolys(polygons_i)
    # polygon_poly_data_j = vtkPolyData()
    # polygon_poly_data_j.SetPoints(boundary_mesh.GetPoints())  # TODO Is this properly working? Is the number of points OK?
    # polygon_poly_data_j.SetPolys(polygons_j)
    extruder_i = vtk.vtkLinearExtrusionFilter()
    # extruder_i.SetExtrusionTypeToNormalExtrusion()
    extruder_i.SetExtrusionTypeToVectorExtrusion()
    extruder_i.SetVector(normal_i)
    extruder_i.SetScaleFactor(face_tolerance)
    # extruder_i.SetScaleFactor(100)
    extruder_i.SetInputData(polygon_poly_data_i)
    extruder_i.Update()
    # e_i = extruder_i.GetOutput()
    # print(boundary_mesh.GetCell(i).GetNumberOfPoints())
    # print(boundary_mesh.GetCell(i).GetCellType())
    # print(e_i.GetNumberOfCells())
    extruder_j = vtk.vtkLinearExtrusionFilter()
    # extruder_j.SetExtrusionTypeToNormalExtrusion()
    extruder_j.SetExtrusionTypeToVectorExtrusion()
    extruder_i.SetVector(normal_j)
    extruder_j.SetScaleFactor(face_tolerance)
    # extruder_j.SetScaleFactor(100)
    extruder_j.SetInputData(polygon_poly_data_j)
    extruder_j.Update()

    # b = vtk.vtkLoopBooleanPolyDataFilter()
    coll = vtk.vtkCollisionDetectionFilter()
    coll.SetCollisionModeToFirstContact()
    coll.SetInputData(0, extruder_i.GetOutput())
    coll.SetInputData(1, extruder_j.GetOutput())
    # m_i = vtk.vtkMatrix4x4()
    # m_j = vtk.vtkMatrix4x4()
    # coll.SetMatrix(0, m_i)
    # coll.SetMatrix(1, m_j)
    m_i = vtk.vtkTransform()
    m_j = vtk.vtkTransform()
    coll.SetTransform(0, m_i)
    coll.SetTransform(1, m_j)
    coll.Update()

    # coll2 = vtk.vtkCollisionDetectionFilter()
    # coll2.SetCollisionModeToAllContacts()
    # coll2.SetInputData(0, extruder_i.GetOutput())
    # coll2.SetInputData(1, extruder_j.GetOutput())
    # # m_i = vtk.vtkMatrix4x4()
    # # m_j = vtk.vtkMatrix4x4()
    # # coll.SetMatrix(0, m_i)
    # # coll.SetMatrix(1, m_j)
    # m_i = vtk.vtkTransform()
    # m_j = vtk.vtkTransform()
    # coll2.SetTransform(0, m_i)
    # coll2.SetTransform(1, m_j)
    # coll2.Update()

    print(coll.GetNumberOfContacts())
    # TODO add a vtkPointLocator and its FindClosestInsertedPoint methods to check that all points have their peer...
    return coll.GetNumberOfContacts()

    # # e_j = extruder_j.GetOutput()
    # intersection = vtk.vtkIntersectionPolyDataFilter()
    # intersection.SetInputConnection(0, extruder_i.GetOutputPort())
    # intersection.SetInputConnection(1, extruder_j.GetOutputPort())
    # intersection.ComputeIntersectionPointArrayOn()
    # intersection.Update()
    # oo = intersection.GetOutput()
    # print(oo.GetNumberOfCells(), i, j, normal_i, normal_j)


def __check(mesh, options: Options) -> Result:
    # boundary_mesh, normals, cell_centers = build_boundary_mesh(mesh)
    bm = BoundaryMesh(mesh)
    cos_theta = abs(math.cos(numpy.deg2rad(options.angle_tolerance)))

    non_conformal_cells = []

    # num_cells = boundary_mesh.GetNumberOfCells()
    num_cells = bm.GetNumberOfCells()

    # Precomputing the bounding boxes of all the boundary cells.
    bounding_boxes = []
    for i in range(num_cells):
        # bb = vtkBoundingBox(boundary_mesh.GetCell(i).GetBounds())
        bb = vtkBoundingBox(bm.bounds(i))
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
        # normal_i, normal_j = normals.GetTuple3(i), normals.GetTuple3(j)
        normal_i, normal_j = bm.normals(i), bm.normals(j)
        if numpy.dot(normal_i, normal_j) > -cos_theta:  # opposite directions only (can be facing or not)
            continue
        # cci, ccj = cell_centers[i], cell_centers[j]
        # direction_i = numpy.dot(ccj - cci, normal_i)  # TODO use closest point to polygon instead (for ccj)
        # direction_j = numpy.dot(cci - ccj, normal_j)  # TODO use
        # # print(f"{i}, {j}, {scalar_product}, {direction_i}, {direction_j}")
        # if direction_i < 0 or direction_j < 0:  # checking direction now.
        #     continue
        # if test(i, j, normal_i, normal_j, boundary_mesh, options.face_tolerance):
        if test(i, j, normal_i, normal_j, bm, options.face_tolerance):
            non_conformal_cells.append((i, j))
        # logging.debug(f"Computing the distance for faces {i} and {j}.")
        # result = scipy.optimize.minimize(lambda x: distance_between_cells(x, i, j, boundary_mesh),  # TODO that somehow sucks
        #                                  (0.5, 0.5, 0.5, 0.5),  # TODO do a better screening!
        #                                  method="Nelder-Mead",
        #                                  bounds=((0, 1), (0, 1), (0, 1), (0, 1))
        #                                  )
        # # TODO deal with error code.
        # logging.debug(f"Result between {i} and {j} is {result.fun}")
        # if abs(result.fun) < options.face_tolerance:
        #     non_conformal_cells.append((i, j))

    # original_cells = get_cell_field_by_name(bm.boundary_mesh, "ORIGINAL_CELLS")  # TODO name is copied
    # assert original_cells
    tmp = []
    for i, j in non_conformal_cells:
        tmp.append((bm.original_cells.GetValue(i), bm.original_cells.GetValue(j)))

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


def reorient_mesh(mesh, cell_indices: Iterator[int]) -> vtkUnstructuredGrid:
    num_cells = mesh.GetNumberOfCells()
    # Building an indicator/predicate from the list
    needs_to_be_reoriented = numpy.zeros(num_cells, dtype=bool)
    for ic in cell_indices:
        needs_to_be_reoriented[ic] = True

    output_mesh = vtkUnstructuredGrid()
    output_mesh.SetPoints(mesh.GetPoints())
    # output.AllocateExact(mesh.GetNumberOfCells())  # TODO what's the size here?
    logging.info("Reorienting the polyhedron cells to enforce normals directed outward.")
    # for ic in tqdm(cell_indices, desc="Reorienting polyhedra"):
    for ic in tqdm(range(num_cells), desc="Reorienting polyhedra"):
        cell = mesh.GetCell(ic)
        cell_type = cell.GetCellType()
        if cell_type == VTK_POLYHEDRON:
            face_stream_ids = vtkIdList()
            mesh.GetFaceStream(ic, face_stream_ids)
            new_face_stream_ids = __reorient_element(mesh.GetPoints(), face_stream_ids) if needs_to_be_reoriented[ic] else face_stream_ids
            output_mesh.InsertNextCell(VTK_POLYHEDRON, new_face_stream_ids)
        else:
            output_mesh.InsertNextCell(cell_type, cell.GetPointIds())
    assert output_mesh.GetNumberOfCells() == mesh.GetNumberOfCells()
    return output_mesh


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
    # tmp_mesh = reorient_mesh(mesh, range(mesh.GetNumberOfCells()))
    # return __check(tmp_mesh, options)

mesh = vtk_utils.read_mesh('/docker-exchange/torn_wedges.vtu')
# mesh = vtk_utils.read_mesh('/docker-exchange/torn_pebi.vtu')
# mesh = vtk_utils.read_mesh('/docker-exchange/torn_pebi_anticline.vtu')
# output = reorient_mesh(mesh, range(mesh.GetNumberOfCells()))
# vtk_utils.write_mesh(output, "/docker-exchange/outward.vtk")
# result = __check(output, Options(angle_tolerance=10, point_tolerance=0.1, face_tolerance=1.))
result = __check(mesh, Options(angle_tolerance=10, point_tolerance=0.1, face_tolerance=1.))

t = []
for pair in result.non_conformal_cells:
    t += pair
t = set(t)

# face_stream = FS
# points = PTS
