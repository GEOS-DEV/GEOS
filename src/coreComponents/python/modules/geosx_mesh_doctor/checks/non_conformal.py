from dataclasses import dataclass
import logging
import math
from typing import List, Tuple, Dict, FrozenSet, Iterator
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
    vtkBoundingBox,
    vtkCellArray,
    vtkPointSet,
    vtkPolyData,
    vtkPolygon,
    vtkPolyhedron,
    vtkStaticCellLocator,
    vtkStaticPointLocator,
    vtkUnstructuredGrid,
    vtkTetra,
)
from vtkmodules.vtkCommonTransforms import (
    vtkTransform,
)
from vtkmodules.vtkFiltersCore import (
    vtkPolyDataNormals,
    vtkTriangleFilter,
)
from vtkmodules.vtkFiltersGeometry import (
    vtkDataSetSurfaceFilter,
)
from vtkmodules.vtkFiltersModeling import (
    vtkCollisionDetectionFilter,
    vtkLinearExtrusionFilter,
)
from vtk import reference as vtk_reference

from . import vtk_utils

from .vtk_polyhedron import (
    FaceStream,
    build_cell_graph,
    to_vtk_id_list,
    _iter,
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


class BoundaryMesh:
    def __init__(self, mesh):
        # Building the boundary meshes
        boundary_mesh, __normals, self.__original_cells = BoundaryMesh.__build_boundary_mesh(mesh)
        cells_to_reorient = filter(lambda c: mesh.GetCell(c).GetCellType() == VTK_POLYHEDRON,
                                   map(self.__original_cells.GetValue,
                                       range(self.__original_cells.GetNumberOfValues())))
        reoriented_mesh = reorient_mesh(mesh, cells_to_reorient)
        self.re_boundary_mesh, re_normals, _ = BoundaryMesh.__build_boundary_mesh(reoriented_mesh, consistency=False)
        num_cells = boundary_mesh.GetNumberOfCells()
        # Precomputing the underlying cell type
        self.__is_underlying_cell_type_a_polyhedron = numpy.zeros(num_cells, dtype=bool)
        for ic in range(num_cells):
            self.__is_underlying_cell_type_a_polyhedron[ic] = mesh.GetCell(self.__original_cells.GetValue(ic)).GetCellType() == VTK_POLYHEDRON
        # Precomputing the normals
        self.__normals = numpy.empty((num_cells, 3), dtype=numpy.double, order='C')  # Do not modify the storage layout
        for ic in range(num_cells):
            if self.__is_underlying_cell_type_a_polyhedron[ic]:
                self.__normals[ic, :] = re_normals.GetTuple3(ic)
            else:
                self.__normals[ic, :] = __normals.GetTuple3(ic)
    @staticmethod
    def __build_boundary_mesh(mesh, consistency=True):
        original_cells_key = "ORIGINAL_CELLS"

        f = vtkDataSetSurfaceFilter()
        f.PassThroughCellIdsOn()
        f.PassThroughPointIdsOff()
        f.FastModeOff()

        f.SetOriginalCellIdsName(original_cells_key)

        boundary_mesh = vtkPolyData()
        f.UnstructuredGridExecute(mesh, boundary_mesh)

        n = vtkPolyDataNormals()
        n.SetConsistency(consistency)
        n.SetAutoOrientNormals(consistency)
        n.FlipNormalsOff()
        n.ComputeCellNormalsOn()
        n.SetInputData(boundary_mesh)
        n.Update()
        normals = get_cell_field_by_name(n.GetOutput(), "Normals")
        assert normals
        assert normals.GetNumberOfComponents() == 3
        assert normals.GetNumberOfTuples() == boundary_mesh.GetNumberOfCells()
        original_cells = get_cell_field_by_name(boundary_mesh, original_cells_key)
        assert original_cells
        return boundary_mesh, normals, original_cells
    def GetNumberOfCells(self):
        return self.re_boundary_mesh.GetNumberOfCells()
    def GetNumberOfPoints(self):
        return self.re_boundary_mesh.GetNumberOfPoints()
    def bounds(self, i):
        return self.re_boundary_mesh.GetCell(i).GetBounds()
    def normals(self, i):
        return self.__normals[i]
    def GetCell(self, i):
        return self.re_boundary_mesh.GetCell(i)
    def GetPoint(self, i):
        return self.re_boundary_mesh.GetPoint(i)
    @property
    def original_cells(self):
        return self.__original_cells


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
    return polygon_poly_data_i, cp_i


def test(i, j, normal_i, normal_j, boundary_mesh, face_tolerance, point_tolerance):
    polygon_poly_data_i, cp_i = build_poly_data_for_extrusion(i, boundary_mesh)
    polygon_poly_data_j, cp_j = build_poly_data_for_extrusion(j, boundary_mesh)

    # TODO deduplicate
    extruder_i = vtkLinearExtrusionFilter()
    extruder_i.SetExtrusionTypeToVectorExtrusion()
    extruder_i.SetVector(normal_i)
    extruder_i.SetScaleFactor(face_tolerance)
    extruder_i.SetInputData(polygon_poly_data_i)
    extruder_i.Update()

    extruder_j = vtkLinearExtrusionFilter()
    extruder_j.SetExtrusionTypeToVectorExtrusion()
    extruder_i.SetVector(normal_j)
    extruder_j.SetScaleFactor(face_tolerance)
    extruder_j.SetInputData(polygon_poly_data_j)
    extruder_j.Update()

    collision = vtkCollisionDetectionFilter()
    collision.SetCollisionModeToFirstContact()
    collision.SetInputData(0, extruder_i.GetOutput())
    collision.SetInputData(1, extruder_j.GetOutput())
    # m_i = vtk.vtkMatrix4x4()
    # m_j = vtk.vtkMatrix4x4()
    # collision.SetMatrix(0, m_i)
    # collision.SetMatrix(1, m_j)
    m_i = vtkTransform()
    m_j = vtkTransform()
    collision.SetTransform(0, m_i)
    collision.SetTransform(1, m_j)
    collision.Update()

    if collision.GetNumberOfContacts() == 0:
        return False

    # In this last step, we check that the nodes are (or not) matching each other.
    if cp_i.GetNumberOfPoints() != cp_j.GetNumberOfPoints():
        return True

    point_locator = vtkStaticPointLocator()
    points = vtkPointSet()
    points.SetPoints(cp_i.GetPoints())
    point_locator.SetDataSet(points)
    point_locator.BuildLocator()
    found_points = set()
    for ip in range(cp_j.GetNumberOfPoints()):
        p = cp_j.GetPoints().GetPoint(ip)
        squared_dist = vtk_reference(0.)
        found_point = point_locator.FindClosestPointWithinRadius(point_tolerance, p, squared_dist)
        found_points.add(found_point)
    return found_points != set(range(cp_i.GetNumberOfPoints()))


def __check(mesh, options: Options) -> Result:
    boundary_mesh = BoundaryMesh(mesh)
    cos_theta = abs(math.cos(numpy.deg2rad(options.angle_tolerance)))
    num_cells = boundary_mesh.GetNumberOfCells()

    # Computing the exact number of cells per node
    num_cells_per_node = numpy.zeros(boundary_mesh.GetNumberOfPoints(), dtype=int)
    for ic in range(boundary_mesh.GetNumberOfCells()):
        c = boundary_mesh.GetCell(ic)
        point_ids = c.GetPointIds()
        for point_id in _iter(point_ids):
            num_cells_per_node[point_id] += 1

    cell_locator = vtkStaticCellLocator()
    cell_locator.Initialize()
    cell_locator.SetNumberOfCellsPerNode(num_cells_per_node.max())
    cell_locator.SetDataSet(boundary_mesh.re_boundary_mesh)
    cell_locator.BuildLocator()

    # Precomputing the bounding boxes.
    bounding_boxes = numpy.empty((boundary_mesh.GetNumberOfCells(), 6), dtype=numpy.double, order="C")
    for i in range(boundary_mesh.GetNumberOfCells()):
        bb = vtkBoundingBox(boundary_mesh.bounds(i))
        bb.Inflate(2 * options.face_tolerance)
        assert bounding_boxes[i, :].data.contiguous  # Do not modify the storage layout since vtk deals with raw memory here.
        bb.GetBounds(bounding_boxes[i, :])

    non_conformal_cells = []
    close_cells = vtkIdList()
    # Looping on all the pairs of boundary cells. We'll hopefully discard most of the pairs.
    for i in tqdm(range(num_cells), desc="Non conformal elements"):
        cell_locator.FindCellsWithinBounds(bounding_boxes[i], close_cells)
        for j in _iter(close_cells):
            if j < i:
                continue
            # Discarding pairs that are not facing each others (with a threshold).
            normal_i, normal_j = boundary_mesh.normals(i), boundary_mesh.normals(j)
            if numpy.dot(normal_i, normal_j) > -cos_theta:  # opposite directions only (can be facing or not)
                continue
            # cci, ccj = cell_centers[i], cell_centers[j]
            # direction_i = numpy.dot(ccj - cci, normal_i)  # TODO use closest point to polygon instead (for ccj)
            # direction_j = numpy.dot(cci - ccj, normal_j)  # TODO use
            # # print(f"{i}, {j}, {scalar_product}, {direction_i}, {direction_j}")
            # if direction_i < 0 or direction_j < 0:  # checking direction now.
            #     continue
            if test(i, j, normal_i, normal_j, boundary_mesh, options.face_tolerance, options.point_tolerance):
                non_conformal_cells.append((i, j))

    tmp = []
    for i, j in non_conformal_cells:
        tmp.append((boundary_mesh.original_cells.GetValue(i), boundary_mesh.original_cells.GetValue(j)))

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
    f = vtkTriangleFilter()
    f.SetInputData(polygon_poly_data)
    f.Update()
    triangles = f.GetOutput()
    # Compute the barycenter TODO extract as a function?
    tmp_barycenter = numpy.empty((face_stream.num_support_points, 3), dtype=float)
    for i, point_id in enumerate(face_stream.support_point_ids):
        tmp_barycenter[i, :] = points.GetPoint(point_id)
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

    output_mesh = vtkUnstructuredGrid()  # TODO use CopyStructure here.
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
                    progress_bar.update(1)
                else:
                    new_face_stream_ids = face_stream_ids
                output_mesh.InsertNextCell(VTK_POLYHEDRON, new_face_stream_ids)
            else:
                output_mesh.InsertNextCell(cell_type, cell.GetPointIds())
                if needs_to_be_reoriented[ic]:
                    progress_bar.update(1)
    assert output_mesh.GetNumberOfCells() == mesh.GetNumberOfCells()
    return output_mesh


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)


# mesh = vtk_utils.read_mesh('/docker-exchange/torn_wedges.vtu')
# # mesh = vtk_utils.read_mesh('/docker-exchange/torn_pebi.vtu')
# # mesh = vtk_utils.read_mesh('/docker-exchange/torn_pebi_anticline.vtu')
# # output = reorient_mesh(mesh, range(mesh.GetNumberOfCells()))
# # vtk_utils.write_mesh(output, "/docker-exchange/outward.vtk")
# # result = __check(output, Options(angle_tolerance=10, point_tolerance=0.1, face_tolerance=1.))
# result = __check(mesh, Options(angle_tolerance=10, point_tolerance=0.1, face_tolerance=1.))
#
# t = []
# for pair in result.non_conformal_cells:
#     t += pair
# t = set(t)

# face_stream = FS
# points = PTS
