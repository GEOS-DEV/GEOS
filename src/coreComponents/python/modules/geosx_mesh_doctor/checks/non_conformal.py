from dataclasses import dataclass
import math
from typing import List, Tuple, Any
import numpy

from tqdm import tqdm

from vtkmodules.vtkCommonCore import (
    vtkIdList,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_POLYHEDRON,
    vtkBoundingBox,
    vtkCell,
    vtkCellArray,
    vtkPointSet,
    vtkPolyData,
    vtkStaticCellLocator,
    vtkStaticPointLocator,
    vtkUnstructuredGrid,
)
from vtkmodules.vtkCommonTransforms import (
    vtkTransform,
)
from vtkmodules.vtkFiltersCore import (
    vtkPolyDataNormals,
)
from vtkmodules.vtkFiltersGeometry import (
    vtkDataSetSurfaceFilter,
)
from vtkmodules.vtkFiltersModeling import (
    vtkCollisionDetectionFilter,
    vtkLinearExtrusionFilter,
)
from vtk import reference as vtk_reference

from .reorient_mesh import reorient_mesh

from . import vtk_utils

from .vtk_polyhedron import (
    vtk_iter,
)

from . import triangle_distance


@dataclass(frozen=True)
class Options:
    angle_tolerance: float
    point_tolerance: float
    face_tolerance: float


@dataclass(frozen=True)
class Result:
    non_conformal_cells: List[Tuple[int, int]]


class BoundaryMesh:
    """
    A BoundaryMesh is the envelope of the 3d mesh on which we want to perform the simulations.
    It is computed by vtk. But we want to be sure that the normals of the envelope are directed outwards.
    The `vtkDataSetSurfaceFilter` does not have the same behavior for standard vtk cells (like tets or hexs),
    and for polyhedron meshes, for which the result is a bit brittle.
    Therefore, we reorient the polyhedron cells ourselves, so we're sure that they point outwards.
    And then we compute the boundary meshes for both meshes, given that the computing options are not identical.
    """
    def __init__(self, mesh: vtkUnstructuredGrid):
        """
        Builds a boundary mesh.
        :param mesh: The 3d mesh.
        """
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
    def __build_boundary_mesh(mesh: vtkUnstructuredGrid, consistency=True) -> Tuple[vtkUnstructuredGrid, Any, Any]:
        """
        From a 3d mesh, build the envelope meshes.
        :param mesh: The input 3d mesh.
        :param consistency: The vtk option passed to the `vtkDataSetSurfaceFilter`.
        :return: A tuple containing the boundary mesh, the normal vectors array,
                 an array that maps the id of the boundary element to the id of the 3d cell it touches.
        """
        f = vtkDataSetSurfaceFilter()
        f.PassThroughCellIdsOn()
        f.PassThroughPointIdsOff()
        f.FastModeOff()

        # Note that we do not need the original points, but we could keep them as well if needed
        original_cells_key = "ORIGINAL_CELLS"
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
        normals = n.GetOutput().GetCellData().GetArray("Normals")
        assert normals
        assert normals.GetNumberOfComponents() == 3
        assert normals.GetNumberOfTuples() == boundary_mesh.GetNumberOfCells()
        original_cells = boundary_mesh.GetCellData().GetArray(original_cells_key)
        assert original_cells
        return boundary_mesh, normals, original_cells

    def GetNumberOfCells(self) -> int:
        """
        The number of cells.
        :return: An integer.
        """
        return self.re_boundary_mesh.GetNumberOfCells()

    def GetNumberOfPoints(self) -> int:
        """
        The number of points.
        :return: An integer.
        """
        return self.re_boundary_mesh.GetNumberOfPoints()

    def bounds(self, i) -> vtkBoundingBox:
        """
        The boundrary box of cell `i`.
        :param i: The boundary cell index.
        :return: The vtk bounding box.
        """
        return self.re_boundary_mesh.GetCell(i).GetBounds()

    def normals(self, i) -> numpy.array:
        """
        The normal of cell `i`. This normal will be directed outwards
        :param i: The boundary cell index.
        :return: The normal as a length-3 numpy array.
        """
        return self.__normals[i]

    def GetCell(self, i) -> vtkCell:
        """
        Cell i of the boundary mesh. This cell will have its normal directed outwards.
        :param i: The boundary cell index.
        :return: The cell instance.
        :warning: This member function relies on the vtkUnstructuredGrid.GetCell member function which is not thread safe.
        """
        return self.re_boundary_mesh.GetCell(i)

    def GetPoint(self, i) -> Tuple[float, float, float]:
        """
        Point i of the boundary mesh.
        :param i: The boundary point index.
        :return: A length-3 tuple containing the coordinates of the point.
        :warning: This member function relies on the vtkUnstructuredGrid.GetPoint member function which is not thread safe.
        """
        return self.re_boundary_mesh.GetPoint(i)

    @property
    def original_cells(self):
        """
        Returns the 2d boundary cell to the 3d cell index of the original mesh.
        :return: A 1d array.
        """
        return self.__original_cells


def build_poly_data_for_extrusion(i: int, boundary_mesh: BoundaryMesh) -> vtkPolyData:
    """
    Creates a vtkPolyData containing the unique cell `i` of the boundary mesh.
    This operation is needed to use the vtk extrusion filter.
    :param i: The boundary cell index that will eventually be extruded.
    :param boundary_mesh:
    :return: The created vtkPolyData.
    """
    cell = boundary_mesh.GetCell(i)
    copied_cell = cell.NewInstance()
    copied_cell.DeepCopy(cell)
    points_ids_mapping = []
    for i in range(copied_cell.GetNumberOfPoints()):
        copied_cell.GetPointIds().SetId(i, i)
        points_ids_mapping.append(cell.GetPointId(i))
    polygons = vtkCellArray()
    polygons.InsertNextCell(copied_cell)
    points = vtkPoints()
    points.SetNumberOfPoints(len(points_ids_mapping))
    for i, v in enumerate(points_ids_mapping):
        points.SetPoint(i, boundary_mesh.GetPoint(v))
    polygon_poly_data = vtkPolyData()
    polygon_poly_data.SetPoints(points)
    polygon_poly_data.SetPolys(polygons)
    return polygon_poly_data


def are_points_conformal(point_tolerance: float, cell_i: vtkCell, cell_j: vtkCell) -> bool:
    """
    Checks if points of cell `i` matches, one by one, the points of cell `j`.
    :param point_tolerance: The point tolerance to consider that two points match.
    :param cell_i: The first cell.
    :param cell_j: The second cell.
    :return: A boolean.
    """
    # In this last step, we check that the nodes are (or not) matching each other.
    if cell_i.GetNumberOfPoints() != cell_j.GetNumberOfPoints():
        return True

    point_locator = vtkStaticPointLocator()
    points = vtkPointSet()
    points.SetPoints(cell_i.GetPoints())
    point_locator.SetDataSet(points)
    point_locator.BuildLocator()
    found_points = set()
    for ip in range(cell_j.GetNumberOfPoints()):
        p = cell_j.GetPoints().GetPoint(ip)
        squared_dist = vtk_reference(0.)  # unused
        found_point = point_locator.FindClosestPointWithinRadius(point_tolerance, p, squared_dist)
        found_points.add(found_point)
    return found_points == set(range(cell_i.GetNumberOfPoints()))


class Extruder:
    """
    Computes and stores all the extrusions of the boundary faces.
    The main reason for this class is to be lazy and cache the extrusions.
    """
    def __init__(self, boundary_mesh: BoundaryMesh, face_tolerance: float):
        self.__extrusions: List[vtkPolyData] = [None, ] * boundary_mesh.GetNumberOfCells()
        self.__boundary_mesh = boundary_mesh
        self.__face_tolerance = face_tolerance

    def __extrude(self, polygon_poly_data, normal) -> vtkPolyData:
        """
        Extrude the polygon data to create a volume that will be used for intersection.
        :param polygon_poly_data: The data to extrude
        :param normal: The (uniform) direction of the extrusion.
        :return: The extrusion.
        """
        extruder = vtkLinearExtrusionFilter()
        extruder.SetExtrusionTypeToVectorExtrusion()
        extruder.SetVector(normal)
        extruder.SetScaleFactor(self.__face_tolerance / 2.)
        extruder.SetInputData(polygon_poly_data)
        extruder.Update()
        return extruder.GetOutput()

    def __getitem__(self, i) -> vtkPolyData:
        """
        Returns the vtk extrusion for boundary element i.
        :param i: The cell index.
        :return: The vtk instance.
        """
        extrusion = self.__extrusions[i]
        if extrusion:
            return extrusion
        extrusion = self.__extrude(build_poly_data_for_extrusion(i, self.__boundary_mesh),
                                   self.__boundary_mesh.normals(i))
        self.__extrusions[i] = extrusion
        return extrusion


def are_faces_conformal_using_extrusions(extrusions: Extruder,
                                         i: int, j: int,
                                         boundary_mesh: vtkUnstructuredGrid,
                                         point_tolerance: float) -> bool:
    """
    Tests if two boundary faces are conformal, checking for intersection between their normal extruded volumes.
    :param extrusions: The extrusions cache.
    :param i: The cell index of the first cell.
    :param j: The cell index of the second cell.
    :param boundary_mesh: The boundary mesh.
    :param point_tolerance: The point tolerance to consider that two points match.
    :return: A boolean.
    """
    collision = vtkCollisionDetectionFilter()
    collision.SetCollisionModeToFirstContact()
    collision.SetInputData(0, extrusions[i])
    collision.SetInputData(1, extrusions[j])
    m_i = vtkTransform()
    m_j = vtkTransform()
    collision.SetTransform(0, m_i)
    collision.SetTransform(1, m_j)
    collision.Update()

    if collision.GetNumberOfContacts() == 0:
        return True

    # Duplicating data not to risk anything w.r.t. thread safety of the GetCell function.
    cell_i = boundary_mesh.GetCell(i)
    copied_cell_i = cell_i.NewInstance()
    copied_cell_i.DeepCopy(cell_i)

    return are_points_conformal(point_tolerance, copied_cell_i, boundary_mesh.GetCell(j))


def are_faces_conformal_using_distances(i: int, j: int,
                                        boundary_mesh: vtkUnstructuredGrid,
                                        face_tolerance: float, point_tolerance: float) -> bool:
    """
    Tests if two boundary faces are conformal, checking the minimal distance between triangulated surfaces.
    :param i: The cell index of the first cell.
    :param j: The cell index of the second cell.
    :param boundary_mesh: The boundary mesh.
    :param face_tolerance: The tolerance under which we should consider the two faces "touching" each other.
    :param point_tolerance: The point tolerance to consider that two points match.
    :return: A boolean.
    """
    cp_i = boundary_mesh.GetCell(i).NewInstance()
    cp_i.DeepCopy(boundary_mesh.GetCell(i))
    cp_j = boundary_mesh.GetCell(j).NewInstance()
    cp_j.DeepCopy(boundary_mesh.GetCell(j))

    def triangulate(cell):
        assert cell.GetCellDimension() == 2
        __points_ids = vtkIdList()
        __points = vtkPoints()
        cell.Triangulate(0, __points_ids, __points)
        __points_ids = tuple(vtk_iter(__points_ids))
        assert len(__points_ids) % 3 == 0
        assert __points.GetNumberOfPoints() % 3 == 0
        return __points_ids, __points

    points_ids_i, points_i = triangulate(cp_i)
    points_ids_j, points_j = triangulate(cp_j)

    def build_numpy_triangles(points_ids):
        __triangles = []
        for __i in range(0, len(points_ids), 3):
            __t = []
            for __pi in points_ids[__i: __i + 3]:
                __t.append(boundary_mesh.GetPoint(__pi))
            __triangles.append(numpy.array(__t, dtype=float))
        return __triangles

    triangles_i = build_numpy_triangles(points_ids_i)
    triangles_j = build_numpy_triangles(points_ids_j)

    min_dist = numpy.inf
    for ti, tj in [(ti, tj) for ti in triangles_i for tj in triangles_j]:
        # Note that here, we compute the exact distance to compare with the threshold.
        # We could improve by exiting the iterative distance computation as soon as
        # we're sure we're smaller than the threshold. No need of the exact solution.
        dist, _, _ = triangle_distance.distance_between_two_triangles(ti, tj)
        if dist < min_dist:
            min_dist = dist
        if min_dist < face_tolerance:
            break
    if min_dist > face_tolerance:
        return True

    return are_points_conformal(point_tolerance, cp_i, cp_j)


def __check(mesh: vtkUnstructuredGrid, options: Options) -> Result:
    """
    Checks if the mesh is "conformal" (i.e. if some of its boundary faces may not be too close to each other without matching nodes).
    :param mesh: The vtk mesh
    :param options: The check options.
    :return: The Result instance.
    """
    boundary_mesh = BoundaryMesh(mesh)
    cos_theta = abs(math.cos(numpy.deg2rad(options.angle_tolerance)))
    num_cells = boundary_mesh.GetNumberOfCells()

    # Computing the exact number of cells per node
    num_cells_per_node = numpy.zeros(boundary_mesh.GetNumberOfPoints(), dtype=int)
    for ic in range(boundary_mesh.GetNumberOfCells()):
        c = boundary_mesh.GetCell(ic)
        point_ids = c.GetPointIds()
        for point_id in vtk_iter(point_ids):
            num_cells_per_node[point_id] += 1

    cell_locator = vtkStaticCellLocator()
    cell_locator.Initialize()
    cell_locator.SetNumberOfCellsPerNode(num_cells_per_node.max())
    cell_locator.SetDataSet(boundary_mesh.re_boundary_mesh)
    cell_locator.BuildLocator()

    # Precomputing the bounding boxes.
    # The options are important to directly interact with memory in C++.
    bounding_boxes = numpy.empty((boundary_mesh.GetNumberOfCells(), 6), dtype=numpy.double, order="C")
    for i in range(boundary_mesh.GetNumberOfCells()):
        bb = vtkBoundingBox(boundary_mesh.bounds(i))
        bb.Inflate(2 * options.face_tolerance)
        assert bounding_boxes[i, :].data.contiguous  # Do not modify the storage layout since vtk deals with raw memory here.
        bb.GetBounds(bounding_boxes[i, :])

    non_conformal_cells = []
    extrusions = Extruder(boundary_mesh, options.face_tolerance)
    close_cells = vtkIdList()
    # Looping on all the pairs of boundary cells. We'll hopefully discard most of the pairs.
    for i in tqdm(range(num_cells), desc="Non conformal elements"):
        cell_locator.FindCellsWithinBounds(bounding_boxes[i], close_cells)
        for j in vtk_iter(close_cells):
            if j < i:
                continue
            # Discarding pairs that are not facing each others (with a threshold).
            normal_i, normal_j = boundary_mesh.normals(i), boundary_mesh.normals(j)
            if numpy.dot(normal_i, normal_j) > -cos_theta:  # opposite directions only (can be facing or not)
                continue
            # At this point, back-to-back and face-to-face pairs of elements are considered.
            if not are_faces_conformal_using_extrusions(extrusions, i, j, boundary_mesh, options.point_tolerance):
                non_conformal_cells.append((i, j))
    # Extracting the original 3d element index (and not the index of the boundary mesh).
    tmp = []
    for i, j in non_conformal_cells:
        tmp.append((boundary_mesh.original_cells.GetValue(i), boundary_mesh.original_cells.GetValue(j)))

    return Result(non_conformal_cells=tmp)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
