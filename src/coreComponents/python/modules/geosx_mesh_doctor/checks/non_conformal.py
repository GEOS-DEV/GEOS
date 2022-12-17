from dataclasses import dataclass
import logging
import math
from typing import List, Tuple
import numpy
import scipy.optimize

from vtkmodules.vtkCommonCore import (
    vtkIdList,
)
from vtkmodules.vtkCommonDataModel import (
    vtkBoundingBox,
    vtkPolyData,
    vtkStaticPointLocator,
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

from . import vtk_utils


@dataclass(frozen=True)
class Options:
    angle_tolerance: float
    point_tolerance: float
    face_tolerance: float


@dataclass(frozen=True)
class Result:
    non_conformal_cells: List[Tuple[int, int]]


# def _iter(id_list: vtkIdList) -> Iterator[int]:
#     """
#     Utility function transforming a vtkIdList into an iterable to be used for building built-ins python containers.
#     :param id_list: the vtkIdList.
#     :return: The iterator.
#     """
#     for i in range(id_list.GetNumberOfIds()):
#         yield id_list.GetId(i)


# def __max_dist(cell_center, cell):
#     distances = []
#     cc = numpy.array(cell_center)
#     cell_points = cell.GetPoints()
#     for ip in range(cell_points.GetNumberOfPoints()):
#         point = numpy.array(cell_points.GetPoint(ip))
#         distances.append(numpy.linalg.norm(cc - point))
#     return max(distances)

# dead code
# BOUNDARY_POINTS = "BOUNDARY_POINTS"
# BOUNDARY_FACES = "BOUNDARY_FACES"
# BOUNDARY_CELLS = "BOUNDARY_CELLS"
# f = vtkMarkBoundaryFilter()
# f.GenerateBoundaryFacesOn()
# f.SetBoundaryPointsName(BOUNDARY_POINTS)
# f.SetBoundaryFacesName(BOUNDARY_FACES)
# f.SetBoundaryCellsName(BOUNDARY_CELLS)
# f.SetInputData(mesh)

# loc = vtkStaticPointLocator()
# loc.SetDataSet(o)
# loc.BuildLocator()
#
# cc = vtkCellCenters()
# cc.VertexCellsOn()
# cc.CopyArraysOff()
# cc.SetInputData(o)
# cc.Update()
# cell_centers = cc.GetOutput()
#
# for ic in range(o.GetNumberOfCells()):
#     cell = o.GetCell(ic)
#     assert cell.GetCellDimension() == 2
#     center = cell_centers.GetPoint(ic) # compute cell.ComputeCentroid?
#     distance = __max_dist(center, cell)  # use cell.ComputeBoundingSphere?
#     result = vtkIdList()
#     loc.FindPointsWithinRadius(distance * 1.5, center, result)
#     result = set(_iter(result))
#     cell_point_ids = set(_iter(cell.GetPointIds()))
#     for close_point_id in filter(lambda r: r not in cell_point_ids, result):
#         print(close_point_id)


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
    n.ConsistencyOn()
    n.AutoOrientNormalsOn()
    n.ComputeCellNormalsOn()
    n.SetInputData(boundary_mesh)
    n.Update()
    # ncd = n.GetOutput().GetCellData()  # TODO build a get_field_by_name
    # for i in range(ncd.GetNumberOfArrays()):
    #     if ncd.GetArrayName(i) == "Normals":
    #         normals = ncd.GetArray(i)
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
    for i in range(num_cells):
        for j in range(i + 1, num_cells):  # TODO use some tree to find the appropriate neighbors?
            bbi = vtkBoundingBox(boundary_mesh.GetCell(i).GetBounds())  # TODO store this
            bbi.Inflate(options.face_tolerance)
            bbj = vtkBoundingBox(boundary_mesh.GetCell(j).GetBounds())
            bbj.Inflate(options.face_tolerance)
            if not bbi.Intersects(bbj):
                continue
            # TODO use vtkCellLocators (but the bounds are weak)
            ni, nj = normals.GetTuple3(i), normals.GetTuple3(j)
            if numpy.dot(ni, nj) > -cos_theta:  # opposite directions (can be facing or not)
                continue
            # if abs(numpy.dot(ni, nj)) > cos_theta:  # opposite directions (can be facing or not)
            #     continue
            cci, ccj = cell_centers[i], cell_centers[j]
            direction_i = numpy.dot(ccj - cci, ni)  # TODO use closest point to polygon instead (for ccj)
            direction_j = numpy.dot(cci - ccj, nj)
            # print(f"{i}, {j}, {scalar_product}, {direction_i}, {direction_j}")
            # if direction_i < 0 or direction_j < 0:
            #     continue
            logging.debug(f"Computing the distance for faces {i} and {j}.")
            result = scipy.optimize.minimize(lambda x: distance_between_cells(x, i, j, boundary_mesh),
                                             (0.5, 0.5, 0.5, 0.5),  # TODO do a better screening!
                                             method="Nelder-Mead",
                                             bounds=((0, 1), (0, 1), (0, 1), (0, 1))
                                             )
            # print(result)
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


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
