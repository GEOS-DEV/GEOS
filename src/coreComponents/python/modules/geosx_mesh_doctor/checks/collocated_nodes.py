from collections import defaultdict
from dataclasses import dataclass
import logging
from typing import (
    Collection,
    Iterable,
)
import numpy

from vtkmodules.vtkCommonCore import (
    reference,
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    vtkIncrementalOctreePointLocator, )

from . import vtk_utils


@dataclass(frozen=True)
class Options:
    tolerance: float


@dataclass(frozen=True)
class Result:
    nodes_buckets: Iterable[Iterable[int]]  # Each bucket contains the duplicated node indices.
    wrong_support_elements: Collection[int]  # Element indices with support node indices appearing more than once.


def __check(mesh, options: Options) -> Result:
    points = mesh.GetPoints()

    locator = vtkIncrementalOctreePointLocator()
    locator.SetTolerance(options.tolerance)
    output = vtkPoints()
    locator.InitPointInsertion(output, points.GetBounds())

    # original ids to/from filtered ids.
    filtered_to_original = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1

    rejected_points = defaultdict(list)
    point_id = reference(0)
    for i in range(points.GetNumberOfPoints()):
        is_inserted = locator.InsertUniquePoint(points.GetPoint(i), point_id)
        if not is_inserted:
            # If it's not inserted, `point_id` contains the node that was already at that location.
            # But in that case, `point_id` is the new numbering in the destination points array.
            # It's more useful for the user to get the old index in the original mesh, so he can look for it in his data.
            logging.debug(
                f"Point {i} at {points.GetPoint(i)} has been rejected, point {filtered_to_original[point_id.get()]} is already inserted."
            )
            rejected_points[point_id.get()].append(i)
        else:
            # If it's inserted, `point_id` contains the new index in the destination array.
            # We store this information to be able to connect the source and destination arrays.
            # original_to_filtered[i] = point_id.get()
            filtered_to_original[point_id.get()] = i

    tmp = []
    for n, ns in rejected_points.items():
        tmp.append((n, *ns))

    # Checking that the support node indices appear only once per element.
    wrong_support_elements = []
    for c in range(mesh.GetNumberOfCells()):
        cell = mesh.GetCell(c)
        num_points_per_cell = cell.GetNumberOfPoints()
        if len({cell.GetPointId(i) for i in range(num_points_per_cell)}) != num_points_per_cell:
            wrong_support_elements.append(c)

    return Result(nodes_buckets=tmp,
                  wrong_support_elements=wrong_support_elements)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
