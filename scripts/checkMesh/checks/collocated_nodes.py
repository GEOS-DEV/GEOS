from collections import defaultdict
from dataclasses import dataclass
import logging
from typing import Tuple

import numpy

import vtk  # TODO use new pyvtk style + deal with ImportError


@dataclass(frozen=True)
class Options:
    tolerance: float


@dataclass(frozen=True)
class Result:
    nodes_buckets: Tuple[Tuple[int]]


def check(vtk_input_file: str, options: Options) -> Result:
    reader = vtk.vtkXMLUnstructuredGridReader()  # TODO Find a generic way to read the vtk mesh.
    reader.SetFileName(vtk_input_file)
    reader.Update()
    mesh = reader.GetOutput()

    points = mesh.GetPoints()

    locator = vtk.vtkIncrementalOctreePointLocator()
    locator.SetTolerance(options.tolerance)
    output = vtk.vtkPoints()
    locator.InitPointInsertion(output, points.GetBounds())

    # original ids to/from filtered ids.
    # original_to_filtered = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1
    filtered_to_original = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1

    rejected_points = defaultdict(list)
    point_id = vtk.reference(0)
    for i in range(points.GetNumberOfPoints()):
        is_inserted = locator.InsertUniquePoint(points.GetPoint(i), point_id)
        if not is_inserted:
            logging.debug(f"Point {i} at {points.GetPoint(i)} has been rejected, point {filtered_to_original[point_id.get()]} is already inserted.")
            rejected_points[point_id.get()].append(i)
        else:
            # original_to_filtered[i] = point_id.get()
            filtered_to_original[point_id.get()] = i

    tmp = []
    for n, ns in rejected_points.items():
        tmp.append((n, *ns))

    return Result(nodes_buckets=tmp)
