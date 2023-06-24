from dataclasses import dataclass
import logging

from typing import (
    Sequence,
)

import numpy

from vtkmodules.vtkCommonDataModel import (
    vtkUnstructuredGrid,
)
from vtkmodules.vtkCommonCore import (
    vtkIdTypeArray,
    vtkPoints,
)
from vtkmodules.vtkIOXML import (
    vtkXMLMultiBlockDataReader,
)
from vtkmodules.util.numpy_support import (
    vtk_to_numpy,
)

# from vtk import *


@dataclass(frozen=True)
class Options:
    tolerance: float
    matrix_name: str
    fracture_name: str
    duplicated_nodes_field_name: str


@dataclass(frozen=True)
class Result:
    # First index is the local index of the fracture mesh.
    # Second is the local index of the matrix mesh.
    # Third is the global index in the matrix mesh.
    errors: Sequence[tuple[int, int, int]]


def __read_multiblock(vtk_input_file: str, matrix_name: str, fracture_name: str) -> (vtkUnstructuredGrid, vtkUnstructuredGrid):
    reader = vtkXMLMultiBlockDataReader()
    reader.SetFileName(vtk_input_file)
    reader.Update()
    multi_block = reader.GetOutput()
    for b in range(multi_block.GetNumberOfBlocks()):
        block_name: str = multi_block.GetMetaData(b).Get(multi_block.NAME())
        if block_name == matrix_name:
            matrix: vtkUnstructuredGrid = multi_block.GetBlock(b)
        if block_name == fracture_name:
            fracture: vtkUnstructuredGrid = multi_block.GetBlock(b)
    assert matrix and fracture
    return matrix, fracture


def __check(vtk_input_file: str, options: Options) -> Result:
    matrix, fracture = __read_multiblock(vtk_input_file, options.matrix_name, options.fracture_name)
    matrix_points: vtkPoints = matrix.GetPoints()
    fracture_points: vtkPoints = fracture.GetPoints()

    duplicated_nodes: vtkIdTypeArray = fracture.GetPointData().GetArray(options.duplicated_nodes_field_name)
    duplicated_nodes = vtk_to_numpy(duplicated_nodes)

    point_ids = vtk_to_numpy(matrix.GetPointData().GetGlobalIds())
    g2l = numpy.ones(len(point_ids), dtype=int)
    for loc, glo in enumerate(point_ids):
        g2l[glo] = loc

    errors = []
    for i, duplicates in enumerate(duplicated_nodes):
        for duplicate in filter(lambda i: i > -1, duplicates):
            p0 = matrix_points.GetPoint(g2l[duplicate])
            p1 = fracture_points.GetPoint(i)
            if numpy.linalg.norm(numpy.array(p1) - numpy.array(p0)) > options.tolerance:
                errors.append((i, g2l[duplicate], duplicate))
    return Result(errors=errors)


def check(vtk_input_file: str, options: Options) -> Result:
    try:
        return __check(vtk_input_file, options)
    except BaseException as e:
        logging.error(e)
        return Result(errors=())


if __name__ == '__main__':
    opt = Options(tolerance=1.e-10, matrix_name="main", fracture_name="fracture", duplicated_nodes_field_name="duplicated_nodes")
    check("/Users/j0436735/CLionProjects/GEOS/inputFiles/tmp_buffer/hi24l/main.vtm", opt)
