from dataclasses import dataclass
import logging
from typing import List

from vtkmodules.vtkFiltersGeneral import (
    vtkCellValidator
)
from vtk.util.numpy_support import (
    vtk_to_numpy,
)

from . import vtk_utils


@dataclass(frozen=True)
class Options:
    tolerance: float


@dataclass(frozen=True)
class Result:
    jumbled_elements: List[int]


def __check(mesh, options: Options) -> Result:
    valid = 0x0
    wrong_number_of_points = 0x01
    intersecting_edges = 0x02
    intersecting_faces = 0x04
    non_contiguous_edges = 0x08
    non_convex = 0x10
    faces_are_oriented_incorrectly = 0x20

    f = vtkCellValidator()
    f.SetTolerance(options.tolerance)

    f.SetInputData(mesh)
    f.Update()
    output = f.GetOutput()

    cd = output.GetCellData()
    for i in range(cd.GetNumberOfArrays()):
        if cd.GetArrayName(i) == "ValidityState":  # TODO change name?
            validity = vtk_to_numpy(cd.GetArray(i))
    # TODO assert validity exists
    jumbled_elements: List[int] = []
    for i, v in enumerate(validity):
        if v & intersecting_faces:
            jumbled_elements.append(i)
    return Result(jumbled_elements=jumbled_elements)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
