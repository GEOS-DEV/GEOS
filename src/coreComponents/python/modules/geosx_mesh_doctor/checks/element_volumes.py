from dataclasses import dataclass
import logging
from typing import List, Tuple

from vtkmodules.vtkFiltersVerdict import (
    vtkCellSizeFilter
)
from vtk.util.numpy_support import (
    vtk_to_numpy,
)


from . import vtk_utils


@dataclass(frozen=True)
class Options:
    min_volume: float


@dataclass(frozen=True)
class Result:
    element_volumes: List[Tuple[int, float]]


def __check(mesh, options: Options) -> Result:
    f = vtkCellSizeFilter()

    f.ComputeAreaOff()
    f.ComputeLengthOff()
    f.ComputeSumOff()
    f.ComputeVertexCountOff()
    f.ComputeVolumeOn()
    volume_array_name = "__MESH_DOCTOR_VOLUME"
    # TODO assert it does not exist
    f.SetVolumeArrayName(volume_array_name)

    f.SetInputData(mesh)
    f.Update()
    output = f.GetOutput()
    cd = output.GetCellData()
    for i in range(cd.GetNumberOfArrays()):
        if cd.GetArrayName(i) == volume_array_name:
            volume = vtk_to_numpy(cd.GetArray(i))
    # TODO assert volume exists
    small_volumes: List[Tuple[int, float]] = []
    for i, v in enumerate(volume):
        if v < options.min_volume:
            small_volumes.append((i, v))
    return Result(element_volumes=small_volumes)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
