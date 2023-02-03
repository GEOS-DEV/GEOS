from dataclasses import dataclass
from typing import List, Tuple
import uuid

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
    volume_array_name = "__MESH_DOCTOR_VOLUME-" + str(uuid.uuid4())  # Making the name unique
    f.SetVolumeArrayName(volume_array_name)

    f.SetInputData(mesh)
    f.Update()
    output = f.GetOutput()
    volume = vtk_utils.get_cell_field_by_name(output, volume_array_name)
    assert volume is not None
    volume = vtk_to_numpy(volume)
    small_volumes: List[Tuple[int, float]] = []
    for i, v in enumerate(volume):
        if v < options.min_volume:
            small_volumes.append((i, v))
    return Result(element_volumes=small_volumes)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
