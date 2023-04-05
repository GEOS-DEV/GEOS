from dataclasses import dataclass
from typing import List, Tuple
import uuid

from vtkmodules.vtkCommonDataModel import (
    VTK_HEXAHEDRON,
    VTK_PYRAMID,
    VTK_TETRA,
    VTK_WEDGE,
)
from vtkmodules.vtkFiltersVerdict import (
    vtkCellSizeFilter,
    vtkMeshQuality,
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
    cs = vtkCellSizeFilter()

    cs.ComputeAreaOff()
    cs.ComputeLengthOff()
    cs.ComputeSumOff()
    cs.ComputeVertexCountOff()
    cs.ComputeVolumeOn()
    volume_array_name = "__MESH_DOCTOR_VOLUME-" + str(uuid.uuid4())  # Making the name unique
    cs.SetVolumeArrayName(volume_array_name)

    cs.SetInputData(mesh)
    cs.Update()

    mq = vtkMeshQuality()

    mq.SetTetQualityMeasureToVolume()
    mq.SetPyramidQualityMeasureToVolume()
    mq.SetWedgeQualityMeasureToVolume()
    mq.SetHexQualityMeasureToVolume()

    mq.SetInputData(mesh)
    mq.Update()

    volume = cs.GetOutput().GetCellData().GetArray(volume_array_name)
    quality = mq.GetOutput().GetCellData().GetArray("Quality")  # Name is imposed by vtk.

    assert volume is not None
    assert quality is not None
    volume = vtk_to_numpy(volume)
    quality = vtk_to_numpy(quality)
    small_volumes: List[Tuple[int, float]] = []
    for i, pack in enumerate(zip(volume, quality)):
        v, q = pack
        vol = q if mesh.GetCellType(i) in (VTK_HEXAHEDRON, VTK_PYRAMID, VTK_TETRA, VTK_WEDGE) else v
        if vol < options.min_volume:
            small_volumes.append((i, vol))
    return Result(element_volumes=small_volumes)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
