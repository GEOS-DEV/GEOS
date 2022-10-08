from dataclasses import dataclass
import logging

from vtkmodules.vtkCommonCore import (
    vtkIdTypeArray,
)

from . import vtk_utils


@dataclass(frozen=True)
class Options:
    output: str


@dataclass(frozen=True)
class Result:
    info: str


def __build_global_ids(mesh) -> None:
    """
    Adds the global ids for cells and points in place into the mesh instance.
    :param mesh:
    :return: None
    """
    # Building GLOBAL_IDS for points and cells.g GLOBAL_IDS for points and cells.
    # First for points...
    if mesh.GetPointData().GetGlobalIds():
        logging.error("Mesh already has globals ids for points; nothing done.")
    else:
        point_global_ids = vtkIdTypeArray()
        point_global_ids.SetName("GLOBAL_IDS_POINTS")
        point_global_ids.Allocate(mesh.GetNumberOfPoints())
        for i in range(mesh.GetNumberOfPoints()):
            point_global_ids.InsertNextValue(i)
        mesh.GetPointData().SetGlobalIds(point_global_ids)
    # ... then for cells.
    if mesh.GetCellData().GetGlobalIds():
        logging.error("Mesh already has globals ids for cells; nothing done.")
    else:
        cells_global_ids = vtkIdTypeArray()
        cells_global_ids.SetName("GLOBAL_IDS_CELLS")
        cells_global_ids.Allocate(mesh.GetNumberOfCells())
        for i in range(mesh.GetNumberOfCells()):
            cells_global_ids.InsertNextValue(i)
        mesh.GetCellData().SetGlobalIds(cells_global_ids)


def __check(mesh, options: Options) -> Result:
    __build_global_ids(mesh)
    vtk_utils.write_mesh(mesh, options.output)
    return Result(info=f"Mesh was written to {options.output}")


def check(vtk_input_file: str, options: Options) -> Result:
    try:
        mesh = vtk_utils.read_mesh(vtk_input_file)
        return __check(mesh, options)
    except BaseException as e:
        logging.error(e)
        return Result(info="Something went wrong")
