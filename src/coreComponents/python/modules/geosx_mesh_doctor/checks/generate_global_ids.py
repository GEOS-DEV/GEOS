from dataclasses import dataclass
import logging

from vtkmodules.vtkCommonCore import (
    vtkIdTypeArray, )

from . import vtk_utils
from .vtk_utils import (
    VtkOutput,
)


@dataclass(frozen=True)
class Options:
    vtk_output: VtkOutput
    generate_cells_global_ids: bool
    generate_points_global_ids: bool


@dataclass(frozen=True)
class Result:
    info: str


def __build_global_ids(mesh,
                       generate_cells_global_ids: bool,
                       generate_points_global_ids: bool) -> None:
    """
    Adds the global ids for cells and points in place into the mesh instance.
    :param mesh:
    :return: None
    """
    # Building GLOBAL_IDS for points and cells.g GLOBAL_IDS for points and cells.
    # First for points...
    if mesh.GetPointData().GetGlobalIds():
        logging.error("Mesh already has globals ids for points; nothing done.")
    elif generate_points_global_ids:
        point_global_ids = vtkIdTypeArray()
        point_global_ids.SetName("GLOBAL_IDS_POINTS")
        point_global_ids.Allocate(mesh.GetNumberOfPoints())
        for i in range(mesh.GetNumberOfPoints()):
            point_global_ids.InsertNextValue(i)
        mesh.GetPointData().SetGlobalIds(point_global_ids)
    # ... then for cells.
    if mesh.GetCellData().GetGlobalIds():
        logging.error("Mesh already has globals ids for cells; nothing done.")
    elif generate_cells_global_ids:
        cells_global_ids = vtkIdTypeArray()
        cells_global_ids.SetName("GLOBAL_IDS_CELLS")
        cells_global_ids.Allocate(mesh.GetNumberOfCells())
        for i in range(mesh.GetNumberOfCells()):
            cells_global_ids.InsertNextValue(i)
        mesh.GetCellData().SetGlobalIds(cells_global_ids)


def __check(mesh, options: Options) -> Result:
    __build_global_ids(mesh, options.generate_cells_global_ids, options.generate_points_global_ids)
    vtk_utils.write_mesh(mesh, options.vtk_output)
    return Result(info=f"Mesh was written to {options.vtk_output.output}")


def check(vtk_input_file: str, options: Options) -> Result:
    try:
        mesh = vtk_utils.read_mesh(vtk_input_file)
        return __check(mesh, options)
    except BaseException as e:
        logging.error(e)
        return Result(info="Something went wrong.")
