from dataclasses import dataclass
import logging
from typing import List, Dict, Set

from vtkmodules.vtkCommonCore import (
    vtkIdList,
)
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
)


from . import vtk_utils

def to_vtk_id_list(data):  # TODO move to utility (and deduplicate)
    result = vtkIdList()
    result.Allocate(len(data))
    for d in data:
        result.InsertNextId(d)
    return result


@dataclass(frozen=True)
class Options:
    output: str
    cell_type_to_ordering: Dict[int, List[int]]


@dataclass(frozen=True)
class Result:
    output: str
    unchanged_cell_types: Set[int]


def copy_fields(input_mesh, output_mesh, copy_field_data=True, copy_point_data=True, copy_cell_data=True):  # TODO move to vtk_utils
    """
    Prefer output_mesh.CopyAttributes(input_mesh) if you want all the fields to be copied
    :param input_mesh:
    :param output_mesh:
    :param copy_field_data:
    :param copy_point_data:
    :param copy_cell_data:
    :return:
    """
    def cp(input_data, output_data):
        num_arrays = input_data.GetNumberOfArrays()
        output_data.AllocateArrays(num_arrays)
        for i in range(num_arrays):
            output_data.AddArray(input_data.GetAbstractArray(i))
        output_data.SetGlobalIds(input_data.GetGlobalIds())  # TODO double check this
    if copy_field_data:
        output_mesh.SetFieldData(input_mesh.GetFieldData())
    if copy_cell_data:
        cp(input_mesh.GetCellData(), output_mesh.GetCellData())
    if copy_point_data:
        cp(input_mesh.GetPointData(), output_mesh.GetPointData())


def __check(mesh, options: Options):
    cell_type_to_ordering: Dict[int, List[int]] = options.cell_type_to_ordering
    unchanged_cell_types = set()
    output_mesh = mesh.NewInstance()  # keeping the same instance type.
    output_mesh.SetPoints(mesh.GetPoints())  # Keeping the same points, obviously.

    new_cells = vtkCellArray()
    new_cells.DeepCopy(mesh.GetCells())

    for cell_idx in range(mesh.GetNumberOfCells()):
        support_point_ids = vtkIdList()
        new_cells.GetCellAtId(cell_idx, support_point_ids)
        cell_type = mesh.GetCell(cell_idx).GetCellType()
        new_ordering = cell_type_to_ordering.get(cell_type)
        if new_ordering:
            tmp = []
            for i, v in enumerate(new_ordering):
                tmp.append(support_point_ids.GetId(new_ordering[i]))
            new_support_point_ids = to_vtk_id_list(tmp)
            new_cells.ReplaceCellAtId(cell_idx, new_support_point_ids)
        else:
            unchanged_cell_types.add(cell_type)

    output_mesh.SetCells(mesh.GetCellTypesArray(), new_cells)  # The cell types are unchanged; we reuse the old cell types!
    output_mesh.CopyAttributes(mesh)
    is_written_error = vtk_utils.write_mesh(output_mesh, options.output)
    return Result(output=options.output if not is_written_error else "",
                  unchanged_cell_types=unchanged_cell_types)


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
