from dataclasses import dataclass
import logging
from typing import List, Dict, Set

from vtkmodules.vtkCommonCore import (
    vtkIdList,
)

from . import vtk_utils
from .vtk_utils import (
    to_vtk_id_list,
)


@dataclass(frozen=True)
class Options:
    output: str
    cell_type_to_ordering: Dict[int, List[int]]


@dataclass(frozen=True)
class Result:
    output: str
    unchanged_cell_types: Set[int]


def __check(mesh, options: Options):
    cell_type_to_ordering: Dict[int, List[int]] = options.cell_type_to_ordering
    unchanged_cell_types = []
    output_mesh = mesh.NewInstance()  # keeping the same instance type.
    output_mesh.CopyStructure(mesh)
    output_mesh.CopyAttributes(mesh)

    cells = output_mesh.GetCells()
    for cell_idx in range(output_mesh.GetNumberOfCells()):
        support_point_ids = vtkIdList()
        cells.GetCellAtId(cell_idx, support_point_ids)
        cell_type = output_mesh.GetCell(cell_idx).GetCellType()
        new_ordering = cell_type_to_ordering.get(cell_type)
        if new_ordering:
            tmp = []
            for i, v in enumerate(new_ordering):
                tmp.append(support_point_ids.GetId(new_ordering[i]))
            new_support_point_ids = to_vtk_id_list(tmp)
            cells.ReplaceCellAtId(cell_idx, new_support_point_ids)
        else:
            unchanged_cell_types.insert(cell_type)
    is_written_error = vtk_utils.write_mesh(output_mesh, options.output)
    return Result(output=options.output if not is_written_error else "",
                  unchanged_cell_types=set(unchanged_cell_types))


def check(vtk_input_file: str, options: Options) -> Result:
    mesh = vtk_utils.read_mesh(vtk_input_file)
    return __check(mesh, options)
