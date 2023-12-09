from dataclasses import dataclass
import logging
from typing import Sequence, Iterable

import numpy

from vtkmodules.vtkCommonCore import (
    vtkPoints,
)
from vtkmodules.vtkCommonDataModel import (
    VTK_HEXAHEDRON,
    vtkCellArray,
    vtkHexahedron,
    vtkRectilinearGrid,
    vtkUnstructuredGrid,
)
from vtkmodules.util.numpy_support import (
    numpy_to_vtk,
)

from . import vtk_utils
from .vtk_utils import (
    VtkOutput,
)

from .generate_global_ids import __build_global_ids


@dataclass(frozen=True)
class Result:
    info: str


@dataclass(frozen=True)
class FieldInfo:
    name: str
    dimension: int
    support: str


@dataclass(frozen=True)
class Options:
    vtk_output: VtkOutput
    generate_cells_global_ids: bool
    generate_points_global_ids: bool
    xs: Sequence[float]
    ys: Sequence[float]
    zs: Sequence[float]
    nxs: Sequence[int]
    nys: Sequence[int]
    nzs: Sequence[int]
    fields: Iterable[FieldInfo]


@dataclass(frozen=True)
class XYZ:
    x: numpy.ndarray
    y: numpy.ndarray
    z: numpy.ndarray


def build_rectilinear_blocks_mesh(xyzs: Iterable[XYZ]) -> vtkUnstructuredGrid:
    """
    Builds an unstructured vtk grid from the `xyzs` blocks. Kind of InternalMeshGenerator.
    :param xyzs: The blocks.
    :return: The unstructured mesh, even if it's topologically structured.
    """
    rgs = []
    for xyz in xyzs:
        rg = vtkRectilinearGrid()
        rg.SetDimensions(len(xyz.x), len(xyz.y), len(xyz.z))
        rg.SetXCoordinates(numpy_to_vtk(xyz.x))
        rg.SetYCoordinates(numpy_to_vtk(xyz.y))
        rg.SetZCoordinates(numpy_to_vtk(xyz.z))
        rgs.append(rg)

    num_points = sum(map(lambda r: r.GetNumberOfPoints(), rgs))
    num_cells = sum(map(lambda r: r.GetNumberOfCells(), rgs))

    points = vtkPoints()
    points.Allocate(num_points)
    for rg in rgs:
        for i in range(rg.GetNumberOfPoints()):
            points.InsertNextPoint(rg.GetPoint(i))

    cell_types = [VTK_HEXAHEDRON] * num_cells
    cells = vtkCellArray()
    cells.AllocateExact(num_cells, num_cells * 8)

    m = (0, 1, 3, 2, 4, 5, 7, 6)  # VTK_VOXEL and VTK_HEXAHEDRON do not share the same ordering.
    offset = 0
    for rg in rgs:
        for i in range(rg.GetNumberOfCells()):
            c = rg.GetCell(i)
            new_cell = vtkHexahedron()
            for j in range(8):
                new_cell.GetPointIds().SetId(j, offset + c.GetPointId(m[j]))
            cells.InsertNextCell(new_cell)
        offset += rg.GetNumberOfPoints()

    mesh = vtkUnstructuredGrid()
    mesh.SetPoints(points)
    mesh.SetCells(cell_types, cells)

    return mesh


def __add_fields(mesh: vtkUnstructuredGrid, fields: Iterable[FieldInfo]) -> vtkUnstructuredGrid:
    for field_info in fields:
        if field_info.support == "CELLS":
            data = mesh.GetCellData()
            n = mesh.GetNumberOfCells()
        elif field_info.support == "POINTS":
            data = mesh.GetPointData()
            n = mesh.GetNumberOfPoints()
        array = numpy.ones((n, field_info.dimension), dtype=float)
        vtk_array = numpy_to_vtk(array)
        vtk_array.SetName(field_info.name)
        data.AddArray(vtk_array)
    return mesh


def __build(options: Options):
    def build_coordinates(positions, num_elements):
        result = []
        it = zip(zip(positions, positions[1:]), num_elements)
        try:
            coords, n = next(it)
            while True:
                start, stop = coords
                end_point = False
                tmp = numpy.linspace(start=start, stop=stop, num=n+end_point, endpoint=end_point)
                coords, n = next(it)
                result.append(tmp)
        except StopIteration:
            end_point = True
            tmp = numpy.linspace(start=start, stop=stop, num=n+end_point, endpoint=end_point)
            result.append(tmp)
        return numpy.concatenate(result)
    x = build_coordinates(options.xs, options.nxs)
    y = build_coordinates(options.ys, options.nys)
    z = build_coordinates(options.zs, options.nzs)
    cube = build_rectilinear_blocks_mesh((XYZ(x, y, z),))
    cube = __add_fields(cube, options.fields)
    __build_global_ids(cube, options.generate_cells_global_ids, options.generate_points_global_ids)
    return cube


def __check(options: Options) -> Result:
    output_mesh = __build(options)
    vtk_utils.write_mesh(output_mesh, options.vtk_output)
    return Result(info=f"Mesh was written to {options.vtk_output.output}")


def check(vtk_input_file: str, options: Options) -> Result:
    try:
        return __check(options)
    except BaseException as e:
        logging.error(e)
        return Result(info="Something went wrong.")
