from dataclasses import dataclass
import logging

from typing import (
    Collection,
    FrozenSet,
    Iterable,
    Sequence,
    Set,
    Tuple,
)

from tqdm import tqdm
import numpy

from vtkmodules.vtkCommonDataModel import (
    vtkUnstructuredGrid,
    vtkCell,
)
from vtkmodules.vtkCommonCore import (
    vtkPoints,
)
from vtkmodules.vtkIOXML import (
    vtkXMLMultiBlockDataReader,
)
from vtkmodules.util.numpy_support import (
    vtk_to_numpy,
)
from vtk_utils import (
    vtk_iter,
)


@dataclass(frozen=True)
class Options:
    tolerance: float
    matrix_name: str
    fracture_name: str
    collocated_nodes_field_name: str


@dataclass(frozen=True)
class Result:
    # First index is the local index of the fracture mesh.
    # Second is the local index of the matrix mesh.
    # Third is the global index in the matrix mesh.
    errors: Sequence[tuple[int, int, int]]


def __read_multiblock(vtk_input_file: str, matrix_name: str, fracture_name: str) -> Tuple[vtkUnstructuredGrid, vtkUnstructuredGrid]:
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


def format_collocated_nodes(fracture_mesh: vtkUnstructuredGrid) -> Sequence[Iterable[int]]:
    """
    Extract the collocated nodes information from the mesh and formats it in a python way.
    :param fracture_mesh: The mesh of the fracture (with 2d cells).
    :return: An iterable over all the buckets of collocated nodes.
    """
    collocated_nodes: numpy.ndarray = vtk_to_numpy(fracture_mesh.GetPointData().GetArray("collocated_nodes"))
    if len(collocated_nodes.shape) == 1:
        collocated_nodes: numpy.ndarray = collocated_nodes.reshape((collocated_nodes.shape[0], 1))
    return tuple(map(lambda bucket: tuple(sorted(filter(lambda i: i != -1, bucket))), collocated_nodes))


def __check_collocated_nodes_positions(matrix_points: Sequence[Tuple[float, float, float]],
                                       fracture_points: Sequence[Tuple[float, float, float]],
                                       g2l: Sequence[int],
                                       collocated_nodes: Iterable[Iterable[int]]) -> Collection[Tuple[int, Iterable[int], Iterable[Tuple[float, float, float]]]]:
    issues = []
    for li, bucket in enumerate(collocated_nodes):
        matrix_nodes = (fracture_points[li], ) + tuple(map(lambda gi: matrix_points[g2l[gi]], bucket))
        m = numpy.array(matrix_nodes)
        rank: int = numpy.linalg.matrix_rank(m)
        if rank > 1:
            issues.append((li, bucket, tuple(map(lambda gi: matrix_points[g2l[gi]], bucket))))
    return issues


def my_iter(ccc):
    car, cdr = ccc[0], ccc[1:]
    for i in car:
        if cdr:
            for j in my_iter(cdr):
                yield i, *j
        else:
            yield (i, )


def __check_neighbors(matrix: vtkUnstructuredGrid,
                      fracture: vtkUnstructuredGrid,
                      g2l: Sequence[int],
                      collocated_nodes: Sequence[Iterable[int]]):
    fracture_nodes: Set[int] = set()
    for bucket in collocated_nodes:
        for gi in bucket:
            fracture_nodes.add(g2l[gi])
    # For each face of each cell,
    # if all the points of the face are "made" of collocated nodes,
    # then this is a fracture face.
    fracture_faces: Set[FrozenSet[int]] = set()
    for c in range(matrix.GetNumberOfCells()):
        cell: vtkCell = matrix.GetCell(c)
        for f in range(cell.GetNumberOfFaces()):
            face: vtkCell = cell.GetFace(f)
            point_ids = frozenset(vtk_iter(face.GetPointIds()))
            if point_ids <= fracture_nodes:
                fracture_faces.add(point_ids)
    # Finding the cells
    for c in tqdm(range(fracture.GetNumberOfCells()), desc="Finding neighbor cell pairs"):
        cell: vtkCell = fracture.GetCell(c)
        cns: Set[FrozenSet[int]] = set()  # subset of collocated_nodes
        point_ids = frozenset(vtk_iter(cell.GetPointIds()))
        for point_id in point_ids:
            bucket = collocated_nodes[point_id]
            local_bucket = frozenset(map(g2l.__getitem__, bucket))
            cns.add(local_bucket)
        found = 0
        tmp = tuple(map(tuple, cns))
        for node_combinations in my_iter(tmp):
            f = frozenset(node_combinations)
            if f in fracture_faces:
                found += 1
        if found != 2:
            logging.warning(f"Something went wrong since we should have found 2 fractures faces (we found {found}) for collocated nodes {cns}.")


def __check(vtk_input_file: str, options: Options) -> Result:
    matrix, fracture = __read_multiblock(vtk_input_file, options.matrix_name, options.fracture_name)
    matrix_points: vtkPoints = matrix.GetPoints()
    fracture_points: vtkPoints = fracture.GetPoints()

    collocated_nodes: Sequence[Iterable[int]] = format_collocated_nodes(fracture)
    assert matrix.GetPointData().GetGlobalIds() and matrix.GetCellData().GetGlobalIds() and \
           fracture.GetPointData().GetGlobalIds() and fracture.GetCellData().GetGlobalIds()

    point_ids = vtk_to_numpy(matrix.GetPointData().GetGlobalIds())
    g2l = numpy.ones(len(point_ids), dtype=int) * -1
    for loc, glo in enumerate(point_ids):
        g2l[glo] = loc
    g2l.flags.writeable = False

    issues = __check_collocated_nodes_positions(vtk_to_numpy(matrix.GetPoints().GetData()),
                                                vtk_to_numpy(fracture.GetPoints().GetData()),
                                                g2l,
                                                collocated_nodes)
    assert len(issues) == 0

    __check_neighbors(matrix, fracture, g2l, collocated_nodes)

    errors = []
    for i, duplicates in enumerate(collocated_nodes):
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
