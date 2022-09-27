import logging
from collections import defaultdict

import numpy

import vtk


def main():
    # vtk_input_file = "/docker-exchange/PEBI_10InjectorWells_ResMin100_TRI_14092022.vtu"
    vtk_input_file = "/docker-exchange/PEBI_10InjectorWells_ResMin100_22092022.vtu"

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtk_input_file)
    reader.Update()
    mesh = reader.GetOutput()

    points = mesh.GetPoints()

    locator = vtk.vtkIncrementalOctreePointLocator()
    locator.SetTolerance(1.e-16)
    output = vtk.vtkPoints()
    locator.InitPointInsertion(output, points.GetBounds())

    # point_mapping = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1
    # original ids to/from filtered ids.
    orig_to_filt = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1
    filt_to_orig = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1

    rejected_points = defaultdict(list)
    point_id = vtk.reference(0)
    for i in range(points.GetNumberOfPoints()):
        is_inserted = locator.InsertUniquePoint(points.GetPoint(i), point_id)
        if not is_inserted:
            print(f"Original point {i}, {points.GetPoint(i)} has been rejected, filtered point {point_id.get()} already exists.")
            rejected_points[point_id.get()].append(i)
        else:
            orig_to_filt[i] = point_id.get()
            filt_to_orig[point_id.get()] = i
        # point_mapping[i] = point_id.get()  # Inserted nodes are considered mapped too.

    # Finding the cells with duplicated nodes
    all_duplicated_nodes = []  # All the nodes that have at least one duplicated node.
    for k, v in rejected_points.items():
        all_duplicated_nodes.append(k)
        all_duplicated_nodes.extend(v)
    all_duplicated_nodes = set(all_duplicated_nodes)

    # Building GLOBAL_IDS for points and cells.
    # First for points...
    point_global_ids = vtk.vtkIdTypeArray()
    point_global_ids.SetName("GLOBAL_IDS_POINTS")
    point_global_ids.Allocate(mesh.GetNumberOfPoints())
    for i in range(mesh.GetNumberOfPoints()):
        point_global_ids.InsertNextValue(i)
    mesh.GetPointData().SetGlobalIds(point_global_ids)
    # ... then for cells.
    cells_global_ids = vtk.vtkIdTypeArray()
    cells_global_ids.SetName("GLOBAL_IDS_CELLS")
    cells_global_ids.Allocate(mesh.GetNumberOfCells())
    for i in range(mesh.GetNumberOfCells()):
        cells_global_ids.InsertNextValue(i)
    mesh.GetCellData().SetGlobalIds(cells_global_ids)

    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileTypeToASCII()
    writer.SetFileName("/docker-exchange/pebi_with_global_ids.vtu")
    writer.SetInputData(mesh)
    writer.Write()


if __name__ == '__main__':
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s', level=logging.INFO)
    logging.info("Cleaning PEBI mesh.")
    main()
