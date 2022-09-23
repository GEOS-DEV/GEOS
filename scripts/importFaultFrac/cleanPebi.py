import logging
import numpy
import vtk

from collections import defaultdict

def main():
    # vtk_input_file = "/docker-exchange/PEBI_10InjectorWells_ResMin100_TRI_14092022.vtu"
    vtk_input_file = "/docker-exchange/PEBI_10InjectorWells_ResMin100_14092022.vtu"

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtk_input_file)
    reader.Update()
    mesh = reader.GetOutput()

    points = mesh.GetPoints()

    locator = vtk.vtkIncrementalOctreePointLocator()
    locator.SetTolerance(1.e-12)
    output = vtk.vtkPoints()
    locator.InitPointInsertion(output, points.GetBounds())

    point_mapping = numpy.ones(points.GetNumberOfPoints(), dtype=int) * -1

    rejected_points = defaultdict(list)
    point_id = vtk.reference(0)
    for i in range(points.GetNumberOfPoints()):
        is_inserted = locator.InsertUniquePoint(points.GetPoint(i), point_id)
        if not is_inserted:
            print(f"Point {i}, {points.GetPoint(i)} has been rejected, {point_id.get()} already exists.")
            rejected_points[point_id.get()].append(i)
        point_mapping[i] = point_id.get()  # Duplicated nodes are considered mapped too.

    # Finding the cells with duplicated nodes
    all_duplicated_nodes = []  # All the nodes that have at least one duplicated node.
    for k, v in rejected_points.items():
        all_duplicated_nodes.append(k)
        all_duplicated_nodes.extend(v)
    all_duplicated_nodes = set(all_duplicated_nodes)

    # Building GLOBAL_IDS for points and cells.
    # First for points...
    point_global_ids = vtk.vtkIntArray()
    point_global_ids.SetName("GLOBAL_IDS")
    point_global_ids.Allocate(points.GetNumberOfPoints())
    for i in range(mesh.GetNumberOfPoints()):
        point_global_ids.InsertNextValue(i)
    mesh.GetPointData().AddArray(point_global_ids)
    # ... then for cells.
    cells_global_ids = vtk.vtkIntArray()
    cells_global_ids.SetName("GLOBAL_IDS")
    cells_global_ids.Allocate(mesh.GetNumberOfCells())
    for i in range(mesh.GetNumberOfCells()):
        cells_global_ids.InsertNextValue(i)
    mesh.GetCellData().AddArray(cells_global_ids)

    # writer = vtk.vtkUnstructuredGridWriter()
    # writer.SetFileTypeToASCII()
    # writer.SetFileName("/docker-exchange/pebi_with_global_ids.vtu")
    # writer.SetInputData(mesh)
    # writer.Write()



if __name__ == '__main__':
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s', level=logging.INFO)
    logging.info("Cleaning PEBI mesh.")
    main()
