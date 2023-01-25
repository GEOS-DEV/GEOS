import os.path
import logging
import sys

from vtkmodules.vtkCommonCore import (
    vtkIdList,
)
from vtkmodules.vtkCommonDataModel import (
    vtkUnstructuredGrid,
)
from vtkmodules.vtkIOLegacy import (
    vtkUnstructuredGridWriter,
    vtkUnstructuredGridReader,
)
from vtkmodules.vtkIOXML import (
    vtkXMLUnstructuredGridReader,
)

def to_vtk_id_list(data):
    result = vtkIdList()
    result.Allocate(len(data))
    for d in data:
        result.InsertNextId(d)
    return result


def get_cell_field_by_name(mesh, field_name):
    cd = mesh.GetCellData()
    for i in range(cd.GetNumberOfArrays()):
        if cd.GetArrayName(i) == field_name:
            return cd.GetArray(i)


def __read_vtk(vtk_input_file: str) -> vtkUnstructuredGrid|None:
    reader = vtkUnstructuredGridReader()
    logging.info(f"Testing file format \"{vtk_input_file}\" using legacy format reader...")
    reader.SetFileName(vtk_input_file)
    if reader.IsFileUnstructuredGrid():
        logging.info(f"Reader matches. Reading file \"{vtk_input_file}\" using legacy format reader.")
        reader.Update()
        return reader.GetOutput()
    else:
        logging.info("Reader did not match the input file format.")
        return None


def __read_vtu(vtk_input_file: str) -> vtkUnstructuredGrid|None:
    reader = vtkXMLUnstructuredGridReader()
    logging.info(f"Testing file format \"{vtk_input_file}\" using XML format reader...")
    if reader.CanReadFile(vtk_input_file):
        reader.SetFileName(vtk_input_file)
        logging.info(f"Reader matches. Reading file \"{vtk_input_file}\" using XML format reader.")
        reader.Update()
        return reader.GetOutput()
    else:
        logging.info("Reader did not match the input file format.")
        return None


def read_mesh(vtk_input_file: str) -> vtkUnstructuredGrid:
    file_extension = os.path.splitext(vtk_input_file)[-1]
    extension_to_reader = {".vtk": __read_vtk,
                           ".vtu": __read_vtu}
    # Testing first the reader that should match
    if file_extension in extension_to_reader:
        output_mesh = extension_to_reader.pop(file_extension)(vtk_input_file)
        if output_mesh:
            return output_mesh
    # If it does not match, then test all the others.
    for reader in extension_to_reader.values():
        output_mesh = reader(vtk_input_file)
        if output_mesh:
            return output_mesh
    # No reader did work. Dying.
    logging.critical(f"Could not find the appropriate VTK reader for file \"{vtk_input_file}\". Dying...")
    sys.exit(1)


def write_mesh(mesh: vtkUnstructuredGrid, output: str) -> int:
    """
    Writes the mesh to disk.
    Nothing will be done if the file already exists.
    :param mesh:
    :param output:
    :return: None
    """
    if os.path.exists(output):
        logging.error(f"File \"{output}\" already exists, nothing done.")
        return 1
    writer = vtkUnstructuredGridWriter()
    writer.SetFileName(output)
    writer.SetInputData(mesh)
    logging.info(f"Writing mesh into file \"{output}\"")
    return 0 if writer.Write() else 2  # the Write member function return 1 in case of success, 0 otherwise.
