import os.path
import logging
import sys

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


def read_mesh(vtk_input_file: str) -> vtkUnstructuredGrid:
    # Testing legacy file format.
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(vtk_input_file)
    if reader.IsFileUnstructuredGrid():
        logging.debug(f"Reading file \"{vtk_input_file}\" using legacy format reader.")
        reader.Update()
        return reader.GetOutput()
    # Now it's xml time!
    reader = vtkXMLUnstructuredGridReader()
    if reader.CanReadFile(vtk_input_file):
        logging.debug(f"Reading file \"{vtk_input_file}\" using XML format reader.")
        reader.SetFileName(vtk_input_file)
        reader.Update()
        return reader.GetOutput()
    logging.critical(f"Could not find the appropriate VTK reader for file \"{vtk_input_file}\". Dying...")
    sys.exit(1)


def write_mesh(mesh: vtkUnstructuredGrid, output: str) -> None:
    """
    Writes the mesh to disk.
    Nothing will be done if the file already exists.
    :param mesh:
    :param output:
    :return: None
    """
    if os.path.exists(output):
        logging.error(f"File \"{output}\" already exists, nothing done.")
        return
    writer = vtkUnstructuredGridWriter()
    writer.SetFileName(output)
    writer.SetInputData(mesh)
    logging.info(f"Writing mesh into file {output}")
    writer.Write()
