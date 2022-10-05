import os.path
import logging

import vtk


def read_mesh(vtk_input_file: str) -> vtk.vtkUnstructuredGrid:
    reader = vtk.vtkXMLUnstructuredGridReader()  # TODO Find a generic way to read the vtk mesh.
    reader.SetFileName(vtk_input_file)
    reader.Update()
    return reader.GetOutput()


def write_mesh(mesh: vtk.vtkUnstructuredGrid, output: str) -> None:
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
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(output)
    writer.SetInputData(mesh)
    logging.info(f"Writing mesh into file {output}")
    writer.Write()


