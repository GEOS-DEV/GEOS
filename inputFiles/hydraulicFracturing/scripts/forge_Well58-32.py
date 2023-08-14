import hdf5_wrapper
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import argparse
from vtk.util.numpy_support import vtk_to_numpy
import xml.etree.ElementTree as ET

from vtkmodules.vtkIOLegacy import (
    vtkUnstructuredGridReader,
)
from vtkmodules.vtkIOXML import (
    vtkXMLUnstructuredGridReader,
    vtkXMLMultiBlockDataReader
)

def get_attribute_array(unstructured_grid, attribute_name):
     # Get the point data
    cell_data = unstructured_grid.GetCellData()

    # Find the index of the attribute by name
    attribute_array = cell_data.GetArray(attribute_name)

    numpy_array = vtk_to_numpy(attribute_array)

    return numpy_array


def read_attributes_from_vtu(file_path, attribute_names):
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_path)
    reader.Update()
    unstructured_grid = reader.GetOutput()

    attibutes_array = {}
    for attribute_name in attribute_names:
        attibutes_array[attribute_name] = get_attribute_array( unstructured_grid, attribute_name)

    return attibutes_array

def readDataFromHDF5( file ):
    data = hdf5_wrapper.hdf5_wrapper(file).get_copy()

    time = data["pressure Time"]
    pressure = data["pressure source"]

    return time, pressure

def readDataFromVTK( file ):
    
    attributes_names = ['pressure', 'elementCenter']

    # Specify the path to the .pvd file
  
    base_directory = os.path.dirname(file)
    vtk_files_dir  = os.path.splitext(os.path.basename(file))[0]
    base_directory = os.path.join(base_directory, vtk_files_dir) 

    tree = ET.parse(file)
    root = tree.getroot()

    time = []
    time_steps = []
    pressure = []
    # Iterate through the Collection elements
    for collection in root.findall("Collection"):
        # Iterate through the DataSet elements within the Collection 
        for dataset in collection.findall("DataSet"):
            time.append( int( dataset.get("timestep") ) )
            vmt_file_path = dataset.get("file")
            time_steps.append( int ( vmt_file_path.split('/')[-1].split('.')[0] ) )
    
    print( time )
    print( time_steps )
    # Loop through time step indices
    for time_step in time_steps:
        vtu_dir_path = os.path.join(base_directory, f"{time_step:06d}/mesh1/Level0/Fracture")

        # List all files in the directory
        all_files = os.listdir(vtu_dir_path)

        # Filter files with the .vtu extension
        vtu_files = [file for file in all_files if file.endswith(".vtu")]

        for vtu_file in vtu_files:
            vtu_file_path = os.path.join(vtu_dir_path, vtu_file)   
       
            attributes_array = read_attributes_from_vtu( vtu_file_path, attributes_names )

            injection_location = [100.0, 98.5, 100.5]
            index = np.where(np.all(attributes_array['elementCenter'] == injection_location, axis=1))
            if index[0].size > 0:
                pressure.append( attributes_array['pressure'][index].tolist() )
    
    # Convert the nested list to a single list of scalars
    pressure = [item for pressure in pressure for item in pressure]
    print(pressure)

    return time, pressure
    

def plot(time, pressure):
    fig, ax = plt.subplots(1, 1, figsize=(24, 18))
        # Plot your data using the ax object
    ax.plot(time, pressure, marker="x")

    # Customize the plot as needed
    ax.set_xlabel('Time')
    ax.set_ylabel('Pressure')
    ax.set_title('Pressure over Time')
    plt.savefig( "forge_58-32.png", dpi=200 )


def main():
    """
    Entry point
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-hd5", "--hdf5-file", help="solution file")
    parser.add_argument("-pvd", "--pvd-file", help="pvd directory")
    args = parser.parse_args()
    
    if args.pvd_file is not None:
        time, pressure = readDataFromVTK( args.pvd_file )
    elif args.hdf5_file is not None:
        time, pressure = readDataFromHDF5( args.hdf5_file )
    else:
        raise Exception("One of --pvd-file or --hdf5-file must be provided.")
    
    plot(time, pressure)
     
if __name__ == '__main__':
    main()
