import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from vtk.util.numpy_support import vtk_to_numpy
import xml.etree.ElementTree as ET
import math
from math import sin,cos,tan,exp
from dataclasses import dataclass, asdict
from typing import Iterable, List
import hdf5_wrapper

from vtkmodules.vtkIOLegacy import (
    vtkUnstructuredGridReader,
)
from vtkmodules.vtkIOXML import (
    vtkXMLUnstructuredGridReader,
    vtkXMLMultiBlockDataReader
)

@dataclass(frozen=False)    
class BHPData:
    time: Iterable[float]
    pressure: Iterable[float]

@dataclass(frozen=False)    
class GfunctionData:
    pressure: Iterable[float]
    Gtime: Iterable[float]
    dPdG: Iterable[float]

Pa2psi = 0.000145038
Pa2MPa = 1e-6  

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
    pressure = [element * Pa2MPa for element in pressure]

    return time, pressure
    
def plot_bhp_vs_time( data_dict, figureName ):
    
    fig, ax = plt.subplots(figsize=(24, 18))
    
    for key, bhpData in data_dict.items():
        # Extract the time and pressure data from the Data instance
        time = bhpData.time
        pressure = bhpData.pressure
        
        # Create a plot for each data set
        ax.plot(time, pressure, marker='o', linestyle='-', label=key)
        
    # Customize the plot
    ax.set_title('BHP vs. Time', size=30)
    ax.set_xlabel('Time (s)', size=30)
    ax.set_ylabel('BHP (Mpa)', size=30)
    ax.legend(fontsize="30")

    ax.xaxis.set_tick_params(labelsize=30)
    ax.yaxis.set_tick_params(labelsize=30)

    plt.savefig( figureName, dpi=200 )

def Gfunction(delt_D):
    gfun = 4./3.*( (1.+ delt_D)**1.5 - delt_D**1.5 )
    g0 = 4./3.
    Gfun = 4./math.pi*( gfun - g0 )
    return Gfun

def compute_GFunction(time, pressure, t_shutIn):
    """
    Calculate the G-function for a given set of time and pressure data after a shut-in event.

    Args:
    - time (list or array): A list of time values.
    - pressure (list or array): A list of corresponding pressure values.
    - t_shutIn (float): The shut-in time.

    Returns:
    - dtD (list): A list of dimensionless time values.
    - G_tD (list): The G-function values.
    - P_tD (list): A list of corresponding pressure values.
    - dPdG
    """
    dtD = []
    G_tD = []
    P_tD = [] 
    for t, p in zip(time, pressure):
        if ( t > t_shutIn ):
            tD = (t - t_shutIn)/t_shutIn
            dtD.append(tD)
            G_tD.append( Gfunction(tD) )
            P_tD.append(p)

    dtD = np.array(dtD)
    G_tD = np.array(G_tD)
    P_tD = np.array(P_tD)

    dPdG = np.zeros(len(G_tD)-1)
    GdPdG = np.zeros(len(G_tD)-1)
    for i in range(1, len(dPdG)):
        dPdG[i] = -(P_tD[i+1] - P_tD[i-1])/(G_tD[i+1] - G_tD[i-1])
        GdPdG[i] = G_tD[i] * dPdG[i]        

    return G_tD, P_tD, dPdG

def plot_Gfunction( data_dict, shutIn_time, figureName):
   
   fig, ax = plt.subplots( figsize=(32, 18) )
   ax2 = ax.twinx()

   for key, bhpData in data_dict.items():
       G, p, dPdG  = compute_GFunction(bhpData.time, bhpData.pressure, shutIn_time)
       # Create a plot for each data set
       ax.plot(G, p, marker='o', linestyle='-', label=key)
       ax2.plot(G[0:len(G)-1], dPdG, marker='x', linestyle='-', color="r")
   
   ax.set_title('G-function plot', size=30)
   ax.set_xlabel('G-function', size=30)
   ax.set_ylabel('BHP (psi)', size=30)
   ax2.set_ylabel('dPdG', size=30)
   ax.legend()
   
   ax.xaxis.set_tick_params(labelsize=30)
   ax.yaxis.set_tick_params(labelsize=30)
   ax2.yaxis.set_tick_params(labelsize=30)
   plt.savefig( figureName, dpi=200 )
