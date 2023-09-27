import os
import importlib.util
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from typing import Iterable, List
import hdf5_wrapper

# Get the current directory of the script
current_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the relative path to the external module
relative_path = os.path.join(current_dir, '../dfit/dfit_analysis.py')

# Use importlib to import the module
module_name = 'dfit_analysis'
spec = importlib.util.spec_from_file_location(module_name, relative_path)
dfit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dfit)

def readFieldData( file ):
    data = hdf5_wrapper.hdf5_wrapper(file).get_copy()["data"]

    time     = np.array(data['block0_values']).flatten().tolist()
    pressure = np.array(data['block1_values'])[:, 0].flatten().tolist()

    return time, pressure 

def main():
    """
    Entry point
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-sol", "--solution-file", help=" solution file")
    parser.add_argument("-fd", "--field-data-file", help="h5 file with field data"  )
    args = parser.parse_args()
    
    data_dict = {}
    if args.solution_file is not None:
        if args.solution_file:
            root, extension = os.path.splitext( args.solution_file )
            if extension == ".pvd":
                time, pressure = dfit.readDataFromVTK( args.solution_file )
            elif extension == ".h5":
                time, pressure = dfit.readDataFromHDF5( args.solution_file )
            else:
                raise Exception("Unsupported solution file extension")
        data_dict["geos"] = dfit.BHPData(time=time, pressure=pressure)    
    if args.field_data_file:
        field_time, field_pressure = readFieldData( args.field_data_file )
        data_dict["field"] = dfit.BHPData(time=field_time, pressure=field_pressure)
    
    dfit.plot_bhp_vs_time( data_dict, "forge_58-32_bhp.png" )

    shutIn_time = 400
    dfit.plot_Gfunction( data_dict, shutIn_time, "forge_58-32_gfunction.png" )
    
     
if __name__ == '__main__':
    main()
