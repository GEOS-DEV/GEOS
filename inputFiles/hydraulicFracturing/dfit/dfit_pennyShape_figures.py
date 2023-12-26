import os
import sys
import argparse

import dfit_analysis as dfit 
def main():
    """
    Entry point
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-sol", "--solution-file", help=" solution file")
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
    
    dfit.plot_bhp_vs_time( data_dict, "dfit_penny_bhp.png" )

    shutIn_time = 360
    dfit.plot_Gfunction( data_dict, shutIn_time, "dfit_penny_gfunction.png")
    
if __name__ == '__main__':
    main()
