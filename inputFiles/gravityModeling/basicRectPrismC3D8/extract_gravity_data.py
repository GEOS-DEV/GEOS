#!/usr/bin/env python3
"""
Extract gravity data from HDF5 file and convert to ASCII format
"""

import h5py
import numpy as np
import sys

def extract_gravity_data(hdf5_file, output_file):
    """Extract gravity data and coordinates from HDF5 file"""
    
    try:
        with h5py.File(hdf5_file, 'r') as f:
            # Read time data
            time_data = f['/gravity Time'][:]
            
            # Read coordinate data
            coord_x = f['/coordinateX'][:]
            coord_y = f['/coordinateY'][:]
            coord_z = f['/coordinateZ'][:]
            
            # Read gravity data
            gravity_data = f['/scatterData'][:]
            
            print(f"Data shape information:")
            print(f"  Time steps: {time_data.shape}")
            print(f"  Coordinates X: {coord_x.shape}")
            print(f"  Coordinates Y: {coord_y.shape}")
            print(f"  Coordinates Z: {coord_z.shape}")
            print(f"  Gravity data: {gravity_data.shape}")
            
            # Write to ASCII file
            with open(output_file, 'w') as out:
                # Write header
                out.write("# GEOS Gravity Data Export\n")
                out.write(f"# Number of time steps: {time_data.shape[0] if len(time_data.shape) > 0 else 1}\n")
                out.write(f"# Number of stations: {coord_x.shape[1] if len(coord_x.shape) > 1 else coord_x.shape[0]}\n")
                out.write("# Format: Station_ID, X, Y, Z, Time, Gravity\n")
                out.write("#\n")
                
                # Get number of stations
                if len(coord_x.shape) > 1:
                    n_stations = coord_x.shape[1]
                    n_times = time_data.shape[0]
                else:
                    n_stations = coord_x.shape[0]
                    n_times = 1
                
                # Write data for each time step
                for t_idx in range(n_times):
                    if n_times > 1:
                        current_time = time_data[t_idx, 0] if len(time_data.shape) > 1 else time_data[t_idx]
                        current_gravity = gravity_data[t_idx, :]
                        current_x = coord_x[t_idx, :]
                        current_y = coord_y[t_idx, :]
                        current_z = coord_z[t_idx, :]
                    else:
                        current_time = time_data[0, 0] if len(time_data.shape) > 1 else time_data[0]
                        current_gravity = gravity_data[0, :] if len(gravity_data.shape) > 1 else gravity_data
                        current_x = coord_x[0, :] if len(coord_x.shape) > 1 else coord_x
                        current_y = coord_y[0, :] if len(coord_y.shape) > 1 else coord_y
                        current_z = coord_z[0, :] if len(coord_z.shape) > 1 else coord_z
                    
                    out.write(f"# Time step {t_idx + 1}, time = {current_time}\n")
                    
                    # Write station data
                    for station_id in range(n_stations):
                        out.write(f"{station_id:6d} {current_x[station_id]:12.3f} {current_y[station_id]:12.3f} {current_z[station_id]:12.3f} {current_time:12.3f} {current_gravity[station_id]:15.6e}\n")
                    
                    out.write("\n")  # Blank line between time steps
                
            print(f"Successfully exported data to: {output_file}")
            
            # Print summary statistics
            print(f"\nSummary:")
            print(f"  Time range: {np.min(time_data):.3f} to {np.max(time_data):.3f}")
            if len(coord_x.shape) > 1:
                print(f"  X coordinate range: {np.min(coord_x):.3f} to {np.max(coord_x):.3f}")
                print(f"  Y coordinate range: {np.min(coord_y):.3f} to {np.max(coord_y):.3f}")
                print(f"  Z coordinate range: {np.min(coord_z):.3f} to {np.max(coord_z):.3f}")
                print(f"  Gravity range: {np.min(gravity_data):.6e} to {np.max(gravity_data):.6e}")
            else:
                print(f"  X coordinate range: {np.min(current_x):.3f} to {np.max(current_x):.3f}")
                print(f"  Y coordinate range: {np.min(current_y):.3f} to {np.max(current_y):.3f}")
                print(f"  Z coordinate range: {np.min(current_z):.3f} to {np.max(current_z):.3f}")
                print(f"  Gravity range: {np.min(current_gravity):.6e} to {np.max(current_gravity):.6e}")
                
    except Exception as e:
        print(f"Error reading HDF5 file: {e}")
        return False
    
    return True

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_gravity_data.py <input.hdf5> [output.txt]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "gravity_data.txt"
    
    print(f"Extracting gravity data from: {input_file}")
    print(f"Output file: {output_file}")
    
    success = extract_gravity_data(input_file, output_file)
    sys.exit(0 if success else 1)
