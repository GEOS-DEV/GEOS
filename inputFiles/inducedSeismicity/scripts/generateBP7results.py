import datetime
import numpy as np
import os
import argparse

MODELER = "Matteo Cusini, Lawrence Livermore National Laboratory"
CODE  = "GEOS"
# Predefined fault station locations (x2, x3 in meters)
fault_stations_locations = [
    ("0 m", "0 m"),
    ("-100 m", "0 m"),
    ("0 m", "100 m"),
    ("100 m", "0 m"),
    ("0 m", "-100 m"),
    ("-100 m", "-100 m"),
    ("-100 m", "100 m"),
    ("100 m", "-100 m"),
    ("100 m", "100 m"),
    ("-300 m", "0 m"),
    ("0 m", "300 m"),
    ("300 m", "0 m"),
    ("0 m", "-300 m"),
]


class FaultStationFile:
    """Handles metadata and file writing for a fault station dataset."""

    def __init__(self, station):
        """
        Initialize the FaultFile with metadata constructed from a FaultStation.
        
        :param station: FaultStation instance containing the numerical data and parameters.
        """
        self.station = station  # Reference to the station for data access
        self.header = self._construct_header()
        self.column_headers = [
            "t",
            "slip_2",
            "slip_3",
            "slip_rate_2 ",
            "Slip rate 3 ",
            "shear_stress_2",
            "shear_stress_3",
            "state"
        ]

    def _construct_header(self):
        """Construct the file header dynamically based on station parameters."""
        return {
            "problem": "SEAS Benchmark BP7-QD-A",
            "code": f'{CODE}',
            "version": "1.0",
            "modeler": f'{MODELER}',
            "date": datetime.date.today().strftime("%Y/%m/%d"),
            "element size": f"{self.station.element_size} m",
            "location": f"on fault, {self.station.location}",
            "minimum time step": f"{self.station.min_time_step:.5e}",
            "maximum time step": f"{self.station.max_time_step:.5e}",
            "num time steps": str(self.station.num_timesteps),
        }

    def write(self, filename="output_data.txt"):
        """Write the metadata and data to a file."""
        with open(filename, "w") as f:
            # Write headers
            for key, value in self.header.items():
                f.write(f"# {key}={value}\n")
            f.write("# The line below lists the names of the data fields\n")  # Blank line separator
            for col in self.column_headers:
                f.write(f"{col} ")
            f.write("\n")     
            f.write("# Here is the time-series data.\n")  # Blank line separator    
            # Write data
            np.savetxt(f, self.station.data, fmt="%.6e", delimiter=" ")

        print(f"Data written to {filename}")

class FaultStation:
    """Stores numerical data for a fault station and provides an interface to generate a FaultFile."""

    def __init__(self, element_size=10, location="0 m, 0 m",
                 min_time_step=1, max_time_step=1, input_data=None, station_number=None):
          """
        Initialize the FaultStation with either synthetic data or input data.

        :param num_timesteps: Number of time steps (default: 100)
        :param element_size: Element size in meters (default: 10)
        :param location: Fault station location as "x2 m, x3 m"
        :param min_time_step: Minimum time step (default: 9.6225e-04)
        :param max_time_step: Maximum time step (default: 5.0636e+04)
        :param input_data: Optional dictionary of input data; if provided, uses real data.
        :param station_number: The ID of the fault station in `input_data` (required if `input_data` is provided).
        """
        self.element_size = element_size
        self.location = location

        if input_data is not None and station_number is not None:
            # If input_data is provided, call process_input_data()
            self.read_input_data(input_data, station_number)
        else:
            # Otherwise, generate synthetic data
            self.generate_data(100)
            self.num_timesteps = num_timesteps
            self.min_time_step = min_time_step
            self.max_time_step = max_time_step     

    def read_input_data(self, input_data, station_number):
        """
        Initialize the FaultStation using simulation input data.

        :param input_data: Dictionary-like dataset containing simulation results.
        :param station_number: ID number of the fault station in the dataset.
        :param location: Fault station location as "x2 m, x3 m".
        """
        # Extract time data
        time = np.asarray(input_data[f'stateVariable Time'])

        # Extract fault station variables
        slip = np.asarray(input_data[f'slip faultStation_{station_number}'])
        slip_velocity = np.asarray(input_data[f'slipVelocity faultStation_{station_number}'])
        shear_traction = np.asarray(input_data[f'shearTraction faultStation_{station_number}'])
        
        # Ensure data has the correct shape
        if slip.shape[1] != 2 or slip_velocity.shape[1] != 2 or shear_traction.shape[1] != 2:
            raise ValueError(f"Data for faultStation_{station_number} is not in expected shape (Nx2).")

        # Assign the components to separate variables
        slip_2, slip_3 = slip[:, 0], slip[:, 1]
        slip_rate_2, slip_rate_3 = np.log10(slip_velocity[:, 0]), np.log10(slip_velocity[:, 1])
        shear_stress_2, shear_stress_3 = shear_traction[:, 0], shear_traction[:, 1]

        # Extract state variable
        state = np.asarray(input_data[f'stateVariable state'])

        self.data = np.column_stack((time, slip_2, slip_3, slip_rate_2, slip_rate_3, shear_stress_2, shear_stress_3, state))
          
    def generate_data(self, num_timesteps):
        """Generate synthetic data with given number of timesteps."""
        time = np.linspace(0, 1, num_timesteps)  # Example time values
        slip_2 = np.sin(2 * np.pi * time)  # Example slip values
        slip_3 = np.cos(2 * np.pi * time)  # Example slip values
        slip_rate_2 = np.log10(np.abs(slip_2) + 1e-10)  # Avoid log10(0)
        slip_rate_3 = np.log10(np.abs(slip_3) + 1e-10)
        shear_stress_2 = np.random.uniform(10, 20, num_timesteps)  # Random shear stress values
        shear_stress_3 = np.random.uniform(15, 25, num_timesteps)
        state = np.log10(time + 1e-10)  # Example state values

        self.data = np.column_stack((time, slip_2, slip_3, slip_rate_2, slip_rate_3, shear_stress_2, shear_stress_3, state))

    def read_data(self, filename):
        """Read data from a file """
        pass

    def get_station_file(self):
        """Return a FaultFile instance linked to this station."""
        return FaultStationFile(self)

def getDataFromHDF5( hdf5FilePath, var_name, set_name):
    # Read HDF5
    data = hdf5_wrapper(f'{hdf5FilePath}').get_copy()
    var = data[f'{var_name} {set_name}']
    var = np.asarray(var)
    time = data[f'{var_name} Time']        

# Example usage
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-dir', type=str, help='input direrctory', default='.')
    parser.add_argument('-o', '--output-dir', type=str, help='output directory', default='.')
    args = parser.parse_args()
    output_dir = os.path.abspath( args.ouput_dir )
    for i in range(13):  # Loop over 13 fault stations
        location_str = f"x2 = {fault_stations_locations[i][0]}, x3 = {fault_stations_locations[i][1]}"
        station = FaultStation(
        num_timesteps=100,  
        element_size=10,  
        location=location_str,  
        min_time_step=9.6225e-04,  
        max_time_step=5.0636e+04  
        )
        file = station.get_station_file()
        filename = output_dir + f"/fault_station_{i}.data"  
        file.write( filename )