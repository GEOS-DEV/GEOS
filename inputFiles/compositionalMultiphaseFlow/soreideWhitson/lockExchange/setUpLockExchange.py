import vtk
import numpy as np
import pandas as pd
import csv
import os
#from generateMesh import *



#### Grid size
nx = 50
ny = 1
nz = nx 

#### Problem Setup
# Physical size
Lx = 50
Ly = 1
Lz = 50


# Initial condition
# CH4, CO2, H2S, H2O
zRIGHT = [0.0, 0.0, 0.0, 1.0]
zLEFT = [0.8, 0.1, 0.1, 0.0]

#### Create the folder for the resolution
folder_path = "InitCond"
os.makedirs(folder_path, exist_ok=True)

#### Generate mesh ####
#generateRectangularMesh(nx, ny, nz, Lx, Ly, Lz, folder_path)

#### Produce x.csv, y.csv and z.csv
data_x = np.linspace(0, Lx, nx, endpoint=False) + (Lx / nx) / 2
data_y = np.linspace(0, Ly, ny, endpoint=False) + (Ly / ny) / 2
data_z = np.linspace(-2100, -2100+Lz, nz, endpoint=False) + (Lz / nz) / 2

with open(os.path.join(folder_path, 'x.csv'), 'w', newline='') as file_x:
    writer_x = csv.writer(file_x)
    for value in data_x:
        writer_x.writerow([value])

with open(os.path.join(folder_path, 'y.csv'), 'w', newline='') as file_y:
    writer_y = csv.writer(file_y)
    for value in data_y:
        writer_y.writerow([value])

with open(os.path.join(folder_path, 'z.csv'), 'w', newline='') as file_z:
    writer_z = csv.writer(file_z)
    for value in data_z:
        writer_z.writerow([value])

#### Produce initial z for all components

# Initialize the components data dictionary
components_data = {
    'zCH4': [],
    'zCO2': [],
    'zH2S': [],
    'zH2O': []
}

# Generate the data for each component
for _ in range(nz):
    for i in range(nx):
        if i < nx // 2:
            components_data['zCH4'].append(zLEFT[0])
            components_data['zCO2'].append(zLEFT[1])
            components_data['zH2S'].append(zLEFT[2])
            components_data['zH2O'].append(zLEFT[3])
        else:
            components_data['zCH4'].append(zRIGHT[0])
            components_data['zCO2'].append(zRIGHT[1])
            components_data['zH2S'].append(zRIGHT[2])
            components_data['zH2O'].append(zRIGHT[3])

# Save each component's data to a CSV file
for component, data in components_data.items():
    file_name = f"{component}.csv"
    with open(os.path.join(folder_path, file_name), 'w', newline='') as file:
        writer = csv.writer(file)
        for value in data:
            writer.writerow([value])




