import vtk
import csv
import os
import xml.etree.ElementTree as ET
import pandas as pd
import matplotlib.pyplot as plt

def read_vtk_file(vtk_filename):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtk_filename)
    reader.Update()
    return reader.GetOutput()

def collect_headers(headers):
    # Collect headers for specified fields
    headers.extend(['globalCompFraction_0', 'globalCompFraction_1', 'globalCompFraction_2', 
                    'phaseVolumeFraction_0', 'phaseVolumeFraction_1', 'pressure'])

def create_table(vtk_data, headers):
    table = []
    num_cells = vtk_data.GetNumberOfCells()
    cell_data = vtk_data.GetCellData()
    
    for cell_id in range(num_cells):
        row = [cell_id]
        for header in headers:
            array_name, component = header.rsplit('_', 1) if '_' in header else (header, None)
            array = cell_data.GetArray(array_name)
            if array:
                if component is not None:
                    row.append(array.GetComponent(cell_id, int(component)))
                else:
                    row.append(array.GetTuple(cell_id)[0])
            else:
                row.append(None)
        table.append(row)
    
    return table

def write_to_csv(table, csv_filename, headers):
    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # Write headers
        writer.writerow(['index'] + headers)
        
        # Write rows from the table
        for row in table:
            writer.writerow(row)
    
    # Log the saved CSV file
    print(f"Saved CSV file: {csv_filename}")

def process_vtm_file(vtm_filename, headers):
    tree = ET.parse(vtm_filename)
    root = tree.getroot()
    
    for dataset in root.findall('.//DataSet'):
        vtu_file = dataset.get('file')
        vtu_path = os.path.join(os.path.dirname(vtm_filename), vtu_file)
        vtk_data = read_vtk_file(vtu_path)
        
        # Create table and return it
        table = create_table(vtk_data, headers)
        return table

def process_pvd_file(pvd_filename, output_dir, Nsnap):
    headers = []
    collect_headers(headers)
    
    tree = ET.parse(pvd_filename)
    root = tree.getroot()
    
    tables_to_plot = []
    time_snapshots = []
    
    total_snapshots = len(root.findall('.//DataSet'))
    snapshots_to_plot = [i * total_snapshots // Nsnap for i in range(1, Nsnap + 1)]

    
    time_snapshot = 0
    for dataset in root.findall('.//DataSet'):
        if time_snapshot in snapshots_to_plot:
            vtm_file = dataset.get('file')
            vtm_path = os.path.join(os.path.dirname(pvd_filename), vtm_file)
            
            timestep_in_seconds = float(dataset.get('timestep'))
            time_in_days = timestep_in_seconds / 86400
            
            table = process_vtm_file(vtm_path, headers)
            
            # Create a unique CSV filename for each time snapshot
            csv_filename = os.path.join(output_dir, f'time_snapshot_{time_snapshot}.csv')
            
            # Write to CSV
            write_to_csv(table, csv_filename, headers)
            
            tables_to_plot.append((table, time_in_days))
        
        time_snapshot += 1
    
    print_headers(headers)
    
    return tables_to_plot

def print_headers(headers):
    print("Headers:")
    for header in headers:
        print(f"  {header}")

def plot_fields(tables_to_plot, output_dir):
    for i, (table, time_in_days) in enumerate(tables_to_plot):
        data = list(zip(*table))
        x = range(100)

        comps = ['CH4', 'CO2', 'H2O']
        sats = ['Soil', 'Sgas']

        # Plot globalCompFraction fields in the same figure
        plt.figure()
        for j in range(3):
            plt.plot(x, data[j+1], label=comps[j])
        
        plt.xlabel('Index')
        plt.ylabel('globalCompFraction')
        plt.ylim(0, 1)
        plt.xlim(0, 99)
        plt.legend()
        plt.title(f'globalCompFraction vs Index (Time: {time_in_days:.2f} days)')
        plt.savefig(os.path.join(output_dir, f'globalCompFraction_{i}.png'))
        plt.close()
        
        # Plot phaseVolumeFraction fields in the same figure
        plt.figure()
        for j in range(3, 5):
            plt.plot(x, data[j+1], label=sats[j-3])
        
        plt.xlabel('Index')
        plt.ylabel('phaseVolumeFraction')
        plt.ylim(0, 1)
        plt.xlim(0, 99)
        plt.legend()
        plt.title(f'phaseVolumeFraction vs Index (Time: {time_in_days:.2f} days)')
        plt.savefig(os.path.join(output_dir, f'phaseVolumeFraction_{i}.png'))
        plt.close()
        
        # Convert pressure from Pa to bars and plot it
        pressure_data = [p / 1e5 for p in data[6]]
        
        plt.figure()
        plt.plot(x, pressure_data)
        plt.xlabel('Index')
        plt.ylabel('Pressure (bars)')
        plt.xlim(0, 99)
        
        # Add time to title (convert from seconds to days by dividing by 86400)
        plt.title(f'Pressure vs Index (Time: {time_in_days:.2f} days)')
        
        plt.savefig(os.path.join(output_dir, f'pressure_{i}.png'))
        plt.close()

def plot_rates(df):
    # Plot Phase0 and Phase1 rates as a function of Time [s]
    plt.figure(figsize=(10, 6))
    plt.plot(df['Time [s]']/86400, df['Phase0 reservoir volumetric rate [rm3/s]']*(-1), label='water rate')
    plt.plot(df['Time [s]']/86400, df['Phase1 reservoir volumetric rate [rm3/s]']*(-1), label='gas rate')

    #plt.plot(df['Time [s]']/86400, df['Phase0 surface volumetric rate [sm3/s]']*(-1), label='water rate')
    #plt.plot(df['Time [s]']/86400, df['Phase1 surface volumetric rate [sm3/s]']*(-1), label='gas rate')

    # Add labels and title
    plt.xlabel('Time [days]')
    plt.ylabel('Production Rate [m3/s]')
    #plt.title('Phase0 and Phase1 Rates as a Function of Time')
    plt.legend()
    #plt.yscale('log')


    # Save the plot as production_rates.png
    plt.savefig('production_rates.png')

def plot_bhp(df1, df2):
    # Plot Phase0 and Phase1 rates as a function of Time [s]
    plt.figure(figsize=(10, 6))
    plt.plot(df1['Time [s]']/86400, df1['BHP [Pa]']/1e5, label='Gas Injector BHP')
    plt.plot(df2['Time [s]']/86400, df2['BHP [Pa]']/1e5, label='Water Injector BHP')

    # Add labels and title
    plt.xlabel('Time [days]')
    plt.ylabel('BHP [bar]')
    #plt.title('Phase0 and Phase1 Rates as a Function of Time')
    plt.legend()

    # Save the plot as production_rates.png
    plt.savefig('injection_bhp.png')

# Print a success message
print("The plot has been saved as production_rates.png")


if __name__ == "__main__":
    pvd_filename = '1D_HCGI/VTK.OUTPUT.pvd'  # Replace with your PVD file path
    output_dir = '1D_HCGI_csv'       # Replace with your desired output directory
    Nsnap = 200 # number of snapshots
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    tables_to_plot = process_pvd_file(pvd_filename, output_dir, Nsnap)
    
    plot_fields(tables_to_plot, output_dir)

    ## Read rates
    df = pd.read_csv('1D_HCGI/WELL.SOLVER_rates/PRODUCER.CONTROL.csv')
    ## Plot rates
    plot_rates(df)
    ## Read bhp
    df1 = pd.read_csv('1D_HCGI/WELL.SOLVER_rates/GAS.INJECTOR.CONTROL.csv')
    df2 = pd.read_csv('1D_HCGI/WELL.SOLVER_rates/WATER.INJECTOR.CONTROL.csv')
    ## Plot rates
    plot_bhp(df1, df2)