import xml.etree.ElementTree as et
import os
import numpy
import vtk
from vtk.util import numpy_support
import h5py
import argparse

def read_input_grid_from_vti(vti_file_path: str) -> vtk.vtkImageData:
    """
    Reads a vtkImageData grid from a .vti file.

    Args:
        vti_file_path (str): The path to the .vti file.

    Returns:
        vtk.vtkImageData: The loaded image data grid.
    """
    if not os.path.exists(vti_file_path):
        raise FileNotFoundError(f"VTI file not found: {vti_file_path}")

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(vti_file_path)
    reader.Update()
    
    grid = reader.GetOutput()
    if grid is None or grid.GetNumberOfPoints() == 0:
        raise ValueError(f"Failed to read valid grid data from {vti_file_path}")
        
    return grid

def read_grid_data(time_step: float, input_grid: vtk.vtkImageData, hdf5_file_name: str, pvd_file_path: str):
    """
    Reads VTK output files for a specific time step, extracts cell data,
    maps it onto a 2D input grid based on element centers, and writes to an HDF5 file.
    
    Args:
        time_step (float): The time step to extract.
        input_grid (vtk.vtkImageData): The 2D grid to map the data onto.
        hdf5_file_name (str): Path to the output HDF5 file.
        pvd_file_path (str): Path to the main .pvd file.
    """
    base_dir = os.path.dirname(pvd_file_path)
    
    tree = et.parse(pvd_file_path)
    root = tree.getroot()
    
    vtm_file_rel = None
    # Assuming standard PVD format where root has a Collection child
    collection = root.find('Collection')
    if collection is not None:
        for dataset in collection.findall('DataSet'):
            ts = float(dataset.get('timestep'))
            if numpy.isclose(ts, time_step):
                vtm_file_rel = dataset.get('file')
                break
    
    if not vtm_file_rel:
        raise ValueError(f"Time step {time_step} not found in PVD file.")
        
    vtm_file_path = os.path.join(base_dir, vtm_file_rel)
    
    vtm_reader = vtk.vtkXMLMultiBlockDataReader()
    vtm_reader.SetFileName(vtm_file_path)
    vtm_reader.Update()
    
    multiblock = vtm_reader.GetOutput()
    if multiblock is None:
        raise ValueError(f"Failed to read valid multiblock data from {vtm_file_path}")

    dims = input_grid.GetDimensions()
    # GetDimensions returns the number of points. We need the number of cells for cell data.
    # For a 2D grid, dimensions will be (Nx, 1, Nz) or (Nx, Ny, 1). Assuming X-Z plane with Ny=1.
    nx_cells = max(1, dims[0] - 1)
    ny_cells = max(1, dims[1] - 1)
    nz_cells = max(1, dims[2] - 1)
    
    if ny_cells != 1 and dims[1] != 1:
         print("Warning: Expected ny=1 for this specific 2D grid orientation (X-Z plane).")
         
    bounds = input_grid.GetBounds()
    xmin, xmax = bounds[0], bounds[1]
    # ymin, ymax = bounds[2], bounds[3] # Y is trivial
    zmin, zmax = bounds[4], bounds[5]
    
    spacing = input_grid.GetSpacing()
    dx, dz = spacing[0], spacing[2]

    # Dictionary to hold the aggregated arrays. Key: Property Name, Value: 2D numpy array
    data_arrays = {}

    input_cell_data = input_grid.GetCellData()
    
    # Read only explicitly the attribute named "attribute" from the input grid
    array_vtk = input_cell_data.GetArray("attribute")
    if array_vtk is not None:
        data = numpy_support.vtk_to_numpy(array_vtk)
        # As it is a fully defined integer property with a single component, directly reshape
        data_arrays["attribute"] = data.reshape((nz_cells, nx_cells))
    else:
        print("Warning: Array 'attribute' not found in input grid.")

    datasets = []
    iterator = multiblock.NewIterator()
    iterator.InitTraversal()
    while not iterator.IsDoneWithTraversal():
        obj = iterator.GetCurrentDataObject()
        if obj is not None and obj.IsA("vtkUnstructuredGrid"):
            datasets.append(obj)
        iterator.GoToNextItem()

    for dataset in datasets:
        cell_data = dataset.GetCellData()
        num_arrays = cell_data.GetNumberOfArrays()
        
        if cell_data.HasArray("elementCenter"):
            element_centers = numpy_support.vtk_to_numpy(cell_data.GetArray("elementCenter"))
        else:
            print("Warning: 'elementCenter' not found in block. Skipping.")
            continue
            
        x_centers = element_centers[:, 0]
        z_centers = element_centers[:, 2]
        
        i_indices = numpy.floor((x_centers - xmin) / dx).astype(int)
        k_indices = numpy.floor((z_centers - zmin) / dz).astype(int)
        
        # Filter out points that fall outside the input grid (optional, for safety)
        valid_mask = (i_indices >= 0) & (i_indices < nx_cells) & (k_indices >= 0) & (k_indices < nz_cells)
        
        if not numpy.any(valid_mask):
            continue
            
        valid_i = i_indices[valid_mask]
        valid_k = k_indices[valid_mask]
        
        for i in range(num_arrays):
            array_vtk = cell_data.GetArray(i)
            original_name = array_vtk.GetName()
            
            if original_name == "elementCenter":
                continue
                
            # Logic to reconcile property names across blocks
            # e.g., SAND_phaseCapPressure -> phaseCapPressure
            prop_name = original_name
            if "_" in original_name:
                parts = original_name.split("_", 1) # Split only on the first underscore
                if len(parts) == 2 and len(parts[1]) >= 5:
                    prop_name = parts[1]
                    
            data = numpy_support.vtk_to_numpy(array_vtk)
            
            # Initialize the 2D array with appropriate fill value and dtype if this is the first time we see this property
            if prop_name not in data_arrays:
                np_dtype = data.dtype
                fill_value = numpy.nan if numpy.issubdtype(np_dtype, numpy.floating) else 0
                num_components = array_vtk.GetNumberOfComponents()
                
                if num_components == 1:
                    data_arrays[prop_name] = numpy.full((nz_cells, nx_cells), fill_value, dtype=np_dtype)
                else:
                     data_arrays[prop_name] = numpy.full((nz_cells, nx_cells, num_components), fill_value, dtype=np_dtype)
                     
            # Map the data onto the 2D grid
            # Note: Numpy arrays are typically indexed [row, col] -> [z, x]
            if array_vtk.GetNumberOfComponents() == 1:
                data_arrays[prop_name][valid_k, valid_i] = data[valid_mask]
            else:
                data_arrays[prop_name][valid_k, valid_i, :] = data[valid_mask]

    SAVE_PROPS = """
attribute
phaseRelPerm
phaseCompFraction
phaseViscosity
phaseMassDensity
phaseCapPressure
phaseFraction
phaseVolumeFraction
globalCompDensity
pressure
globalCompFraction
""".split()
    with h5py.File(hdf5_file_name, 'w') as h5f:
        # Save grid metadata (optional but good practice)
        h5f.attrs['time_step'] = time_step
        h5f.attrs['grid_dimensions'] = [nx_cells, nz_cells]
        h5f.attrs['grid_spacing'] = [dx, dz]
        h5f.attrs['grid_origin'] = [xmin, zmin]
        
        for prop_name, array_2d in data_arrays.items():
            if prop_name not in SAVE_PROPS: continue
            # HDF5 handles NaNs naturally for float datasets
            h5f.create_dataset(prop_name, data=array_2d, compression="gzip")
            
    print(f"Successfully wrote data for time step {time_step} to {hdf5_file_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert VTK/PVD output to HDF5 mapped on a 2D grid.")
    parser.add_argument("--input", type=str, required=True, help="Path to the input .vti grid file.")
    parser.add_argument("--output", type=str, required=True, help="Path to the output .h5 file.")
    parser.add_argument("--pvd", type=str, required=True, help="Path to the main .pvd file.")
    parser.add_argument("--time", type=float, required=True, help="Time step to extract from the PVD file.")
    
    args = parser.parse_args()

    try:
        # Read the input grid from a VTI file
        input_grid = read_input_grid_from_vti(args.input)

        # Convert data for a specific time step
        read_grid_data(args.time, input_grid, args.output, args.pvd)
    except Exception as e:
        print(f"Error: {e}")
