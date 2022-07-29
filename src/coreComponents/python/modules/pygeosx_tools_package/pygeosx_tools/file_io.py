import os
import numpy as np


def save_tables(axes, properties, table_root='./tables', axes_names=[]):
    """
    Saves a set of tables in GEOSX format
    
    The shape of these arrays should match the length of each axis in the specified order.
    The output directory will be created if it does not exist yet.
    If axes_names are not supplied, then they will be selected based
    on the dimensionality of the grid: 1D=[t]; 3D=[x, y, z]; 4D=[x, y, z, t].

    Args:
        axes (list): A list of numpy ndarrays defining the table axes
        properties (dict): A dict of numpy ndarrays defning the table values
        table_root (str): The root path for the output directory
        axes_names (list): A list of names for each potential axis (optional)
    """
    # Check to see if the axes, properties have consistent shapes
    axes_size = tuple([len(x) for x in axes])
    axes_dimension = len(axes_size)
    for k, p in properties.items():
        property_size = np.shape(p)
        if (property_size != axes_size):
            print('Property:', k)
            print('Grid size:', axes_size)
            print('Property size', property_size)
            raise Exception('Table dimensions do not match proprerties')

    # Check the axes names
    if axes_names:
        if (axes_dimension != len(axes_names)):
            print('Axes dimensions:', axes_dimension)
            print('Number of axis names provided:', len(axes_names))
            raise Exception('The grid dimensions and axes names do not match')
    else:
        if (axes_dimension == 1):
            axes_names = ['t']
        elif (axes_dimension == 3):
            axes_names = ['x', 'y', 'z']
        elif (axes_dimension == 4):
            axes_names = ['x', 'y', 'z', 't']
        else:
            axes_names = ['x%i' % (ii) for ii in range(axes_dimension)]

    # Write the axes
    os.makedirs(table_root, exist_ok=True)
    for g, a in zip(axes, axes_names):
        np.savetxt('%s/%s.csv' % (table_root, a), g, fmt='%1.5f', delimiter=',')

    for k, p in properties.items():
        np.savetxt('%s/%s.csv' % (table_root, k), np.reshape(p, (-1), order='F'), fmt='%1.5e', delimiter=',')


def load_tables(axes_names, property_names, table_root='./tables', extension='csv'):
    """
    Load a set of tables in GEOSX format

    Args:
        axes_names (list): Axis file names in the target directory (with no extension)
        property_names (list): Property file names in the target directory (with not extension)
        table_root (str): Root path for the table directory
        extension (str): Table file extension (default = 'csv')

    Returns:
        tuple: List of axes values, and dictionary of table values
    """
    # Load axes
    axes = [np.loadtxt('%s/%s.%s' % (table_root, axis, extension), unpack=True) for axis in axes_names]
    N = tuple([len(x) for x in axes])

    # Load properties
    properties = {
        p: np.reshape(np.loadtxt('%s/%s.%s' % (table_root, p, extension)), N, order='F')
        for p in property_names
    }

    return axes, properties
