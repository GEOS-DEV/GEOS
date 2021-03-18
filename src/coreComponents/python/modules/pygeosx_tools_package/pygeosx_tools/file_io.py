
import os
import numpy as np


def save_tables(grid, properties, table_root='./tables', axes=['x', 'y', 'z', 't']):
    """
    Saves a set of tables in GEOSX format
    Notes: The shape of these arrays should match the length of each axis in the specified order
           The output directory will be created if it does not exist yet

    Args:
        grid (list): A list of numpy ndarrays defining the table axes
        properties (dict): A dict of numpy ndarrays defning the table values
        table_root (str): The root path for the output directory
        axes (list): A list of names for each potential axis (default = ['x', 'y', 'z', 't'])
    """
    os.makedirs(table_root, exist_ok=True)
    for ii, d in zip(range(len(grid)), axes):
        np.savetxt('%s/%s.csv' % (table_root, d),
                   grid[ii],
                   fmt='%1.5f',
                   delimiter=',')

    for ka in properties:
        np.savetxt('%s/%s.csv' % (table_root, ka),
                   np.reshape(properties[ka], (-1), order='F'),
                   fmt='%1.5e',
                   delimiter=',')


def load_tables(axes, properties, table_root='./tables', extension='csv'):
    """
    Load a set of tables in GEOSX format

    Args:
        axes (list): Axis file names in the target directory (with no extension)
        properties (list): Property file names in the target directory (with not extension)
        table_root (str): Root path for the table directory
        extension (str): Table file extension (default = 'csv')

    Returns:
        tuple: List of axes values, and dictionary of table values
    """
    # Load axes
    X = [np.loadtxt('%s/%s.%s' % (table_root, f, extension), unpack=True) for f in axes]
    N = tuple([len(tmp) for tmp in X])

    # Load properties
    p = {f: np.reshape(np.loadtxt('%s/%s.%s' % (table_root, f, extension)), N, order='F') for f in properties}

    return X, p

