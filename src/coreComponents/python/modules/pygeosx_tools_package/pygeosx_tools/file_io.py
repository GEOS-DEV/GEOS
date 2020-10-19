
import os
import numpy as np


def save_tables(grid, properties, table_root='./tables', axes=['x', 'y', 'z', 't']):
  """
  @brief saves a set of tables in GEOSX format
  @param grid a list of numpy ndarrays defining the table axes
  @param properties a dict of numpy ndarrays defning the table values
         Note: the shape of these arrays should match the
         length of each axis in the specified order
  @param table_root the root path for the output directory
         Note: this directory will be created if it does not exist yet
  @param axes a list of names for each potential axis (default = ['x', 'y', 'z', 't'])
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
  @brief load a set of tables in GEOSX format
  @param axes a list of axis file names in the target directory (with no extension)
  @param properties a list of property file names in the target directory (with not extension)
  @param table_root the root path for the table directory
  @param extension the table file extension (default = '.csv')
  """
  # Load axes
  X = [np.loadtxt('%s/%s.%s' % (table_root, f, extension), unpack=True) for f in axes]
  N = tuple([len(tmp) for tmp in X])

  # Load properties
  p = {f: np.reshape(np.loadtxt('%s/%s.%s' % (table_root, f, extension)), N, order='F') for f in properties}

  return X, p

