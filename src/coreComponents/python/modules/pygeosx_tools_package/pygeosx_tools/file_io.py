
import os
import numpy as np


def save_tables(grid, properties, table_root='./tables', axes=['x', 'y', 'z']):
  os.makedirs(table_root, exist_ok=True)
  for ii, d in zip(range(len(axes)), axes):
    np.savetxt('%s/%s.csv' % (table_root, d),
               grid[ii],
               fmt='%1.5f',
               delimiter=',')

  for ka in properties:
    np.savetxt('%s/%s.csv' % (table_root, ka),
               np.reshape(properties[ka], (-1), order='F'),
               fmt='%1.5e',
               delimiter=',')

