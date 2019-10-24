
import numpy as np


# Write an GEOS-compatible ascii table
#   axes_values = A list of the axes values in order
#   properties = A dictionary of properties given as an ndarray
#                The shape of the array must be consistent with the axes
#   axes_names = A list of axes names (optional)
#   string_format = Format for table values
def write_GEOS_table(axes_values, properties, axes_names=['x', 'y', 'z', 't'], string_format='%1.5e'):
  # Check to make sure the axes/property files have the correct shape
  axes_shape = tuple([len(x) for x in axes_values])
  for k in properties.keys():
    if (np.shape(properties[k]) != axes_shape):
      raise Exception("Shape of parameter %s is incompatible with given axes" % (k))

  # Write axes files
  for ii in range(0, len(axes_values)):
    np.savetxt('%s.geos' % (axes_names[ii]), axes_values[ii], fmt=string_format, delimiter=',')

  # Write property files
  for k in properties.keys():
    tmp = np.reshape(properties[k], (-1), order='F')
    np.savetxt('%s.geos' % (k), tmp, fmt=string_format, delimiter=',')


# Read an GEOS-compatible ascii table
#   axes_files = A list of the axes file names in order
#   property_files = A list of the property file names
def read_GEOS_table(axes_files, property_files):
  # Open spatial files
  axes_values = []
  for f in axes_files:
    axes_values.append(np.loadtxt('%s.geos' % (f), unpack=True, delimiter=','))
  axes_shape = tuple([len(x) for x in axes_values])

  # Open property files
  properties = {}
  for f in property_files:
    tmp = np.loadtxt('%s.geos' % (f), unpack=True, delimiter=',')
    properties[f] = np.reshape(tmp, axes_shape, order='F')

  return axes_values, properties


# Example of how to read/write GEOS tables using the above functions
def write_read_GEOS_table_example():
  # Define table axes
  a = np.array([0.0, 1.0])
  b = np.array([0.0, 0.5, 1.0])
  axes_values = [a, b]

  # Generate table values (note: the indexing argument is important)
  A, B = np.meshgrid(a, b, indexing='ij')
  properties = {'c': A + 2.0*B}

  # Write, then read tables
  write_GEOS_table(axes_values, properties, axes_names=['a', 'b'])
  axes_b, properties_b = read_GEOS_table(['a', 'b'], ['c'])



