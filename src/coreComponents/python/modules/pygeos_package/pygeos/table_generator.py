
from numpy import array, loadtxt


def writeGEOSTable(spatial, properties):
  # Open spatial files
  spatial_names = ('x', 'y', 'z', 't')
  for ii in range(0, len(spatial)):
    with open('%s.txt' % spatial_names[ii], 'w') as f:
      for p in spatial[ii]:
        f.write('%s\n' % (p))

  # Open property files
  for k in properties.keys():
    with open('%s.txt' % (k), 'w') as f:
      # Header
      for ii in range(0, len(spatial)):
        f.write('%s\n' % (len(spatial[ii])))

      # Table
      values = array(properties[k]).reshape(-1, order='F')
      for v in values:
        f.write('%1.3e\n' % (v))


def readGEOSTable(sfiles, pfiles):
  # Open spatial files
  spatial = []
  for s in sfiles:
    spatial.append(loadtxt('%s.txt' % (s), unpack=True))

  # Open property files
  properties = {}
  for p in pfiles:
    tmp = loadtxt('%s.txt' % (p), unpack=True)
    properties[p] = tmp[3:].reshape(tmp[:3], order='F')

  return spatial, properties
