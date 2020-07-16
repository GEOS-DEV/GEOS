
import h5py
import numpy as np
from numpy.core.defchararray import encode, decode


class hdf5_wrapper():
  def __init__(self, fname='', target='', mode='r'):
    self.mode = mode
    self.target = target
    if fname:
      self.target = h5py.File(fname, self.mode)

  def __getitem__(self, k):
    # Test to see if the value is in the database, create it if opened in read/write mode
    if (k not in self.target):
      if (self.mode in ['w', 'a']):
        self.target.create_group(k)
      else:
        raise ValueError('Entry does not exist in database: %s' % (k))

    tmp = self.target[k]

    if isinstance(tmp, h5py._hl.group.Group):
      return hdf5_wrapper(target=tmp, mode=self.mode)
    elif isinstance(tmp, h5py._hl.dataset.Dataset):
      tmp = np.array(tmp)

      # Decode any string types
      if (tmp.dtype.kind in ['S', 'U', 'O']):
        tmp = decode(tmp)

      # Convert any 0-length arrays to native types
      if not tmp.shape:
        tmp = tmp[()]

      return tmp
    else:
      return tmp

  def __setitem__(self, k, value):
    if (self.mode in ['w', 'a']):
      if isinstance(value, dict):
        # Recursively add groups and their children
        if (k not in self.target):
          self.target.create_group(k)
        new_group = self[k]
        for x in value:
          new_group[x] = value[x]
      else:
        # Delete the old copy if necessary
        if (k in self.target):
          del(self.target[k])

        # Add everything else as an ndarray
        tmp = np.array(value)
        if (tmp.dtype.kind in ['S', 'U', 'O']):
          tmp = encode(tmp)
        self.target[k] = tmp
    else:
      raise ValueError('Cannot write to an hdf5 opened in read-only mode!  This can be changed by overriding the default mode argument for the wrapper.')

  def link(self, k, target):
    self.target[k] = h5py.ExternalLink(target, '/')

  def keys(self):
    if isinstance(self.target, h5py._hl.group.Group):
      return list(self.target)
    else:
      raise ValueError('Object not a group!')

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.target.close()

  def __del__(self):
    try:
      if isinstance(self.target, h5py._hl.files.File):
        self.target.close()
    except:
      pass

  def close(self):
    if isinstance(self.target, h5py._hl.files.File):
      self.target.close()

  def get_copy(self):
    tmp = {}
    self.copy(tmp)
    return tmp

  def copy(self, output):
    for k in self.keys():
      tmp = self[k]

      if isinstance(tmp, hdf5_wrapper):
        output[k] = {}
        tmp.copy(output[k])
      else:
        output[k] = tmp


