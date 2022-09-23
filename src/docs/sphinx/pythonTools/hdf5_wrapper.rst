
HDF5 Wrapper
--------------------------

The `hdf5_wrapper` python package adds a wrapper to `h5py` that greatly simplifies reading/writing to/from hdf5-format files.


Usage
^^^^^^^

Once loaded, the contents of a file can be navigated in the same way as a native python dictionary.

.. code-block:: python

  import hdf5_wrapper

  data = hdf5_wrapper.hdf5_wrapper('data.hdf5')

  test = data['test']
  for k, v in data.items():
    print('key: %s, value: %s' % (k, str(v)))


If the user indicates that a file should be opened in write-mode (`w`) or read/write-mode (`a`), then the file can be created or modified.
Note: for these changes to be written to the disk, the wrapper may need to be closed or deleted.

.. code-block:: python

  import hdf5_wrapper
  import numpy as np

  data = hdf5_wrapper.hdf5_wrapper('data.hdf5', mode='w')
  data['string'] = 'string'
  data['integer'] = 123
  data['array'] = np.random.randn(3, 4, 5)
  data['child'] = {'float': 1.234}


Existing dictionaries can be placed on the current level:

.. code-block:: python

  existing_dict = {'some': 'value'}
  data.insert(existing_dict)


And external hdf5 format files can be linked together:

.. code-block:: python

  for k in ['child_a', 'child_b']:
    data.link(k, '%s.hdf5' % (k))



API
^^^^^

.. automodule:: hdf5_wrapper.wrapper
    :members:
