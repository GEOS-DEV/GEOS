import h5py
import numpy as np
from numpy.core.defchararray import encode, decode


class hdf5_wrapper():
    """
    @brief a class for reading/writing hdf5 files, which behaves similar to a native dict
    """

    def __init__(self, fname='', target='', mode='r'):
        """
        @brief initialize the hdf5_wrapper class
        @param fname the filename of a new or existing hdf5 database
        @param target the handle of an existing hdf5 dataset
        @param mode the read/write behavior of the database (default='r')

        @details If the fname is supplied (either by a positional or keyword argument),
                 the wrapper will open a hdf5 database from the filesystem.  The reccomended
                 options for the mode flag include 'r' for read-only and 'a' for read/write
                 access.  If write mode is enabled, and the fname does not point
                 to an existing file, a new database will be created.

                 If the target is supplied, then a new instance of the wrapper will
                 be created using an existing database handle.
        """

        self.mode = mode
        self.target = target
        if fname:
            self.target = h5py.File(fname, self.mode)

    def __getitem__(self, k):
        """
        @brief get a target from the database
        @param k name of target group or array

        @return the returned value depends on the type of the target:
                  - An existing hdf5 group will return an instance of hdf5_wrapper
                  - An existing array will return an numpy ndarray
                  - If the target is not present in the datastructure and the
                    database is open in read/write mode, the wrapper will create a
                    new group and return an hdf5_wrapper
                  - Otherwise, this will throw an error
        """
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
        """
        @brief write an object to the database if write-mode is enabled
        @param k the name of the object
        @param value the object to be written
        """

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
                    del (self.target[k])

                # Add everything else as an ndarray
                tmp = np.array(value)
                if (tmp.dtype.kind in ['S', 'U', 'O']):
                    tmp = encode(tmp)
                self.target[k] = tmp
        else:
            raise ValueError(
                'Cannot write to an hdf5 opened in read-only mode!  This can be changed by overriding the default mode argument for the wrapper.'
            )

    def link(self, k, target):
        """
        @brief link an external hdf5 file to this location in the database
        @param k the name of the new link in the database
        @param target the path to the external database
        """
        self.target[k] = h5py.ExternalLink(target, '/')

    def keys(self):
        """
        @brief get a list of groups and arrays located at the current level
        @return a list of strings 
        """
        if isinstance(self.target, h5py._hl.group.Group):
            return list(self.target)
        else:
            raise ValueError('Object not a group!')

    def __enter__(self):
        """
        @brief entry point for an iterator
        """
        return self

    def __exit__(self, type, value, traceback):
        """
        @brief end point for an iterator
        """
        self.target.close()

    def __del__(self):
        """
        @brief closes the database on wrapper deletion
        """
        try:
            if isinstance(self.target, h5py._hl.files.File):
                self.target.close()
        except:
            pass

    def close(self):
        """
        @brief closes the database
        """
        if isinstance(self.target, h5py._hl.files.File):
            self.target.close()

    def get_copy(self):
        """
        @brief copy the entire database into memory
        @return a dictionary holding the database contents
        """
        tmp = {}
        self.copy(tmp)
        return tmp

    def copy(self, output):
        """
        @brief pack the contents of the current database level onto the target dict
        @param output the dictionary to pack objects into
        """
        for k in self.keys():
            tmp = self[k]

            if isinstance(tmp, hdf5_wrapper):
                output[k] = {}
                tmp.copy(output[k])
            else:
                output[k] = tmp
