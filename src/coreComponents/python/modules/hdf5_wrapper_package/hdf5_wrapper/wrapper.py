import h5py  # type: ignore[import]
import numpy as np
from numpy.core.defchararray import encode, decode
from typing import Union, Dict, Any, Iterable, Optional, Tuple


# Note: I would like to replace Any here with str, float, int, np.ndarray, etc.
#       However, this heterogeneous pattern causes issues with mypy indexing
hdf5_key_types = Union[str, int]
hdf5_get_types = Union['hdf5_wrapper', Any]
nested_dict_type = Dict[hdf5_key_types, Any]
hdf5_set_types = Union['hdf5_wrapper', nested_dict_type, Any]


class hdf5_wrapper():
    """
    A class for reading/writing hdf5 files, which behaves similar to a native dict
    """

    def __init__(self, fname: str = '', target: Optional[h5py.File] = None, mode: str = 'r') -> None:
        """
        Initialize the hdf5_wrapper class

        If the fname is supplied (either by a positional or keyword argument),
        the wrapper will open a hdf5 database from the filesystem.
        The recommended options for the mode flag include 'r' for read-only
        and 'a' for read/write access.
        If write mode is enabled, and the fname does not point
        to an existing file, a new database will be created.

        If the target is supplied, then a new instance of the wrapper will
        be created using an existing database handle.

        Args:
            fname (str): the filename of a new or existing hdf5 database
            target (hdf5_wrapper): the handle of an existing hdf5 dataset
            mode (str): the read/write behavior of the database (default='r')
        """
        self.mode: str = mode
        self.target: h5py.File = target
        if fname:
            self.target = h5py.File(fname, self.mode)

    def __getitem__(self, k: hdf5_key_types) -> hdf5_get_types:
        """
        Get a target from the database

        If the target is not present in the datastructure and the
        database is open in read/write mode, the wrapper will create a
        new group and return an hdf5_wrapper.  Otherwise it will throw an error

        Args:
            k (str): name of target group or array

        Returns:
            hdf5_wrapper/np.ndarray: The returned value
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

    def __setitem__(self, k: hdf5_key_types, value: hdf5_set_types):
        """
        Write an object to the database if write-mode is enabled

        Args:
            k (str): the name of the object
            value (dict, np.ndarray, float, int, str): the object to be written
        """
        if (self.mode in ['w', 'a']):
            if isinstance(value, (dict, hdf5_wrapper)):
                # Recursively add groups and their children
                if (k not in self.target):
                    self.target.create_group(k)
                new_group = self[k]
                for kb, x in value.items():
                    new_group[kb] = x
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

    def link(self, k: hdf5_key_types, target: str) -> None:
        """
        Link an external hdf5 file to this location in the database

        Args:
            k (str): the name of the new link in the database
            target (str): the path to the external database
        """
        self.target[k] = h5py.ExternalLink(target, '/')

    def keys(self) -> Iterable[hdf5_key_types]:
        """
        Get a list of groups and arrays located at the current level

        Returns:
            list: a list of key names pointing to objects at the current level
        """
        if isinstance(self.target, h5py._hl.group.Group):
            return list(self.target)
        else:
            raise ValueError('Object not a group!')

    def values(self) -> Iterable[hdf5_get_types]:
        """
        Get a list of values located on the current level
        """
        return [self[k] for k in self.keys()]

    def items(self) -> Iterable[Tuple[hdf5_key_types, hdf5_get_types]]:
        return zip(self.keys(), self.values())

    def __enter__(self):
        """
        Entry point for an iterator
        """
        return self

    def __exit__(self, type, value, traceback) -> None:
        """
        End point for an iterator
        """
        self.target.close()

    def __del__(self) -> None:
        """
        Closes the database on wrapper deletion
        """
        try:
            if isinstance(self.target, h5py._hl.files.File):
                self.target.close()
        except:
            pass

    def close(self) -> None:
        """
        Closes the database
        """
        if isinstance(self.target, h5py._hl.files.File):
            self.target.close()

    def get_copy(self) -> nested_dict_type:
        """
        Copy the entire database into memory

        Returns:
            dict: a dictionary holding the database contents
        """
        result: Dict[Union[str, int], Any] = {}
        for k in self.keys():
            tmp = self[k]
            if isinstance(tmp, hdf5_wrapper):
                result[k] = tmp.get_copy()
            else:
                result[k] = tmp

        return result

    def copy(self) -> nested_dict_type:
        """
        Copy the entire database into memory

        Returns:
            dict: a dictionary holding the database contents
        """
        return self.get_copy()

    def insert(self, x: Union[nested_dict_type, 'hdf5_wrapper']) -> None:
        """
        Insert the contents of the target object to the current location

        Args:
            x (dict, hdf5_wrapper): the dictionary to insert
        """
        for k, v in x.items():
            self[k] = v
