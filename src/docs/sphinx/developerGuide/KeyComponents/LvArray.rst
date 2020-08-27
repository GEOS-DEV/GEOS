.. _LvArray:

LvArray
##################################################

LvArray overview
=============================

LvArray is a collection of container classes designed for performance portability in that they are usable on the host and device and provide performance similar to direct pointer manipulation.
It consists of six main template classes:

- ``LvArray::Array``: A multidimensional array.
- ``LvArray::SortedArray``: A sorted unique collection of values stored in a contiguous allocation.
- ``LvArray::ArrayOfArrays``: Provides functionality similar to ``std::vector< std::vector< T > >`` except all the values are stored in one allocation.
- ``LvArray::ArrayOfSets``: Provides functionality similar to ``std::vector< std::set< T > >`` except all the values are stored in one allocation.
- ``LvArray::SparsityPattern``: A compressed row storage sparsity pattern, or equivalently a boolean CRS matrix.
- ``LvArray::CRSMatrix``: A compressed row storage matrix.

These template classes all take three common template arguments:

- The type of the value stored in the container.
- The integral type used for indexing calculations, the recommended type is ``std::ptrdiff_t``.
- The buffer type to use for allocation and de-allocation.


Use in GEOSX
=============================

``LvArray`` containers are used in GEOSX as primary storage mechanism for mesh topology, field data and any other type of "large" data sets (i.e. ones that scale with the size of the problem).
When allocating a new field, using one of ``LvArray`` containers is mandatory if the data is meant to be used in any computational kernel.
The file ``common/DataTypes.hpp`` provides shorthand aliases for commonly used containers:

.. literalinclude:: ../../../../coreComponents/common/DataTypes.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00


LvArray documentation
=====================

Please refer to the full `LvArray documentation <https://lvarray.readthedocs.io/en/latest/>`_ for details on each of the classes.
