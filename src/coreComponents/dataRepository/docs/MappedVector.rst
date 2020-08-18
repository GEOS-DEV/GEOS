.. _MappedVector:


MappedVector
============

Defines a container that provides an interface which is an superposition of a contiguous array and a hash map.
The types contained in the MappedVector are stored in a ``std::vector``, and thus are accessible through an integral index lookup.
In addition, there is a ``std::unordered_map<KEY,INDEX>`` member that provides a key lookup capability to the container if that is the preferred interface.
Note that the lookup will have some overhead associated with the hash map lookup.

API documentation
-----------------

`MappedVector <../../../doxygen_output/html/classgeosx_1_1_mapped_vector.html>`_

.. todo::
   This page needs to be significantly expanded.
