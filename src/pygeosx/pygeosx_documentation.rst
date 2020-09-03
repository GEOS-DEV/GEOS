GEOSX in Python
=================
GEOSX can be manipulated and executed through a Python script.

High-level control of GEOSX is managed through the top-level ``pygeosx`` functions,
like ``initialize`` and ``run``. GEOSX's data can be manipulated by getting
:ref:`pylvarray <pylvarray>` views of LvArray objects living in GEOSX's data repository.
These ``pylvarray`` views are fetched by calling ``Wrapper.value()`` after getting
a ``Wrapper`` out of the data repository.

.. warning:: 
	The ``pygeosx`` module provides plenty of opportunites to crash Python. 
	See the Segmentation Faults section below.

Module Classes
--------------

.. automodule:: pygeosx
   :members:


Segmentation Faults
-------------------
Improper use of this module and associated programs can easily cause Python to crash. 
There are two main causes of crashes. The first is under user control, the second is not.

Stale Numpy Views
^^^^^^^^^^^^^^^^^
The ``pylvarray`` classes (which may be returned from ``Wrapper.value()``) 
provide various ways to get Numpy views of
their data. However, those views are only valid as long as the 
LvArray object's buffer is not reallocated. The buffer may be reallocated
by invoking methods (the ones that require 
the ``pylarray.RESIZEABLE`` permission) or by calls into ``pygeosx``. 
It is strongly recommended that you do not keep Numpy views of 
LvArray objects around after calls to ``pygeosx``.

.. code:: python

	my_array = pygeosx.get_wrapper("path").value()
	view = my_array.to_numpy()
	my_array.resize(1000)
	print(view) # segfault

Destroyed LvArray C++ objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As mentioned earlier, the classes defined in this
module cannot be created in Python; ``pygeosx`` must create an LvArray 
object in C++, then create a
``pylvarray`` view of it. However, the Python view will only be
valid as long as the underlying LvArray C++ object is kept around. If 
that is destroyed, the Python object will be left holding an invalid
pointer and subsequent attempts to use the Python object will cause 
undefined behavior. Unfortunately, ``pygeosx`` may destroy LvArray
objects without warning, and there is nothing users can do.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

