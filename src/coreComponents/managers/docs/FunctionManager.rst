###############################################################################
Function Manager
###############################################################################

GEOSX functions are specified in the ``Functions`` block in the input xml file.  There are three types of tables: ``TableFunction``, ``SymbolicFunction``, and ``CompositeFunction``.  For example:

.. code-block:: xml

  <Functions>
    <TableFunction name="f_a"
                   inputVarNames="time"
                   coordinates="-10 10"
                   values="-10 10" />

    <SymbolicFunction name="f_b"
                      inputVarNames="ReferencePosition time"
                      variableNames="x y z t"
                      expression="x+y*z"/>

    <CompositeFunction name="composite_ab"
                       inputVarNames="ReferencePosition time"
                       functionNames="f_a f_b"
                       variableNames="a b"
                       expression="a+b"/>
  </Events>





FunctionBase
========================

Each function type derived from FunctionBase requires the following input attribute:

* ``inputVarName`` - Specifies the names of input parameters to the function.  These can be either the name of an array, such as ``Pressure``, or the special keyword, ``time``.


A function can be used in one of three ways:

1. Single mode - Using this option requires the developer to manually construct the input array to the function, of type ``real64 const * const input``.  This method is useful for evaluating functions of single inputs, such as ``time``.

2. Object mode - This method is used to evaluate a function for arrays stored on an object.  When called, it will iterate over the inputVarName list, construct the appropriate array, evaluate the function in sequence, and store the results in an array.  If any of the inputs point to a vector, the individual components will be appended to the input in sequence (e.g.: if ``inputVarName="time, ReferencePosition"``, then ``input[0]=[time, x[0], y[0], z[0]]``).

3. Statistics mode - This method behaves in the same manner as the object mode, except that instead of returning an array with the function results, it will return a structure containing the minmumum, average, and maximum values of the results.



TableFunction
========================

This function type is used to compute an arbitrary-dimensional function based upon table data.  The table axes can be tied to spatial dimensions, or any other set of parameters, such as a phase diagram.  The available xml attributes for this function are:

* ``coordinates`` - The coordinates for the table (1D function only)
* ``values`` - The values for the table (1D function only)
* ``coordinateFiles`` - A list of comma-delimited text files that contain the coordinate definitions for each axis of the table (e.g.: "a.csv, b.csv, c.csv")
* ``voxelFile`` - A comma-delimited file that contains the table values.  For multi-dimensional tables, the values should be given in Fortran order (where the first index changes fastest).
* ``interpolation`` - This defines the interpolation method for the function (only linear interpolation is currently implemented).



SymbolicFunction
=========================

This function type is used to compute a symbolic function from a user-defined string.  It uses a just-in-time compiler to construct the function, so is is nearly as fast as natively compiled C++ code.  The available xml attributes for this function are:

* ``variableNames`` - This is a list of short-hand names for the inputs to the symbolic function.  Each name should map to a single alphabetic character.  There should be a definition for each scalar input and for each component of a vector input.  For example if ``inputVarName="time, ReferencePosition"``, then ``variableNames="t, x, y, z"``
* ``expression`` = This is the definition of the symbolic function.  For the most part, this follows a python-esque format.  However, the function string cannot contain any spaces, and uses the ``pow(x, 3)`` to represent power instead of the ``**`` operator.



CompositeFunction
==============================

This function is derived from the symbolic function.  However, instead of using the time or object as inputs, it is used to combine the outputs of other functions using a symbolic expression.  The available xml attribures are:

* ``functionNames`` - This is a list of the input function names to use
* ``variableNames`` - This is a list of short-hand names for the *functions* in the symbolic expression.    Each name should map to a single alphabetic character.
* ``expression`` = This is the definition of the symbolic function.  For the most part, this follows a python-esque format.  However, the function string cannot contain any spaces, and uses the ``pow(x, 3)`` to represent power instead of the ``**`` operator.

For example, the composite function in this example, would compute the expression ``f = x^2 + t^3``:

.. code-block:: xml

  <Functions>
    <SymbolicFunction name="f_a"
                      inputVarNames="ReferencePosition"
                      variableNames="x y z"
                      expression="pow(x,2)"/>

    <SymbolicFunction name="f_b"
                      inputVarNames="time"
                      variableNames="t"
                      expression="pow(t,3)"/>

    <CompositeFunction name="composite_ab"
                       inputVarNames="ReferencePosition time"
                       functionNames="f_a f_b"
                       variableNames="a b"
                       expression="a+b"/>
  </Events>




