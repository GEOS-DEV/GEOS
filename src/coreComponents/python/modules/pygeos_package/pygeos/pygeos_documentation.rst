
.. _advanced_xml_features:

###############################################################################
Advanced XML Features (pygeos)
###############################################################################

Pygeos is a python module that enables advanced .xml features in GEOSX.
It allows users to include child files, parameters, units, and symbolic math. 


Setup
=================================
The source of the pygeos module is located in the GEOSX repository: `src/coreComponents/python/modules/pygeos_package`.
It can be installed within the GEOSX build directory via the `make pygeos` command.
The pygeos executable will be located in the bin directory.

Notes: The above command creates a virtual python environment derived from the version of python used with cmake.
To build the environment, the `virtualenv` package must be installed.
The environment can also be constructed by calling `scripts/setupVirtualPythonEnvironment.bash`.



Usage
=================================

Pygeos can be used either as a command-line tool or be imported into a python script.  Note: pygeos will be on your path if you are in a virtual environment.  Otherwise, it can be called directly from the bin directory of your python distribution.



Command-line Example
------------------------------

Pygeos can be called from the command line to process .xml files via the script.
The following will read a raw .xml file, generate a processed version, and return the new file name:

``pygeos input_file.xml``

Optional arguments include:

- ``-o/--output``: The desired name for the output file (otherwise, it is randomly generated)
- ``-s/--schema``: The location of a schema to validate the final .xml file
- ``-v/--verbose``: Increase module verbosity

For convenience, this script can be embedded within GEOSX arguments:

``srun -n 16 geosx -i `pygeos input_file.xml```


Script-based Example
------------------------------

The pygeos module can also be called from within a python script.  For example:

.. code-block:: python

  import sys
  import pygeos

  new_filename = pygeos.preprocessGEOSXML(sys.argv[1])
  print(new_filename)



Advanced XML Features
=================================

The xml preprocessor in pygeos is designed to take a raw input file, and generate a new file that can be directly read by GEOSX.
The syntax for the advanced XML format is given below.
During the processing the order of operations are:

1) Merging any included XML files into the root structure
2) Substituting in any parameters
3) Evaluating unit strings
4) Evaluating symbolic math
5) Error checking and validation


Including Child XML Files
------------------------------
XML inputs can point to included children (these children can then include grandchildren and so on).
During processing, these are recursively inserted into the root XML structure by the following scheme:

- Merge two objects if:

    - At the root level an object with the matching tag exists.
    - If the "name" attribute is present and an object with the matching tag and name exists.
    - Any preexisting attributes are overwritten by the donor.
- Otherwise append the XML structure with the target.


.. code-block:: xml

  <Included>
    <File name='/path/to/included_a.xml'/>
    <File name='/path/to/included_b.xml'/>
  </Included>



Parameters
------------------------------
Parameters are a convenient way to build a configurable and human-readable input XML.
They are defined via a block in the XML structure.
Parameter names may only include upper/lower case letters and underscores (to avoid conflicts with symbolic math).
Parameters may have any value:

- Path to a file
- Numbers
- A symbolic expression
- Other parameters
- Etc.


They can be used in any field within in the XML file (except in Includes) as follows:

- $x_par
- $:x_par
- $x_par$ 
- $:x_par$


For Example:

.. code-block:: xml

  <Parameters>
    <Parameter name='x' value='5'/>
    <Parameter name='y' value='5'/>
  </Parameters>
  <Partition>
    <SpatialPartition xPar='$x$' yPar='$y$' zPar='1'/>
  </Partition>


Units
------------------------------
By default, input values are specified using SI units.
In some cases, it is useful to override this behavior by explicitly specifying the units of the input.
These are specified by appending a valid number with a unit definition in square braces.
The unit manager supports most common units and SI prefixes, using both long- and abbreviated names (e.g.: c, centi, k, kilo, etc.).
Units may include predefined composite units (dyne, N, etc.) or may be built up from sub-units using a python syntax (e.g.: [N], [kg*m/s**2]).
Any (or no) amount of whitespace is allowed between the number and the unit bracket.


Examples:

.. code-block:: xml

  <Parameters>
    <Parameter name='a' value='2[m]'/>
    <Parameter name='b' value='1.2 [cm]'/>
    <Parameter name='c' value='1.23e4 [bbl/day]'/>
    <Parameter name='d' value='1.23E-4 [km**2]'/>
  </Parameters>



Symbolic Math
------------------------------
Input XML files can also include symbolic mathematical expressions.
These are placed within pairs of backticks (\`), and use a python syntax.
Parameters and units are evaluated before symbolic expressions.
Note: symbolic expressions are sanitized by removing any residual alpha characters, but this can be relaxed if more complicated function are needed.


Examples:

.. code-block:: xml

  <Parameters>
    <Parameter name='a' value='2[m]'/>
    <Parameter name='b' value='1.2 [cm]'/>
    <Parameter name='c' value='1.23e4 [bbl/day]'/>
    <Parameter name='d' value='1.23E-4 [km**2]'/>
  </Parameters>
  <Nodesets>
    <Nodeset name='perf' xmin='`$a$ - 0.2*$b$` -1e6 -1e6' xmax='`$c$**2 / $d$` 1e6 1e6' />
  </Nodesets>


Validation
------------------------------
Unmatched special characters ($, [, \`, etc.) mean that parameters, units, or symbolic math were not specified correctly.  
If the code detects these, it will throw an error.
The XML is validated against the input schema to check if all of the requireds field are present, and that input parameters match their expected types.


