###############################################################################
Python Input Processing: pygeos
###############################################################################

Pygeos is a python-based module designed to construct input files for GEOSX.
It allows users to create xml files that include child files, parameters, units, and symbolic math.
It also includes helper-functions for validating input files and constructing tables. 


Setup
=================================
The pygeos module is located here: `src/coreComponents/python/modules/pygeos_package` .
The module is compabitable with python 2/3, and its dependencies include the `numpy` and `lxml` packages.
It can easily be installed using tools such as pip:

`pip install src/coreComponents/python/modules/pygeos_package`


Virtual Python Environment
---------------------------------

On systems with a shared python installation, we reccomend that users install the package within a virtual python environment.
A script designed to automatically setup a virtual environment, install the dependencies, and install pygeos is included in `scripts/setupVirtualPythonEnvironment.bash` .
The options for the script include:

- -p/--python_root : This specifies the root path for the parent python environment
- -o/--output_path : This specifies location where the virtual environment will be installed

The default settings for the script will build a virtual environment for the python-3.6.4 installation on LC systems, and will place the environment in `~/Python/virtual/geosx`

To use the load the python environment, run the following: `source ~/Python/virtual/geosx/bin/activate` .
Within the environment, the commands `python` and `pip` will point to the correct versions.
To exit the virtual environment, run the command: `deactivate` .



XML Preprocessing
=================================

The xml preprocessor is designed to take an raw input file, and generate an new file that can be directly read by GEOSX.
The syntax for the advanced xml format is given below.
During the processing the order of operations are:

1) Merging any included xml files into the root structure
2) Substituting in any parameters
3) Evaluating unit strings
4) Evaluating symbolic math
5) Error checking and validation



Use Example
------------------------------

`pygeos.preprocessGEOSXML(input, schema='')` is used to process the input xml.
The name of the newly generated, preprocessed xml file will be returned by this function call.
By default it will have the name 'prep_' + a randomly generated string.

The following is a example python script that will read process an xml file specified via the command line.

.. code-block:: python

  import sys
  import pygeos

  new_filename = pygeos.preprocessGEOSXML(sys.argv[1])
  print(new_filename)



Including Child XML Files
------------------------------
XML inputs can point to included children (these children can then include grandchildren and so on).
During processing, these are recursively inserted into the root XML structure by the following scheme:

- Merge two objects if:
    - At the root level and an object with the matching tag exists.
    - If the “name” attribute is present and a object with the matching tag and name exist.
    - Any preexisting attributes are overwritten by the donor.
- Otherwise append the xml structure with the target.


.. code-block:: xml

  <Included>
    <File name='/path/to/included_a.xml'/>
    <File name='/path/to/included_b.xml'/>
  </Included>



Parameters
------------------------------
Parameters are a convenient way to build a configurable and human-readable input xml.
They are defined via a block in the xml structure.
Parameter names may only include upper/lower case letters and underscores (to avoid conflicts with symbolic math).
Parameters may have any value:

- Path to a file
- Numbers
- A symbolic expression
- Other parameters
- Etc.


They can be used in any field within in the xml file (except in Includes) as follows:

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
The unit manager supports most common units and SI prefixes, using both long- and abbreviated names (e.g.: c, centi, k, kilo, etc.)
Units may include predefined composite units (dyne, N, etc.) or may be built up from sub-units using a python syntax (e.g.: [N], [kg*m/s**2].
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
Input xml files can also include symbolic mathematical expressions.
These are indicated with curly braces, and use a python syntax.
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
    <Nodeset name='perf' xmin='{$a$ - 0.2*$b$} -1e6 -1e6' xmax='{$c$**2 / $d$} 1e6 1e6' />
  </Nodesets>


Validation
------------------------------
Unmatched special characters ($, [, }, etc.) mean that parameters, units, or symbolic math were not specified correctly.  
If the code detects these, it will throw an error.
The XML is validated against the input schema to check if all of the required field are present, and that input parameters match their expected types.



Input Table Generation
=================================

The pygeos package also includes some tools to convert numpy arrays into GEOSX input tables, and vice-versa.

The following example shows how to write a series of structrued tables using a list of spatial axes and a dictionary of table values:

.. code-block:: python

  import numpy as np
  import pygeos

  # Config
  N = (10, 20, 30)

  # Generate coordinate axes
  # These correspond to each axis of the numpy tables
  # The example function accepts up to four axes, and
  # will assign them the names x, y, z, and t in order
  spatial = [np.linspace(0, 1, N[0]),
             np.linspace(0, 1, N[1]),
             np.linspace(0, 1, N[2])]

  # Generate the property dictionary
  properties = {'random_variable_A': np.randn(N),
                'random_variable_B': np.randn(N),
                'random_variable_C': np.randn(N)}

  # Write the tables
  # The files will be written to the current directory,
  # and will have the .txt extension
  pygeos.writeGEOSTable(spatial, properties)


The following shows how to read the geos tables written in the previous example:

. code-block:: python

  import numpy as np
  import pygeos

  spatial_names = ['x', 'y', 'z']
  property_names = ['random_variable_A', 'random_variable_B', 'random_variable_C']

  spatial, properties = pygeos.readGEOSTable(spatial_names, property_names)




