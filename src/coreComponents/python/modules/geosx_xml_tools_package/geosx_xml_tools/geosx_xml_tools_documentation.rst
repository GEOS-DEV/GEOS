
.. _advanced_xml_features:

###############################################################################
Advanced XML Features (geosx_xml_tools)
###############################################################################

geosx_xml_tools is a python module that enables advanced xml features in GEOSX (parameters, units, symbolic math, etc.), and is used to format xml files.


Setup
=================================

The `geosx_xml_tools` package can be installed using the following command within the GEOSX build directory:

`make geosx_xml_tools`

During the installation step, two console scripts will be created: `preprocess_xml` and `format_xml`.
These will be located within the GEOSX build/bin directory.
Additional things to consider: 

- The above make command will create a new virtual python environment to install the package.  By default, the source python environment will be the the version used to run the `config-build.py` command.  The version of python can be changed by specifying the `PYTHON_POST_EXECUTABLE` variable.  To build the environment, the `virtualenv` package must be installed in the source distribution.
- The `geosx_xml_tools` package depends on the `lxml` package. If `lxml` is missing from the parent environment, the install script will attempt to fetch an appropriate version from the internet.
- The package may also be manually installed within an existing python distribution (this required administrative priviliges) via pip: `pip install src/coreComponents/python/modules/geosx_xml_tools_package`.  In this case, the console scripts will be located in the python/bin directory


Usage
=================================

geosx_xml_tools can be used via the command-line or can be imported into a python script.


Command-line xml formatting
------------------------------

The following command will update the formatting of an existing xml file in-place:

`format_xml input_file.xml`

To update the formatting for all xml files located within the GEOSX repository, you can run the following command within the GEOSX build directory:

`make geosx_format_all_xml_files`


Command-line xml preprocessing
------------------------------

The following command will read an xml file, process any advanced xml features located within it, and write a new file that can be read by GEOSX:

`preprocess_xml input_file.xml`

The script returns the (randomly generated) name of the new xml file.
Optional arguments for this script include:

- `-o/--output`: The desired name for the output file (otherwise, it is randomly generated)
- `-s/--schema`: The location of a schema to validate the final .xml file
- `-v/--verbose`: Increase module verbosity

For convenience, this script can be embedded within a call to GEOSX:

`srun -n 16 geosx -i \`preprocess_xml input_file.xml\` -x 4 -z 4`


Script-based Example
------------------------------

The geosx_xml_tools module can also be called from within a python script.  For example:

.. code-block:: python

  from geosx_xml_tools import xml_processor

  initial_filename = 'input.xml'
  new_filename = output_name = xml_processor.process(initial_filename)



Advanced XML Features
=================================

The xml preprocessor in geosx_xml_tools is designed to take a raw input file, and generate a new file that can be directly read by GEOSX.
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

- $x_par$  (preferred)
- $x_par
- $:x_par
- $:x_par$


For Example:

.. code-block:: xml

  <Parameters>
    <Parameter
      name="mu"
      value="0.005"/>
    <Parameter
      name="table_root"
      value="/path/to/table/root"/>
  </Parameters>
  
  <Constitutive>
    <CompressibleSinglePhaseFluid
      name="water"
      defaultDensity="1000"
      defaultViscosity="$mu$"
      referencePressure="0.0"
      referenceDensity="1000"
      compressibility="5e-10"
      referenceViscosity="$mu$"
      viscosibility="0.0"/>
  </Constitutive>

  <Functions>
    <TableFunction
      name="flow_rate"
      inputVarNames="{time}"
      coordinateFiles="{$table_root$/time_flow.geos}"
      voxelFile="$table_root$/flow.geos"
      interpolation="linear"/>
  </Functions>



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
Note: symbolic expressions are sanitized before by removing any residual alpha characters, but this can be relaxed if more complicated function are needed.
Also, while symbolic expressions are allowed within parameters, errors may occur if these are used in a way that leads to nested symbolic expressions.


Examples:

.. code-block:: xml

  <Parameters>
    <Parameter name='a' value='2[m]'/>
    <Parameter name='b' value='1.2 [cm]'/>
    <Parameter name='c' value='1.23e4 [bbl/day]'/>
    <Parameter name='d' value='1.23E-4 [km**2]'/>
  </Parameters>
  <Geometry>
    <Box
      name='perf'
      xMin='`$a$ - 0.2*$b$`, -1e6, -1e6'
      xMax='`$c$**2 / $d$`, 1e6, 1e6' />
  </Geometry>


Validation
------------------------------
Unmatched special characters ($, [, \`, etc.) mean that parameters, units, or symbolic math were not specified correctly.  
If the code detects these, it will throw an error.
The XML is validated against the input schema to check if all of the requireds field are present, and that input parameters match their expected types.


