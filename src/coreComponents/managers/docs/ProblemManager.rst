###############################################################################
Problem Manager
###############################################################################


Command Line Interface
=========================

The command-line options for GEOSX include:

* ``-i, --input`` - Input .xml filename (required)
* ``-r, --restart`` - Target restart filename
* ``-x, --x-partitions`` - Number of partitions in the x-direction
* ``-y, --y-partitions`` - Number of partitions in the y-direction
* ``-z, --z-partitions`` - Number of partitions in the z-direction
* ``-s, --schema`` - Name of the schema file to generate
* ``-l, --schema-level`` - Verbosity level of output schema (default=0)
* ``-o, --output`` - Directory to place output files


Input Schema Generation
===========================

A schema file is a useful tool for validating input .xml files and constructing user-interfaces.  Rather than manually maintaining the schema during development, GEOSX is designed to automatically generate one by traversing the documentation structure.

To generate the schema, run GEOSX with the input, schema, and the (optional) schema_level arguments, i.e.: ``geosx -i input.xml -s schema.xsd``.  There are two ways to limit the scope of the schema:

1. Setting the verbosity flag for an object in the documentation structure.  If the schema-level argument is used, then only objects (and their children) and attributes with ``(verbosity < schema-level)`` will be output.

2. By supplying a limited input xml file.  When GEOSX builds its data structure, it will only include objects that are listed within the xml (or those that are explicitly appended when those objects are initialized).  The code will add all available *attributes* for these objects to the schema.

To take advantage of this design it is necessary to use the automatic xml parsing approach that relies upon the documentation node.  If values are read in manually, then the schema can not be used to validate xml those inputs.  

Note: the lightweight xml parser that is used in GEOSX cannot be used to validate inputs with the schema directly.  As such, it is necessary to use an external tool for validation, such as the geosx_tools python module.



