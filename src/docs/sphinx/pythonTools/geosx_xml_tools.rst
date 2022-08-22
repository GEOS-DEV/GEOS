
GEOSX XML Tools
--------------------------

The `geosx_xml_tools` python package adds a set of advanced features to the GEOSX xml format: units, parameters, and symbolic expressions.
See :ref:`PythonToolsSetup` for details on setup instructions, and :ref:`AdvancedXMLFeatures` for a detailed description of the input format.
The available console scripts for this package and its API are described below.


preprocess_xml
^^^^^^^^^^^^^^

Compile an xml file with advanced features into a single file that can be read by GEOSX.
Note: this script will be automatically applied if you use `geosx_preprocessed` instead of `geosx` when calling the code.

.. argparse::
   :module: geosx_xml_tools.main
   :func: build_preprocessor_input_parser
   :prog: preprocess_xml


format_xml
^^^^^^^^^^^^^^

Formats an xml file.

.. argparse::
   :module: geosx_xml_tools.xml_formatter
   :func: build_xml_formatter_input_parser
   :prog: format_xml


check_xml_attribute_coverage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Checks xml attribute coverage for files in the GEOSX repository.

.. argparse::
   :module: geosx_xml_tools.attribute_coverage
   :func: build_attribute_coverage_input_parser
   :prog: check_xml_attribute_coverage


check_xml_redundancy
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Checks for redundant attribute definitions in an xml file, such as those that duplicate the default value.

.. argparse::
   :module: geosx_xml_tools.xml_redundancy_check
   :func: build_xml_redundancy_input_parser
   :prog: check_xml_redundancy


API
^^^

.. automodule:: geosx_xml_tools.main
    :members:

.. automodule:: geosx_xml_tools.xml_processor
    :members:

.. automodule:: geosx_xml_tools.xml_formatter
    :members:

.. automodule:: geosx_xml_tools.unit_manager
    :members:

.. automodule:: geosx_xml_tools.regex_tools
    :members:

.. automodule:: geosx_xml_tools.xml_redundancy_check
    :members:

.. automodule:: geosx_xml_tools.attribute_coverage
    :members:

.. automodule:: geosx_xml_tools.table_generator
    :members:

