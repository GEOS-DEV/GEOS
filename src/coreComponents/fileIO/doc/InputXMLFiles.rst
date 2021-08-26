###############################################################################
Input Files
###############################################################################


XML
=================================

GEOSX is configured via one (or more) `Extensible Markup Language <https://en.wikipedia.org/wiki/XML>`_ (XML) files.
These files contain a set of elements and attributes that closely follow the internal datastructure of GEOSX.
When running GEOSX, these files are specified using the `-i` argument:

.. code-block:: bash

    geosx -i input.xml



XML Components
------------------------------

The following illustrates some of the key features of a GEOSX-format xml file:

.. code-block:: xml

    <?xml version="1.0" ?>

    <Problem>
        <BlockA
            someAttribute="1.234">

            <!-- Some comment -->
            <BlockB
              name="firstNamedBlock"
              anotherAttribute="0"/>
            <BlockB
              name="secondNamedBlock"
              anotherAttribute="1"/>
        </BlockA>
    </Problem>


The two basic components of an xml file are blocks, which are specified using angle brackets ("<BlockA>  </BlockA>"), and attributes that are attached to blocks (attributeName="attributeValue").
Block and attributes can use any Ascii character aside from `<`, `&`, `'`, and `"` (if necessary, use `&lt;`, `&amp;`, `&apos;`, or `&quot;`).
Comments are indicated as follows: `<!-- Some comment -->`.

At the beginning of a GEOSX input file, you will find an optional xml declaration (`<?xml version="1.0" ?>`) that is used to indicate the format to certain text editors.
You will also find the root `Problem` block, where the GEOSX configuration is placed.
Note that, aside from these elements and commented text, the xml format requires that no other objects exist at the first level.

In the example above, there is a single element within the `Problem` block: `BlockA`.
`BlockA` has an attribute `someAttribute`, which has a value of 1.234, and has three children: a comented string "Some comment" and two instances of `BlockB`.
The `name` attribute is required for blocks that allow multiple instances, and should include a unique string to avoid potential errors.
Where applicable these blocks will be executed in the order in which they are specified in input file.


Input Validation
=================================

The optional `xmlns:xsi` and `xsi:noNamespaceSchemaLocation` attributes in the Problem block can be used to indicate the type of document and the location of the xml schema to the text editor:

.. code-block:: xml

    <Problem
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xsi:noNamespaceSchemaLocation="/path/to/schema.xsd" />

The schema contains a list of xml blocks and attributes that are supported by GEOSX, indicates whether a given object is optional or required, and defines the format of the object (string, floating point number, etc.).
A copy of the schema is included in the GEOSX source code (/path/to/GEOSX/src/coreComponents/schema/schema.xsd).
It can also be generated using GEOSX: ``geosx -s schema.xsd``

Many text editors can use the schema to help in the construction of an xml file and to indicate whether it is valid.
Using a validation tool is highly reccomended for all users.
The following instructions indicate how to turn on validation for a variety of tools:


xmllint
---------------------------------

xmllint is a command-line tool that is typically pre-installed on unix-like systems.
To check whether an input file is valid, run the following command:

xmllint --schema /path/to/schema.xsd input_file.xml


Sublime Text
------------------------------

We reccomend using the `Exalt <https://github.com/eerohele/exalt>`_ or `SublimeLinter_xmllint <https://github.com/SublimeLinter/SublimeLinter-xmllint>`_ plug-ins to validate xml files within sublime.
If you have not done so already, install the sublime `Package Control <https://packagecontrol.io/installation>`_.
To install the package, press ``ctrl + shift + p``, type and select ``Package Controll: Install Package``, and search for ``exalt`` or ``SublimeLinter`` / ``SublimeLinter-xmllint``.
Note that, depending on the circumstances, these tools may indicate only a subset of the validation errors at a given time.
Once resolved, the tools should re-check the document to look for any additional errors.

As an additional step for SublimLinter-xmllint, you will need to add a linter configuration.
To do so, go to Preferences/Package Settings/SublimeLinter/Settings.
In the right-hand side of the new window, add the xmllint configuration:

.. code-block:: python

    {
        "linters": {
            "xmllint":
            {
                "args": "--schema /path/to/schema.xsd",
                "styles": [
                    {
                        "mark_style": "fill",
                        "scope": "region.bluish",
                        "types": ["error"],
                        "icon": "stop",
                    }
                ]
            },
        }
    }



Eclipse
------------------------------

To Eclipse Web Develop Tools includes features for validating xml files.
To install them, go to Help -> Eclipse Marketplace, search for the Eclipse Web Developer Tools, install the package, and restart eclipse.
Finally, configure the xml validation preferences under Window -> Preferences -> XML -> XML Files -> Validation.
Eclipse will automatically fetch the schema, and validate an active xml file.
The editor will highlight any lines with errors, and underline the specific errors.


GEOSX XML Tools
------------------------------

The geosx_xml_tools package, which is used to enable advanced features such as parameters, symbolic math, etc., contains tools for validating xml files.
To do so, call the command-line script with the -s argument, i.e.: `preprocess_xml input_file.xml -s /path/to/schema.xsd`.
After compiling the final xml file, pygeos will fetch the designated schema, validate, and print any errors to the screen.

Note: Attributes that are using advanced xml features will likely contain characters that are not allowed by their corresponding type pattern.
As such, file editors that are configured to use other validation methods will likely identify errors in the raw input file.


XML Schema
=================================

An XML schema definition (XSD) file lays out the expected structure of an input XML file.
During the build process, GEOSX automatically constructs a comprehensive schema from the code's data structure, and updates the version in the source (GEOSX/src/coreComponents/schema/schema.xsd).


Schema Components
------------------------------

The first entry in the schema are a set of headers the file type and version.
Following this, the set of available simple types for attributes are layed out.
Each of these include a variable type name, which mirrors those used in the main code, and a regular expression, which is designed to match valid inputs.
These patterns are defined and documented in ``DataTypes::typeRegex``.
The final part of the schema is the file layout, beginning with the root ``Problem``.
Each complex type defines an element, its children, and its attributes.
Each attribute defines the input name, type, default value, and/or usage.
Comments preceding each attribute are used to relay additional information to the user.


Automatic Schema Generation
------------------------------

A schema may be generated by calling the main code with the -s argument , e.g.: ``geosx -s schema.xsd`` (Note: this is done automatically during the bulid process).
To do this, GEOSX does the following:

  1) Initialize the GEOSX data structure.
  2) Initialize objects that are registered to catalogs via ``ManagedGroup::ExpandObjectCatalogs()``.
  3) Recursively write element and attribute definitions to the schema using information stored in GEOSX groups and wrappers.
  4) Define any expected deviations from the schema via ``ManagedGroup::SetSchemaDeviations()``.



