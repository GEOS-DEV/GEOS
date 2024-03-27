################################################################################
Sphinx Documentation
################################################################################

Generating the documentation
====================================

- To generate the documentation files, you will need to install Sphinx using :

  .. code-block:: sh

    pip install sphinx
    pip install sphinx-design sphinx-argparse sphinxcontrib-plantuml sphinxcontrib.programoutput sphinx_rtd_theme
    pip install scipy

- Then you can generate the documentation files with the following command :

  .. code-block:: sh

    cd /path/to/GEOS/build-your-platform-release
    make geosx_docs

- That will create a new folder

  .. code-block:: sh

    /path/to/GEOS/build-your-platform-release/html/docs/sphinx

which contains all the html files generated.

Documenting the code
====================================

The documentation is generated from restructured text files (``.rst``). Most files
can be found in ``src/docs/sphinx``. Files which are specific to parts of the code,
like those describing a specific class, can instead be found in ``docs`` subdirectory
in the folder containing the source code.

Information about how to write ``rst`` files can be found `here <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_ .

Generated ``.rst`` tables
====================================

To describe the GEOS classes, some parts of the Sphinx documentation are automatically filled with
generated ``.rst`` tables. Those tables are filled by the registered ``Group`` & ``Wrapper``
configuration (with ``Group::registerWrapper()``, ``Wrapper::setDescription()``, ...).

When modifying a ``Wrapper`` configuration, one should call ``make geosx_update_rst_tables`` to
generate the tables, and should then include the updated / added tables in its PR, exactly the
same way the ``schema.xsd`` is maintained with ``make geosx_generate_schema``.

To make these commands work, the folowing line is needed in the host-config :

  .. code-block:: sh

    set( ENABLE_XML_UPDATES ON CACHE BOOL "" FORCE ) 