
Python Tools
==========================


.. _PythonToolsSetup:

Python Tools Setup
---------------------------------

The preferred method to setup the GEOSX python tools is to run the following command in the build directory:

.. code-block:: bash

    make geosx_python_tools


This will attempt to install the required packages into the python distribution indicated via the `Python3_EXECUTABLE` cmake variable (also used by pygeosx).

If the user does not have write access for the target python distribution, the installation will attempt to create a new virtual python environment (Note: this requires that the virtualenv package be installed).
If any package dependencies are missing, then the install script will attempt to fetch them from the internet using pip.
After installation, these packages will be available for import within the associated python distribution, and a set of console scripts will be available within the GEOSX build bin directory.

Alternatively, these packages can be installed manually into a python environment using pip:

.. code-block:: bash

    cd GEOSX/src/coreComponents/python/modules/geosx_mesh_tools_package
    pip install --upgrade .

    cd ../geosx_xml_tools_package
    pip install --upgrade .

    # Etc.


Packages
-----------------------


.. toctree::
    :maxdepth: 1

    hdf5_wrapper

    geosx_mesh_tools

    geosx_xml_tools

    pygeosx_tools

    timehistory

    mesh_doctor
