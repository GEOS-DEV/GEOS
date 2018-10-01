###############################################################################
Getting Started with GEOSX
###############################################################################

Getting Ready
=================================
GEOSX resides in a git repository hosted at https://github.com/GEOSX/GEOSX.
It is suggested that you setup ssh keys, and use ssh for your clones as discussed 
`here <https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/>`__.
It is suggested that a directory be created to host the various clones that one may require for an effective development workflow, so the first step is to make a directory

1. Setup working directory

.. code-block:: sh

  mkdir geosx
  cd geosx

Downloading the Code
=================================
There are currently two separate repositories that should be downloaded.
The first is the main repository, which may be cloned and initialized by the following steps: 

2. Clone the main repository

.. code-block:: sh

   git clone git@github.com:GEOSX/GEOSX.git
   
   cd GEOSX
   
   git submodule init
   
   git submodule update
   
   cd ..


3. Clone the third-party libraries

.. code-block:: sh

   git clone git@github.com:GEOSX/thirdPartyLibs.git
   cd thirdPartyLibs
   git submodule init
   git submodule update
   cd ..


Compiling the Code
=================================

GEOSX compilations are typically driven by a hostconfig file, which reside in GEOSX/host-configs.
If your platform does not have a host-config in the repository, you are encouraged to maintain one.
If you are running on an LC system, there is already a hostconfig and copy of the thirdPartyLibs installed, and you can skip step 4.

The first step in compiling GEOSX is to run cmake and generate the makefiles.
Starting with the third-party libraries, the config-build.script will run cmake for you.
Note that the 'make' step should be run serially, as the indiviudal package builds are run in parallel by default.

4. Configure and make the third party libraries

.. code-block:: sh

   cd thirdPartyLibs
   python scripts/config-build.py -hc ../GEOSX/host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make hdf5
   make -j1

The next step is to compile the main code. 
Again, the config-build sets up cmake for you.

5. Configure and make the main code

.. code-block:: sh

   cd ../../GEOSX
   python scripts/config-build.py -hc host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make -j4

   
Running the Code
=================================

GEOSX executables read in a XML input file. A simple example XML is located
`here <https://github.com/GEOSX/GEOSX/blob/develop/src/components/core/tests/PhysicsSolvers/LaplaceFEM.xml/>`__. 
To execute a serial run enter the following command from a working directory:

.. code-block:: sh

    path-to-geosx-bin/geosx -i path-to-xml/LaplaceFEM.xml
