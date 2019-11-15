Getting and compiling the codes
###############################

Of the repositories organisation of GEOSX
=========================================

In order to benefit from other's expertise, GEOSX makes extensive use of third part developments.
Through both *git submodules* and *external libraries*.

The source code of GEOSX is hosted in a git repository at https://github.com/GEOSX/GEOSX.
It contains links (through `git submodules <https://git-scm.com/docs/git-submodule>`_) to multiple other sources codes.
Let us for example cite cmake build extension for HPC `blt <https://github.com/LLNL/blt>`_.
The up to date list of course being the `.gitmodules <https://github.com/GEOSX/GEOSX/blob/develop/.gitmodules>`_ file itself.

Other kinds of dependencies (see :doc:`Dependencies<./Dependencies>`) are available in the `third party repository <https://github.com/GEOSX/thirdPartyLibs.git>`_.
This repository mirrors most of other third part libraries.
It uses git lfs (large file storage) to store and serve tarball of some dependencies (see the `tplMirror folder <https://github.com/GEOSX/thirdPartyLibs/tree/master/tplMirror>`_).
The git lfs settings are in the `.gitatributes <https://github.com/GEOSX/thirdPartyLibs/blob/master/.gitattributes>`_ file.

These elements explain why you should be careful about git submodules.
And for git lfs too, while it's more transparent for a daily usage.

Getting your acces working
==========================

It is suggested that you setup ssh keys, and use ssh for your clones as discussed
`on this github article <https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/>`_.
Before going further, you should `test your ssh connection <https://help.github.com/en/github/authenticating-to-github/testing-your-ssh-connection>`_.

If it fails (maybe because of your company proxy...), you may consider the `personnal access token option <https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line>`_.

The next command lines will use the ssh connection pattern.
Note that you'll have to adapt the commands if you adopt https connections.
The clone `git clone git@github.com:GEOSX/GEOSX.git` becomes `git clone https://${USER}:${TOKEN}@github.com/GEOSX/GEOSX.git`.
(You may or may not include your username and passwords in the URL).

Downloading the code
====================

You need to clone the GEOSX and the thirdPartyLibs repositories.
If you do not specify anything specific, the building command lines expect the GEOSX and thirdPartyLibs to be next to each other.
You can of course make another choice, but as a first attempt we do not recommend it.

First (for example), GEOSX itself...

.. code-block:: sh

   git clone git@github.com:GEOSX/GEOSX.git
   cd GEOSX
   git lfs install
   git submodule init
   git submodule update
   cd ..

... and second the third party libraries.

.. code-block:: sh

   git clone git@github.com:GEOSX/thirdPartyLibs.git
   cd thirdPartyLibs
   git lfs install
   git pull
   git submodule init
   git submodule update
   cd ..

Note that git-lfs may not funct-on properly (or may be very slow) if version of git and git-lfs are not current.
If you are using an older version of git/git-lfs you may need to add "git lfs pull" after "git pull" in the above procedures.

Compiling the Code
==================

GEOSX compilations are typically driven by a hostconfig file, which reside in GEOSX/host-configs.
If your platform does not have a host-config in the repository, you are encouraged to maintain one.
A good starting point shall be the `host-configs/environment.cmake` for both `GEOSX <https://github.com/GEOSX/thirdPartyLibs/blob/master/host-configs/environment.cmake>`_ and `third party libs <https://github.com/GEOSX/GEOSX/blob/develop/host-configs/environment.cmake>`_.
They are very similar, and you may whether define the dedicated (classical) environment variables (`CC`, `CXX`, ...) or replace the field directly in the file. 

(If you are running on an LC system, there is already a hostconfig and copy of the thirdPartyLibs installed, so you can skip the third party libraries compilation step).

Starting with the third-party libraries, the config-build.script will run cmake for you.
More details about compiling the thirdPartyLibs are provided in :ref:`Third party libraries build management pattern`.

.. code-block:: sh

   cd thirdPartyLibs
   python scripts/config-build.py -hc host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make

The next step is to compile the main code.
Again, the config-build sets up cmake for you.

.. code-block:: sh

   cd ../../GEOSX
   python scripts/config-build.py -hc host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make -j4
   make install


Running the tests
=================

It is wise to run the unit tests. You can run them in the build folder GEOSX: `ctest -V`.
Running the integrated tests shall be considered to. This task it a little more complex, so please refer to :ref:`Integrated Tests` for further information.

Running GEOSX
=============

GEOSX executables read in a XML input file. A simple Laplace example XML is located `in the integratedTests submodule <https://github.com/GEOSX/integratedTests/blob/develop/update/run/SimpleSolvers/10x10x10_LaplaceFEM.xml>`_
(if you've used git lfs properly, it should be in the `integratedTests/update/run/SimpleSolvers` folder).

To execute a serial run enter the following command from a working directory:

.. code-block:: sh

    path-to-geosx-bin/geosx -i path-to-xml/LaplaceFEM.xml

If you want to run the same test case in parallel, run the following mpi command:

.. code-block:: sh

    mpirun -np 8 path-to-geosx-bin/geosx -i path-to-xml/LaplaceFEM.xml -x 2 -y 2 -z 2

where you split the domain by 2 il all x, y and z directions (the total amount of mpi cores being the number of plit domains).

The temperature field should look like

.. image:: Laplace10x10x10-temperature.png
   :width: 400

if you want to write VTK output files, you can add

.. code-block:: xml

    <PeriodicEvent name="outputs"
                   timeFrequency="1.0"
                   targetExactTimestep="1"
                   target="/Outputs/vtkOutput"/>

to the `Events` node and 

.. code-block:: xml

    <VTK name="vtkOutput"/>

to the `Outputs` node.

Command Line options
--------------------

The following command-line options can be specified:

  * ``-?, --help`` - Provides a list of command-line options
  * ``-i, --input`` - Input filename.xml  (required)
  * ``-l, --schema-level`` - Verbosity level of output schema (default=0)
  * ``-n, --name`` - Name of the problem (used for output)
  * ``-o, --output`` - Directory to place output files
  * ``-r, --restart`` - Target restart filename
  * ``-s, --schema`` - Name of the schema file to generate
  * ``-x, --x-partitions`` - Number of partitions in the x-direction
  * ``-y, --y-partitions`` - Number of partitions in the y-direction
  * ``-z, --z-partitions`` - Number of partitions in the z-direction
