.. _QuickStart:

###############################
Quick Start Guide
###############################

This page takes you from zero to a first GEOS simulation:

#. get the source code,
#. obtain the third-party libraries (TPLs) GEOS depends on,
#. configure and compile GEOS,
#. run a first simulation and the unit tests.

GEOS has no graphical user interface: you describe a problem in an XML input file (optionally together with a mesh file),
run the ``geos`` executable from a terminal or a batch job, and visualize the results with
`ParaView <https://www.paraview.org/>`_ or `VisIt <https://visit-dav.github.io/visit-website/>`_.
GEOS is currently distributed as source code only, so you will need to compile it.

.. note::
   The executable is called ``geos``. The build directory contains it as ``bin/geosx``;
   ``make install`` copies it to ``<install-dir>/bin/geosx`` and also creates an identical ``<install-dir>/bin/geos``.
   Both names refer to the same program.

At a glance
===========

For an experienced user on a Linux workstation with the :ref:`Prerequisites` installed, the whole process is:

.. code-block:: sh

   mkdir codes && cd codes
   git clone https://github.com/GEOS-DEV/GEOS.git
   git clone https://github.com/GEOS-DEV/thirdPartyLibs.git
   (cd GEOS && git submodule update --init --recursive)
   (cd thirdPartyLibs && git submodule update --init --recursive)

   # Write a host-config describing your compilers and libraries (see below)
   cp GEOS/host-configs/quick-start-template.cmake GEOS/host-configs/my-machine.cmake
   # ... edit GEOS/host-configs/my-machine.cmake ...

   # Third-party libraries (long: 30 minutes to a few hours)
   cd thirdPartyLibs
   python3 scripts/config-build.py -hc ../GEOS/host-configs/my-machine.cmake -bt Release
   make -C build-my-machine-release

   # GEOS
   cd ../GEOS
   python3 scripts/config-build.py -hc host-configs/my-machine.cmake -bt Release
   make -C build-my-machine-release -j 8 install

   # First run
   install-my-machine-release/bin/geos -i inputFiles/singlePhaseFlow/3D_10x10x10_compressible_smoke.xml -o first_run

The rest of this page explains each step.

Choose how to get the third-party libraries
===========================================

GEOS relies on several open-source third-party libraries (hypre, VTK, RAJA, HDF5, ...) that must be built with
compatible versions and options. There are three ways to get them:

* **Use a Docker image or a GitHub Codespace (fastest).**
  The GEOS team publishes Docker images with the TPLs already compiled, which are used by the continuous integration.
  You only compile GEOS itself. See :ref:`QuickStartDocker`.
* **Build the TPLs yourself (most common on workstations).**
  The `thirdPartyLibs <https://github.com/GEOS-DEV/thirdPartyLibs>`_ repository builds every library with tested versions in a single command.
  Follow the sections below, starting at :ref:`QuickStartDownload`.
* **Use a shared TPL installation (HPC centers).**
  On clusters with other GEOS users, the TPLs are often already installed in a shared location, and a matching host-config
  is provided in ``GEOS/host-configs/`` (for example for LLNL, ORNL, Stanford Sherlock and TotalEnergies systems).
  Ask your local GEOS contact, then skip directly to :ref:`QuickStartCompileGEOS`.

.. _QuickStartDocker:

Fastest path: Docker or GitHub Codespaces
=========================================

The images are named after their operating system and compiler, e.g. ``geosx/ubuntu24.04-gcc13``, and are hosted on
`Docker Hub <https://hub.docker.com/u/geosx>`_. The image tag that matches your version of GEOS is stored in
``.devcontainer/devcontainer.json`` (``build.args.GEOS_TPL_TAG``).

.. important::
   These images contain the compiled TPLs and the compilers, **not** a compiled ``geos`` executable.
   You still need to build GEOS inside the container (about 10 to 30 minutes).

**GitHub Codespaces / VS Code Dev Containers.**
The repository contains a ready-to-use ``.devcontainer`` configuration.
Open the GEOS repository on GitHub and choose *Code > Codespaces > Create codespace*,
or open your local clone in VS Code and select *Reopen in Container*.
The container comes with the CMake Tools extension preconfigured: configure and build from the VS Code CMake panel.

**Plain Docker.** From your GEOS clone (see :ref:`QuickStartDownload`):

.. code-block:: sh

   cd GEOS
   TAG=$(python3 -c "import json; print(json.load(open('.devcontainer/devcontainer.json'))['build']['args']['GEOS_TPL_TAG'])")
   docker run -it --rm -v "$PWD":/workspace -w /workspace \
              -e OMPI_ALLOW_RUN_AS_ROOT=1 -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
              geosx/ubuntu24.04-gcc13:${TAG} bash

   # Inside the container: the image provides a host-config at /spack-generated.cmake
   git config --global --add safe.directory '*'
   python3 scripts/config-build.py -hc /spack-generated.cmake -bt Release \
           -bp /workspace/build-docker -ip /workspace/install-docker \
           -DENABLE_HYPRE=ON -DENABLE_TRILINOS=OFF -DGEOS_LA_INTERFACE=Hypre
   make -C /workspace/build-docker -j $(nproc) install

The executable is then ``/workspace/install-docker/bin/geos``. Continue with :ref:`QuickStartFirstRun`.

.. _QuickStartDownload:

Download the source
===================

You need `git <https://git-scm.com>`_. Clone GEOS and thirdPartyLibs next to each other; the default
configuration expects this layout:

.. code-block:: sh

  codes/
  ├── GEOS/
  └── thirdPartyLibs/

.. code-block:: sh

   mkdir codes && cd codes
   git clone https://github.com/GEOS-DEV/GEOS.git
   cd GEOS
   git submodule update --init --recursive
   cd ..
   git clone https://github.com/GEOS-DEV/thirdPartyLibs.git
   cd thirdPartyLibs
   git submodule update --init --recursive
   cd ..

The GEOS submodules are ``src/cmake/blt`` (build system), ``src/coreComponents/LvArray`` (array library),
``src/coreComponents/fileIO/coupling/hdf5_interface`` and ``src/coreComponents/constitutive/HPCReact`` (chemistry package).
The thirdPartyLibs repository contains the build scripts; the library sources themselves are downloaded from the internet
when you build them, so the TPL build machine needs internet access.

.. note::
   If you use SSH keys with GitHub, replace ``https://github.com/`` by ``git@github.com:`` in the clone commands.
   ``git-lfs`` is not needed to build GEOS.

.. _QuickStartHostConfig:

Write a host-config file
========================

Before you start, install the :ref:`Prerequisites` (compilers, MPI, CMake 3.24 or newer, BLAS/LAPACK, Python 3).
It is useful to record the tools you use; this is the first thing the developers will ask for if something goes wrong:

.. code-block:: sh

  g++ --version
  mpicxx --version
  cmake --version

GEOS and its TPLs are configured with a CMake *host-config* file: a short CMake script that tells the build system
which compilers to use, where MPI and BLAS/LAPACK are, and which options to enable.
**Use the same host-config for the TPLs and for GEOS**, since options such as ``ENABLE_MPI`` or ``ENABLE_CUDA`` must match.

Start from ``host-configs/quick-start-template.cmake``:

.. literalinclude:: ../../../host-configs/quick-start-template.cmake
   :language: cmake

For example, on Ubuntu 24.04 (x86_64) with the packages listed in :ref:`Prerequisites`, copy the template to
``host-configs/my-machine.cmake`` and change the following lines:

.. code-block:: cmake

   set( CONFIG_NAME "my-machine" )

   set( CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "" )
   set( CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "" )

   set( MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "" )
   set( MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "" )
   set( MPIEXEC "/usr/bin/mpirun" CACHE PATH "" )

   set( BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so" CACHE PATH "" FORCE )
   set( LAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/liblapack.so" CACHE PATH "" FORCE )

.. warning::
   Keep ``CONFIG_NAME`` identical to the host-config file name (without ``.cmake``).
   ``config-build.py`` names the build and install directories after the *file name*
   (``build-my-machine-release``, ``install-my-machine-release``), while the template locates the TPLs using ``CONFIG_NAME``
   (``../thirdPartyLibs/install-${CONFIG_NAME}-release``). If the two differ, GEOS will not find the TPLs
   and CMake prints ``'GEOS_TPL_DIR' does not exist``.

Ready-made host-configs are provided in ``host-configs/`` for macOS (``apple/``, see :ref:`AppleMacOS`) and for several HPC systems
(``LLNL/``, ``ORNL/``, ``Stanford/``, ``TOTAL/``). Start from the closest one.
If you write a host-config for a platform that may be useful to others, please consider contributing it.

Compile the third-party libraries
=================================

.. code-block:: sh

   cd thirdPartyLibs
   python3 scripts/config-build.py -hc ../GEOS/host-configs/my-machine.cmake -bt Release
   cd build-my-machine-release
   make

The libraries are installed in ``thirdPartyLibs/install-my-machine-release``.
Depending on your machine, this takes from 30 minutes to a few hours.

.. note::
   Do not pass ``-j`` to this ``make``: the top-level build only launches the individual library builds, each of which already
   uses all the cores of the machine. To limit the number of cores, add ``-DNUM_PROC=4`` to the ``config-build.py`` command.

   A warning ``'GEOS_TPL_DIR' does not exist`` during this step is expected: the TPLs are not installed yet.

.. _QuickStartCompileGEOS:

Compile GEOS
============

.. code-block:: sh

   cd ../../GEOS
   python3 scripts/config-build.py -hc host-configs/my-machine.cmake -bt Release
   cd build-my-machine-release
   make -j 8
   make install

Adjust ``-j 8`` to the number of cores of your machine. A build takes 10 to 30 minutes.

If your TPLs are not in the default location (for example a shared installation), pass their path explicitly
(use a full path):

.. code-block:: sh

   python3 scripts/config-build.py -hc host-configs/my-machine.cmake -bt Release -DGEOS_TPL_DIR=/full/path/to/tpl/install

``config-build.py`` accepts other options (build type, build and install directories, Ninja generator, ...),
and any ``-D`` option is passed to CMake. See ``python3 scripts/config-build.py -h`` and :ref:`BuildProcess`.

.. note::
   Running ``config-build.py`` again **deletes** the existing build directory before configuring.
   To change a single option in an existing build, run ``cmake -D<OPTION>=<value> .`` in the build directory instead.

If all goes well, the executable is available at

.. code-block:: sh

  GEOS/install-my-machine-release/bin/geos

Check the executable
====================

.. code-block:: sh

  cd GEOS/install-my-machine-release
  ./bin/geos --help

This prints the available command-line options:

.. code-block:: none

    USAGE: geos -i input.xml [options]
           geos -s schema-output.xml

    Options:
    -?, --help
    -i, --input,             Input xml filename (required)
    -r, --restart,           Target restart filename
    -x, --xpartitions,       Number of partitions in the x-direction
    -y, --ypartitions,       Number of partitions in the y-direction
    -z, --zpartitions,       Number of partitions in the z-direction
    -s, --schema,            Name of the output schema
    -v, --validate-input,    Only do the loading phase, and not actual simulation. Useful to validate 'input'.
    -b, --use-nonblocking,   Use non-blocking MPI communication
    -n, --name,              Name of the problem, used for output
    --suppress-pinned,       Suppress usage of pinned memory for MPI communication buffers
    -o, --output,            Directory to put the output files
    -t, --timers,            String specifying the type of timer output
    --trace-data-migration,  Trace host-device data migration
    -m, --memory-usage,      Minimum threshold for printing out memory allocations in a member of the data repository.
    --pause-for,             Pause geos for a given number of seconds before starting execution
    -e, --errorsOutput,      Output path for the errors file (".yaml" supported)

.. _QuickStartFirstRun:

Run your first simulation
=========================

The repository contains many ready-to-run input files in ``inputFiles/``. From the ``GEOS`` directory, run a small
single-phase flow problem (a 10x10x10 mesh, a few seconds on a laptop):

.. code-block:: sh

   cd GEOS
   install-my-machine-release/bin/geos -i inputFiles/singlePhaseFlow/3D_10x10x10_compressible_smoke.xml -o first_run

GEOS prints a summary of the problem, then one block per time step with the nonlinear (Newton) and linear solver iterations,
and ends with timing information. The results are written to the directory given with ``-o`` (here ``first_run``;
the default is the current directory). This input file requests Silo output, written to ``first_run/siloFiles``
and readable with VisIt; most other examples write VTK files (``.pvd``) that open directly in ParaView.

To check an input file without running the simulation, use ``-v`` (``--validate-input``).

To run in parallel with MPI, split the domain into ``-x``, ``-y`` and ``-z`` partitions. The product of the three
numbers must equal the number of MPI ranks:

.. code-block:: sh

   mpirun -np 8 install-my-machine-release/bin/geos -i inputFiles/singlePhaseFlow/3D_10x10x10_compressible_smoke.xml \
          -x 2 -y 2 -z 2 -o first_run_parallel

See :ref:`RunningGEOS` for all command-line options, parallel runs, restarts and reading the log.

Run the unit tests
==================

Running the unit tests is a good check of your build:

.. code-block:: sh

  cd GEOS/build-my-machine-release
  ctest --output-on-failure -j 8

The full suite takes several minutes. If *all* tests fail, something is likely wrong with your installation
(see the troubleshooting section below). If only a few fail, your platform may have exposed an issue that our continuous
integration does not cover: please report it on the issue tracker. Developers can also run the integrated tests, see :ref:`IntegratedTests`.

Next steps
==========

* Work through the :ref:`Tutorials`, starting with :ref:`TutorialSinglePhaseFlowWithInternalMesh`: they explain the
  structure of an input file step by step.
* Browse the :ref:`BasicExamples` organized by physics (multiphase flow, wells, CO\ :sub:`2` injection, poromechanics, hydraulic fracturing).
* Look up any XML element or attribute in the :ref:`UserGuide` and in the complete input schema (:doc:`Datastructure Index </docs/sphinx/datastructure/CompleteXMLSchema>`).

Frequently asked questions
==========================

**Do I need to be a code developer to use GEOS?**
You need to be comfortable with compiling software from source, but not with C++ development.
Once installed, GEOS is driven entirely by XML input files.
GEOS is an open-source project, and contributions are welcome (see the :ref:`Contributing` section of the developer guide).

**What are the system requirements?**
GEOS runs on Linux and macOS, from laptops to the largest supercomputers, on CPUs and on NVIDIA or AMD GPUs.
On Windows, use a virtual machine, WSL, or a Docker container (see :ref:`InstallWin`).
GPU builds require additional software (CUDA or ROCm), see the :ref:`BuildGuide`.

**Where is the source code?**
All repositories are hosted in the `GEOS-DEV GitHub organization <https://github.com/GEOS-DEV>`_:
`GEOS <https://github.com/GEOS-DEV/GEOS>`_ itself, `thirdPartyLibs <https://github.com/GEOS-DEV/thirdPartyLibs>`_ (TPL build scripts),
first-party libraries included as git submodules (`LvArray <https://github.com/GEOS-DEV/LvArray>`_, `HPCReact <https://github.com/GEOS-DEV/HPCReact>`_, ...),
and `geosPythonPackages <https://github.com/GEOS-DEV/geosPythonPackages>`_ (Python tools for pre- and post-processing).

**I get errors while downloading, compiling or running.**
Unfortunately, no set of instructions is foolproof. We recommend:

#. Re-read the instructions and the error message. When a compilation fails, look for the *first* error: later ones are often consequences.
#. Search for the error message online; a library-specific build problem may be answered by that library's support pages.
#. Search the `GEOS issue tracker <https://github.com/GEOS-DEV/GEOS/issues>`_, including closed issues (remove the default ``is:open`` filter).
#. If you are still stuck, open a new issue with your operating system, compiler, MPI and CMake versions, your host-config,
   and the full error message. We may ask you to run additional tests.
