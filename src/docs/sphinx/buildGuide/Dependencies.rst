.. _Dependencies:

Third-party dependencies
========================

GEOS makes use of multiple third-party libraries (TPLs) and tools, some of which are mandatory and some optional.
We only test against specific versions, and sometimes even require development snapshots (specific git commits).
Not all of these guarantee backwards compatibility, so we strongly recommend building with these specific versions.

List of third-party libraries and tools
---------------------------------------

The tables below list the dependencies and the relevant CMake variables.
The exact versions that GEOS is tested with are defined in the
`thirdPartyLibs CMake script <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/CMakeLists.txt>`_
(and, for Spack builds, in ``scripts/spack_configs/versions.yaml`` of the same repository), which is the reference.
Some of these libraries may have their own system prerequisites.

Libraries
~~~~~~~~~

The following libraries are linked to by GEOS:

============= =========================== ============================= =====================================
Name          Enable option               Path variable                 Description
============= =========================== ============================= =====================================
Adiak_        :code:`ENABLE_CALIPER`      :code:`ADIAK_DIR`             Library for collecting metadata from HPC application runs, and distributing that metadata to subscriber tools.
Caliper_      :code:`ENABLE_CALIPER`      :code:`CALIPER_DIR`           Instrumentation and performance profiling library.
conduit_      *mandatory*                 :code:`CONDUIT_DIR`           Simplified Data Exchange for HPC Simulations.
CHAI_         *mandatory*                 :code:`CHAI_DIR`              Copy-hiding array abstraction to automatically migrate data between memory spaces (includes Umpire).
RAJA_         *mandatory*                 :code:`RAJA_DIR`              Collection of C++ software abstractions that enable architecture portability for HPC applications.
fmt_          *mandatory*                 :code:`FMT_DIR`               Formatting library.
hdf5_         *mandatory*                 :code:`HDF5_DIR`              High-performance data management and storage suite.
pugixml_      *mandatory*                 :code:`PUGIXML_DIR`           Light-weight, simple and fast XML parser for C++ with XPath support.
mathpresso_   :code:`ENABLE_MATHPRESSO`   :code:`MATHPRESSO_DIR`        Mathematical Expression Parser and JIT Compiler (used by ``SymbolicFunction``).
metis_        :code:`ENABLE_METIS`        :code:`METIS_DIR`             Serial graph partitioning library.
parmetis_     :code:`ENABLE_PARMETIS`     :code:`PARMETIS_DIR`          Parallel graph partitioning library. Should be built with 64-bit :code:`idx_t` type.
scotch_       :code:`ENABLE_SCOTCH`       :code:`SCOTCH_DIR`            Graph partitioning library (alternative to ParMETIS).
suitesparse_  :code:`ENABLE_SUITESPARSE`  :code:`SUITESPARSE_DIR`       A suite of sparse matrix software.
superlu_dist_ :code:`ENABLE_SUPERLU_DIST` :code:`SUPERLU_DIST_DIR`      General purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations.
hypre_        :code:`ENABLE_HYPRE`        :code:`HYPRE_DIR`             Library of high performance preconditioners and solvers for large sparse linear systems on massively parallel computers.
hypredrive_   :code:`ENABLE_HYPREDRV`     :code:`HYPREDRV_DIR`          Driver for configuring hypre solvers from YAML files.
PETSc_        :code:`ENABLE_PETSC`        :code:`PETSC_DIR`             Suite of data structures and routines for the scalable (parallel) solution of scientific applications.
Trilinos_     :code:`ENABLE_TRILINOS`     :code:`TRILINOS_DIR`          Collection of reusable scientific software libraries, known in particular for linear solvers.
silo_         :code:`ENABLE_SILO`         :code:`SILO_DIR`              A Mesh and Field I/O Library and Scientific Database.
VTK_          :code:`ENABLE_VTK`          :code:`VTK_DIR`               Open source software for manipulating and displaying scientific data (VTK mesh input and output).
============= =========================== ============================= =====================================

Tools
~~~~~~~~~

The following tools are used as part of the build process to support GEOS development:

============= =========================== ============================= =====================================
Name          Enable option               Path variable                 Description
============= =========================== ============================= =====================================
doxygen_      :code:`ENABLE_DOXYGEN`      :code:`DOXYGEN_EXECUTABLE`    De facto standard tool for generating documentation from annotated C++ sources.
sphinx_       :code:`ENABLE_SPHINX`       :code:`SPHINX_EXECUTABLE`     A tool that makes it easy to create intelligent and beautiful documentation.
uncrustify_   :code:`ENABLE_UNCRUSTIFY`   :code:`UNCRUSTIFY_EXECUTABLE` A source code beautifier for C, C++, C#, ObjectiveC, D, Java, Pawn and VALA.
============= =========================== ============================= =====================================

.. _Adiak : https://github.com/LLNL/Adiak
.. _Caliper: https://github.com/LLNL/Caliper
.. _conduit: https://github.com/LLNL/conduit
.. _CHAI : https://github.com/LLNL/CHAI
.. _RAJA : https://github.com/LLNL/RAJA
.. _hdf5 : https://www.hdfgroup.org/solutions/hdf5
.. _mathpresso : https://github.com/kobalicek/mathpresso
.. _pugixml : https://pugixml.org
.. _fmt : https://github.com/fmtlib/fmt
.. _metis : https://github.com/KarypisLab/METIS
.. _scotch : https://gitlab.inria.fr/scotch/scotch
.. _hypredrive : https://github.com/hypre-space/hypredrive
.. _parmetis : https://github.com/KarypisLab/ParMETIS
.. _silo : https://github.com/LLNL/Silo
.. _VTK : https://vtk.org/
.. _suitesparse : https://people.engr.tamu.edu/davis/suitesparse.html
.. _superlu_dist : https://portal.nersc.gov/project/sparse/superlu
.. _hypre : https://github.com/hypre-space/hypre
.. _PETSc : https://www.mcs.anl.gov/petsc
.. _Trilinos : https://trilinos.github.io
.. _doxygen : https://www.doxygen.nl/index.html
.. _sphinx : https://www.sphinx-doc.org/en/master/
.. _uncrustify : http://uncrustify.sourceforge.net
.. _GoogleTest : https://github.com/google/googletest
.. _GoogleBenchmark : https://github.com/google/benchmark
.. _BLT : https://github.com/LLNL/blt

Some other dependencies (GoogleTest_, GoogleBenchmark_) are provided through BLT_ build system which is embedded in GEOS source.
No actions are needed to build them.

Docker images with all dependencies already compiled are published for every TPL update (see :ref:`QuickStartDocker`).
The `Dockerfiles <https://github.com/GEOS-DEV/thirdPartyLibs/tree/master/docker>`_ used to create them are part of the thirdPartyLibs repository.

Building bundled dependencies
-----------------------------

To simplify the process of building TPLs, we provide a git repository `thirdPartyLibs <https://github.com/GEOS-DEV/thirdPartyLibs>`_.
It contains a CMake script that downloads the exact TPL versions required and builds all of them in a single command
(the build machine therefore needs internet access). A step-by-step example is given in the :ref:`QuickStart`.

The recommended steps to build TPLs are:

- Create a host-config file that sets all system-specific CMake variables (compiler and library paths, configuration flags, etc.)
  Take a look at `host-config examples <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs>`_.
- Configure via ``config-build.py`` script:

  .. code-block:: console

     cd thirdPartyLibs
     python3 scripts/config-build.py --hostconfig=/path/to/host-config.cmake --buildtype=Release --installpath=/path/to/install/dir -DNUM_PROC=8

  The options are the same as for GEOS (see :ref:`BuildProcess`). ``-DNUM_PROC=N`` sets the number of cores used to build
  each library (default: all cores of the machine).

- Run the build:

  .. code-block:: console

     cd <buildpath>
     make

  .. warning::
     Do not provide ``-j`` argument to ``make`` here, since the top-level make only launches sub-project builds.
     Instead use the ``-DNUM_PROC`` option above, which is passed to each sub-project's ``make`` command.

You may also run the CMake configure step manually instead of relying on ``config-build.py``.
The full TPL build may take anywhere between 30 minutes and a few hours, depending on your machine, number of threads and libraries enabled.

.. note::
   An exception from the above pattern, ``sphinx`` is currently not a part of the TPL bundle and must be installed with your Python or package manager.

Installing dependencies individually
------------------------------------

You may also install each individual TPL separately, either manually or through a package manager.
This is a more difficult route, since you are responsible for configuring dependencies in a compatible manner.
Again, we strongly recommend using the exact versions used by the thirdPartyLibs repository, to avoid possible build problems.

You may look at `our TPL CMake script <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/CMakeLists.txt>`_ to see how we configure TPL builds.
