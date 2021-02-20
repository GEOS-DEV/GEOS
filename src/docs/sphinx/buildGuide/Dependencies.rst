.. _Dependencies:

Third-party library dependencies
================================

GEOSX makes use of multiple third-party libraries (TPLs) and tools, some of which are mandatory and some optional.
We only test against specific versions, and sometimes even require development snapshots (specific git commits).
Not all of these guarantee backwards compatibility, so we strongly recommend building with these specific versions.

List of third-party libraries
-----------------------------

The table below lists the dependencies with their specific versions and relevant CMake variables.
Some of these libraries may have their own system prerequisites.
For example, Doxygen depends on ``flex``/``bison``, ``latex``, ``ghostscript`` and ``graphviz``.

============= ========== =========================== ============================= ========================================
Name          Version    Enable option               Path variable                 Other options and comments
============= ========== =========================== ============================= ========================================
axom_         0.3.1      *mandatory*                 :code:`AXOM_DIR`
Adiak_        0.2.0      :code:`ENABLE_CALIPER`      :code:`ADIAK_DIR`
Caliper_      2.4.0      :code:`ENABLE_CALIPER`      :code:`CALIPER_DIR`
conduit_      0.5.0      *mandatory*                 :code:`CONDUIT_DIR`
CHAI_         2.2.2      *mandatory*                 :code:`CHAI_DIR`
RAJA_         0.12.1     *mandatory*                 :code:`RAJA_DIR`
hdf5_         1.10.5     *mandatory*                 :code:`HDF5_DIR`
mathpresso_   2015-12-15 :code:`ENABLE_MATHPRESSO`   :code:`MATHPRESSO_DIR`
pugixml_      1.8.0      *mandatory*                 :code:`PUGIXML_DIR`
parmetis_     4.0.3      *mandatory* (with MPI)      :code:`PARMETIS_DIR`          built with 64-bit :code:`idx_t`
suitesparse_  5.8.1      :code:`ENABLE_SUITESPARSE`  :code:`SUITESPARSE_DIR`
superlu_dist_ 0f6efc3    :code:`ENABLE_SUPERLU_DIST` :code:`SUPERLU_DIST_DIR`
hypre_        2186a8f    :code:`ENABLE_HYPRE`        :code:`HYPRE_DIR`
PETSc_        3.13.0     :code:`ENABLE_PETSC`        :code:`PETSC_DIR`
Trilinos_     12.18.1    :code:`ENABLE_TRILINOS`     :code:`TRILINOS_DIR`
silo_         4.10.3     *mandatory*                 :code:`SILO_DIR`
VTK_          9.0.0-rc3  :code:`ENABLE_VTK`          :code:`VTK_DIR`               only a small subset of modules required
doxygen_      1.8.20     :code:`ENABLE_DOXYGEN`      :code:`DOXYGEN_EXECUTABLE`
sphinx_       1.8.20     :code:`ENABLE_SPHINX`       :code:`SPHINX_EXECUTABLE`
uncrustify_   401a409    :code:`ENABLE_UNCRUSTIFY`   :code:`UNCRUSTIFY_EXECUTABLE`
============= ========== =========================== ============================= ========================================

.. _axom : https://github.com/LLNL/axom
.. _Adiak : https://github.com/LLNL/Adiak
.. _Caliper: https://github.com/LLNL/Caliper
.. _conduit: https://github.com/LLNL/conduit
.. _CHAI : https://github.com/LLNL/CHAI
.. _RAJA : https://github.com/LLNL/RAJA
.. _hdf5 : https://portal.hdfgroup.org/display/HDF5/HDF5
.. _mathpresso : https://github.com/kobalicek/mathpresso
.. _pugixml : https://pugixml.org
.. _parmetis : http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
.. _silo : https://wci.llnl.gov/simulation/computer-codes/silo
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

Some other dependencies (GoogleTest_, GoogleBenchmark_) are provided through BLT_ build system which is embedded in GEOSX source.
No actions are needed to build them.

If you would like to create a Docker image with all dependencies, take a look at
`Dockerfiles <https://github.com/GEOSX/thirdPartyLibs/tree/master/docker>`_
that are used in our CI process.

Building bundled dependencies
-----------------------------

To simplify the process of building TPLs, we provide a git repository `thirdPartyLibs <https://github.com/GEOSX/thirdPartyLibs>`_.
It contains source copies of exact TPL versions required and is updated periodically.
It also contains a CMake script for building all TPLs in a single command.

The recommended steps to build TPLs are:

- Create a host-config file that sets all system-specific CMake variables (compiler and library paths, configuration flags, etc.)
  Take a look at `host-config examples <https://github.com/GEOSX/GEOSX/blob/develop/host-configs>`_.
- Configure via ``config-build.py`` script:

  .. code-block:: console

     cd thirdPartyLibs
     python scripts/config-build.py --hostconfig=/path/to/host-config.cmake --buildtype=Release --installpath=/path/to/install/dir -DNUM_PROC=8

  where

  * ``--buildpath`` or ``-bp`` is the build directory (by default, created under current).
  * ``--installpath`` or ``-ip`` is the installation directory(wraps ``CMAKE_INSTALL_PREFIX``).
  * ``--buildtype`` or ``-bt`` is a wrapper to the ``CMAKE_BUILD_TYPE`` option.
  * ``--hostconfig`` or ``-hc`` is a path to host-config file.
  * all unrecognized options are passed to CMake.

- Run the build:

  .. code-block:: console

     cd <buildpath>
     make

  .. warning::
     Do not provide ``-j`` argument to ``make`` here, since the top-level make only launches sub-project builds.
     Instead use ``-DNUM_PROC`` option above, which is passed to each sub-project's ``make`` command.

You may also run the CMake configure step manually instead of relying on ``config-build.py``.
The full TPL build may take anywhere between 15 minutes and 2 hours, depending on your machine, number of threads and libraries enabled.

.. note::
   An exception from the above pattern, ``sphinx`` is currently not a part of the TPL bundle and must be installed with your Python or package manager.

.. note::
   PETSc build currently downloads `pt-scotch <https://www.labri.fr/perso/pelegrin/scotch/scotch_en.html>`_ from the internet.
   If you do not have access to internet, modify the `./configure` step of petsc in `CMakeLists.txt` and change the ``--download-ptscotch`` option accordingly.
   `pt-scotch` also relies on `bison` and `flex`.

Installing dependencies individually
------------------------------------

You may also install each individual TPL separately, either manually or through a package manager.
This is a more difficult route, since you are responsible for configuring dependencies in a compatible manner.
Again, we strongly recommend using the exact versions listed above, to avoid possible build problems.

You may look at `our TPL CMake script <https://github.com/GEOSX/thirdPartyLibs/blob/master/CMakeLists.txt>`_ to see how we configure TPL builds.
