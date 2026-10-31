.. _BuildProcess:

Building GEOS
==============

This page is a reference for the GEOS build. For a step-by-step walkthrough, see the :ref:`QuickStart`.

Build steps
---------------------

- Create a host-config file that sets all system-specific CMake variables.
  Start from ``host-configs/quick-start-template.cmake`` or one of the
  `host-config examples <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs>`_.
  Use the same host-config for both TPL and GEOS builds:
  certain options (such as ``ENABLE_MPI`` or ``ENABLE_CUDA``) need to match between the two.

- Provide paths to all enabled TPLs.
  This can be done in one of two ways:

  * If you built the TPLs with the `thirdPartyLibs <https://github.com/GEOS-DEV/thirdPartyLibs>`_ repository,
    set the ``GEOS_TPL_DIR`` variable in your host-config to the TPL installation path, and

    .. code-block:: cmake

       include("/path/to/GEOS/host-configs/tpls.cmake")

    which will set all the individual TPL paths for you.
  * Otherwise, provide each path via a separate CMake variable (see :ref:`Dependencies` for path variable names).

- Configure via the ``config-build.py`` script:

  .. code-block:: console

     cd GEOS
     python3 scripts/config-build.py --hostconfig=/path/to/host-config.cmake --buildtype=Release --installpath=/path/to/install/dir

  where

  * ``--hostconfig`` or ``-hc`` is a path to the host-config file.
  * ``--buildtype`` or ``-bt`` is a wrapper to the ``CMAKE_BUILD_TYPE`` option (``Release``, ``Debug``, ``RelWithDebInfo``).
  * ``--buildpath`` or ``-bp`` is the build directory (by default, created in the current working directory).
  * ``--installpath`` or ``-ip`` is the installation directory (wraps ``CMAKE_INSTALL_PREFIX``).
  * ``--ninja`` or ``-n`` uses the Ninja generator instead of Makefiles.
  * all unrecognized options (typically ``-D<OPTION>=<value>``) are passed to CMake.

  If ``--buildpath`` is not used, the build directory is named ``build-<config-filename-without-extension>-<buildtype>``.
  It is possible to keep automatic naming and change the build root directory with ``--buildrootdir``.
  ``--buildpath`` and ``--buildrootdir`` cannot be used at the same time.
  The same pattern applies to the install path, with ``--installpath`` and ``--installrootdir``.
  Run ``python3 scripts/config-build.py -h`` for the full list.

  .. warning::
     ``config-build.py`` deletes the build directory if it already exists.

- Run the build:

  .. code-block:: console

     cd <buildpath>
     make -j $(nproc)
     make install

You may also run the CMake configure step manually instead of relying on ``config-build.py``.
A full build typically takes between 10 and 30 minutes, depending on chosen compilers, options and number of cores.

Configuration options
---------------------

Below is a list of the most useful CMake configuration options, in addition to the TPL options listed in :ref:`Dependencies`.
Some options, when enabled, require additional settings (e.g. ``ENABLE_CUDA``).
Please see the `host-config examples <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs>`_.
The authoritative list, with default values, is in ``src/cmake/GeosOptions.cmake``.

General options
~~~~~~~~~~~~~~~

=============================== ========= ==============================================================================
Option                          Default   Explanation
=============================== ========= ==============================================================================
``ENABLE_MPI``                  ``ON``    Build with MPI (also applies to TPLs)
``ENABLE_OPENMP``               ``ON``    Build with OpenMP (also applies to TPLs). Defaults to ``OFF`` with Apple Clang.
``ENABLE_CUDA``                 ``OFF``   Build with CUDA (also applies to TPLs)
``ENABLE_CUDA_NVTOOLSEXT``      ``OFF``   Enable CUDA NVTX user instrumentation (via GEOS_MARK_SCOPE or GEOS_MARK_FUNCTION macros)
``ENABLE_HIP``                  ``OFF``   Build with HIP/ROCm (also applies to TPLs)
``ENABLE_HYPRE_DEVICE``         ``CPU``   Where hypre runs: ``CPU``, ``CUDA`` or ``HIP``
``GEOS_LA_INTERFACE``           ``Hypre`` Choice of linear algebra backend (``Hypre``, ``Trilinos`` or ``Petsc``)
``ENABLE_PYGEOSX``              ``OFF``   Build the ``pygeosx`` Python interface (see :ref:`pygeosxInterface`)
``ENABLE_DOCS``                 ``ON``    Build documentation (Sphinx and Doxygen)
``ENABLE_WARNINGS_AS_ERRORS``   ``ON``    Treat all warnings as errors (``OFF`` by default for ``Debug`` builds)
``ENABLE_COVERAGE``             ``OFF``   Enables code coverage
``ENABLE_TOTALVIEW_OUTPUT``     ``OFF``   Enables TotalView debugger custom view of GEOS data structures
``GEOS_ENABLE_TESTS``           ``ON``    Enables unit testing targets
``GEOS_ENABLE_FPE``             ``ON``    Enables floating point exceptions (traps on invalid operations, divisions by zero and overflows)
``GEOS_ENABLE_BOUNDS_CHECK``    ``OFF``   Enables array bounds checking (always ``ON`` for ``Debug`` builds)
``GEOS_BUILD_OBJ_LIBS``         ``OFF``   Build coreComponents modules as CMake object libraries
``GEOS_BUILD_SHARED_LIBS``      ``ON``    Build ``geosx_core`` as a shared library instead of static
``GEOS_MAX_FLUID_COMPONENTS``   ``5``     Maximum component count instantiated in compositional flow solvers (3 to 20)
``GEOS_LA_INTERFACE``           ``Hypre`` Choiсe of Linear Algebra backend (Hypre/Petsc/Trilinos)
``GEOS_BUILD_OBJ_LIBS``         ``ON``    Use CMake Object Libraries build
``GEOS_PARALLEL_COMPILE_JOBS``            Max. number of compile jobs (when using Ninja), in addition to ``-j`` flag
``GEOS_PARALLEL_LINK_JOBS``               Max. number of link jobs (when using Ninja), in addition to ``-j`` flag
``GEOS_INSTALL_SCHEMA``         ``ON``    Enables schema generation and installation
``GEOS_USE_GIT_VERSION_INFO``   ``ON``    Query git for version metadata printed at startup
=============================== ========= ==============================================================================

Physics packages
~~~~~~~~~~~~~~~~

Each family of physics solvers can be disabled to reduce build time and binary size.
All packages are enabled by default.

=================================== ====================================================================
Option                              Package
=================================== ====================================================================
``GEOS_ENABLE_FLUIDFLOW``           Single-phase, multiphase and compositional flow, wells
``GEOS_ENABLE_SOLIDMECHANICS``      Solid mechanics, including contact mechanics on fractures and faults
``GEOS_ENABLE_MULTIPHYSICS``        Coupled solvers (poromechanics, hydraulic fracturing, ...)
``GEOS_ENABLE_SURFACEGENERATION``   Fracture surface generation
``GEOS_ENABLE_WAVEPROPAGATION``     Acoustic and elastic wave propagation
``GEOS_ENABLE_INDUCEDSEISMICITY``   Induced seismicity
``GEOS_ENABLE_SIMPLEPDE``           Simple PDE solvers (Laplace, phase-field damage)
=================================== ====================================================================

.. note::
   The multiphysics package links against the fluid flow, solid mechanics, simple PDE and surface generation packages,
   so these must stay enabled when ``GEOS_ENABLE_MULTIPHYSICS`` is ``ON``.
