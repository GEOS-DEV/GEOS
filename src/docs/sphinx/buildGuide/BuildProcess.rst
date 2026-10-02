.. _BuildProcess:

Building GEOS
==============

Build steps
---------------------

- Create a host-config file that sets all system-specific CMake variables.
  Take a look at `host-config examples <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs>`_.
  We recommend the same host-config is used for both TPL and GEOS builds.
  In particular, certain options (such as ``ENABLE_MPI`` or ``ENABLE_CUDA``) need to match between the two.

- Provide paths to all enabled TPLs.
  This can be done in one of two ways:

  * Provide each path via a separate CMake variable (see :ref:`Dependencies` for path variable names).
  * If you built TPLs from the ``tplMirror`` repository, you can set ``GEOSX_TPL_DIR`` variable in your host-config to point to the TPL installation path, and

    .. code-block:: cmake

       include("/path/to/GEOS/host-configs/tpls.cmake")

    which will set all the individual TPL paths for you.

- Configure via ``config-build.py`` script:

  .. code-block:: console

     cd GEOS
     python scripts/config-build.py --hostconfig=/path/to/host-config.cmake --buildtype=Release --installpath=/path/to/install/dir

  where

  * ``--buildpath`` or ``-bp`` is the build directory (by default, created under current working dir).
  * ``--installpath`` or ``-ip`` is the installation directory(wraps ``CMAKE_INSTALL_PREFIX``).
  * ``--buildtype`` or ``-bt`` is a wrapper to the ``CMAKE_BUILD_TYPE`` option.
  * ``--hostconfig`` or ``-hc`` is a path to host-config file.
  * all unrecognized options are passed to CMake.

  If ``--buildpath`` is not used, build directory is automatically named ``build-<config-filename-without-extension>-<buildtype>``.
  It is possible to keep automatic naming and change the build root directory with ``--buildrootdir``.
  In that case, build path will be set to ``<buildrootdir>/<config-filename-without-extension>-<buildtype>``.
  Both ``--buildpath`` and ``--buildrootdir`` are incompatible and cannot be used in the same time.
  Same pattern is applicable to install path, with ``--installpath`` and ``--installrootdir`` options.

- Run the build:

  .. code-block:: console

     cd <buildpath>
     make -j $(nproc)

You may also run the CMake configure step manually instead of relying on ``config-build.py``.
A full build typically takes between 10 and 30 minutes, depending on chosen compilers, options and number of cores.

Configuration options
---------------------

Below is a list of CMake configuration options, in addition to TPL options above.
Some options, when enabled, require additional settings (e.g. ``ENABLE_CUDA``).
Please see `host-config examples <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs>`_.

=============================== ========= ==============================================================================
Option                          Default   Explanation
=============================== ========= ==============================================================================
``ENABLE_MPI``                  ``ON``    Build with MPI (also applies to TPLs)
``ENABLE_OPENMP``               ``OFF``   Build with OpenMP (also applies to TPLs)
``ENABLE_CUDA``                 ``OFF``   Build with CUDA (also applies to TPLs)
``ENABLE_CUDA_NVTOOLSEXT``      ``OFF``   Enable CUDA NVTX user instrumentation (via GEOS_MARK_SCOPE or GEOS_MARK_FUNCTION macros)
``ENABLE_HIP``                  ``OFF``   Build with HIP/ROCM (also applies to TPLs)
``ENABLE_DOCS``                 ``ON``    Build documentation (Sphinx and Doxygen)
``ENABLE_WARNINGS_AS_ERRORS``   ``ON``    Treat all warnings as errors
``ENABLE_TOTALVIEW_OUTPUT``     ``OFF``   Enables TotalView debugger custom view of GEOS data structures
``ENABLE_COV``                  ``OFF``   Enables code coverage
``GEOS_ENABLE_TESTS``           ``ON``    Enables unit testing targets
``GEOS_MAX_FLUID_COMPONENTS``   ``5``     Maximum component count instantiated in compositional flow solvers (2 to 20)
``GEOS_LA_INTERFACE``           ``Hypre`` Choiсe of Linear Algebra backend (Hypre/Petsc/Trilinos)
``GEOS_BUILD_OBJ_LIBS``         ``ON``    Use CMake Object Libraries build
``GEOS_BUILD_SHARED_LIBS``      ``OFF``   Build ``geosx_core`` as a shared library instead of static
``GEOS_PARALLEL_COMPILE_JOBS``            Max. number of compile jobs (when using Ninja), in addition to ``-j`` flag
``GEOS_PARALLEL_LINK_JOBS``               Max. number of link jobs (when using Ninja), in addition to ``-j`` flag
``GEOS_INSTALL_SCHEMA``         ``ON``    Enables schema generation and installation
=============================== ========= ==============================================================================

Compositional component limit
-----------------------------

Set ``GEOS_MAX_FLUID_COMPONENTS`` in the host configuration or pass, for example,
``-DGEOS_MAX_FLUID_COMPONENTS=20`` to CMake to build compositional flow solvers for up to
twenty components. The default of five preserves the existing solver range.
Values must be integers from two to twenty; smaller limits reduce compilation
work and the number of kernel instantiations. One-component kernels remain
available at every supported setting. The count includes every modeled fluid
species, including water when present.

This setting controls both runtime dispatch and explicit instantiations,
including thermal, hybrid finite-volume, aquifer, CFL and well kernels. Solvers
reject models exceeding the configured limit during setup,
with a message identifying the required setting. Changing the option requires
reconfiguring and rebuilding GEOS.

The constitutive fluid working arrays have capacity for at least nine components
and grow with ``GEOS_MAX_FLUID_COMPONENTS`` above nine. Smaller builds retain existing
nine-component standalone fluid-property calculations, such as the PVT driver.

The table-based reactive OBL solver uses the smaller of the configured limit and
seven components. Its per-thread interpolation workspace grows exponentially with
the component count and would exceed the CUDA local-memory limit at eight, so
increasing ``GEOS_MAX_FLUID_COMPONENTS`` above seven extends the EOS compositional
solvers without increasing the OBL limit.
