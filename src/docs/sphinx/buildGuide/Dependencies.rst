.. _Dependencies:

Dependencies
============

GEOSX relies on multiple dependencies

The build management process it self uses

- `git <https://git-scm.com/>`_ and its large file extension plugin `git-lfs <https://git-lfs.github.com/>`_ (git 2.20 and git-lfs 2.7.1 are tested and working versions)
- `cmake <https://cmake.org/>`_ While GEOSX may surely be compiled with version 3.9, some of our third party libraries require a more recent version of cmake (cmake 3.13 is tested and working). We advise you to consider a rather recent version of cmake.
- C/C++/FORTRAN compilers (GEOSX uses c++14).
- `python <https://www.python.org/>`_ (version 2.7 tested and validated, only standard modules are used).

**It is the developer's responsibility to provide these tools.**
On Debian flavored distributions, consider the following command line.

.. code-block:: console

    apt-get install git git-lfs gcc g++ gfortran python2.7

The GEOSX software makes use of multiple libraries.
**Most of them are mirrored in the** `thidPartyLib/tplMirror <https://github.com/GEOSX/thirdPartyLibs/tree/master/tplMirror>`__ **folder or can be downloaded using the** `thirdPartyLibs/scripts/download_prerequisites.py <https://github.com/GEOSX/thirdPartyLibs/blob/master/scripts/download_prerequisites.py>`__ ** (see :ref:`QuickStart_download`). They will be configured and build by the GEOSX TPL building process.**
You may want to check the `CMakeLists.txt <https://github.com/GEOSX/thirdPartyLibs/blob/master/CMakeLists.txt>`_ that contains the versions of the dependencies.

LLNL HPC libraries and tool boxes...

- `axom <https://github.com/LLNL/axom>`_
- `Adiak <https://github.com/LLNL/Adiak>`_ (Library for collecting metadata from HPC application runs, and distributing that metadata to subscriber tools.)
- `caliper <https://github.com/LLNL/Caliper>`_ (Performance analysis toolbox.)
- `conduit <https://github.com/LLNL/conduit>`_ (Simplified data exchange tool for HPC simulations.)
- `CHAI <https://github.com/LLNL/CHAI>`_ (Library that handles automatic data migration to different memory spaces behind an array-style interface.)
- `RAJA <https://github.com/LLNL/RAJA>`_ (Collection of C++ software abstractions that enable architecture portability for HPC applications.)

... and more mainstream C++ stuff.

- `asmjit <https://github.com/asmjit/asmjit>`_ (Lightweight library for machine code just in time generation.)
- `camp <https://github.com/llnl/camp>`_ (Macros and metaprogramming facilities for C++.)
- `fparser <http://warp.povusers.org/FunctionParser>`_ (Parse and evaluate a mathematical function from a string.)
- `mathpresso <https://github.com/kobalicek/mathpresso>`_ (Mathematical expression parser and JIT compiler.)
- `pugixml <https://pugixml.org>`_ (Light-weight, simple and fast XML parser for C++ with XPath support.)

Linear (or more) algebra solvers

- `hypre <https://github.com/hypre-space/hypre>`_ (Library of high performance preconditioners and solvers featuring multigrid methods for the solution of large, sparse linear systems of equations on massively parallel computers.)
- `petsc <https://www.mcs.anl.gov/petsc>`_ (Suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations)
- `suitesparse <https://people.engr.tamu.edu/davis/suitesparse.html>`_ (A suite of sparse matrix software, offering different algorithms)
- `superlu_dist <https://portal.nersc.gov/project/sparse/superlu>`_ (General purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations).
- `trilinos <https://trilinos.github.io>`_ (Collection of reusable scientific software libraries, known in particular for linear solvers, non-linear solvers, transient solvers, optimization solvers, and uncertainty quantification (UQ) solvers.)

Note that petsc currently downloads `pt-scotch <https://www.labri.fr/perso/pelegrin/scotch/scotch_en.html>`_ from the internet.
If you do not have access to internet, you should modify the `./configure` step of petsc in the `CMakeLists.txt` file,
and change the ``--download-ptscotch`` option accordingly. 

A graph partitioning tool

- `parmetis <http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview>`_

A large and complex data collection library

- `hdf5 <https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git>`_ or its `clone on github <https://github.com/live-clones/hdf5>`_ and `its project page <https://portal.hdfgroup.org/display/knowledge>`_

Visualisation tools and libraries

- `silo <https://wci.llnl.gov/simulation/computer-codes/silo>`_
- `VTK <https://vtk.org/>`_

Code formatters & linters

- `astyle <http://astyle.sourceforge.net>`_
- `uncrustify <http://uncrustify.sourceforge.net>`_

Some of these libraries rely on blas, lapack & MPI implementations.
As well as the zlib compression library.
**It is the developer's responsibility to provide them.**

On Debian flavored distribution, consider installing (apt-get install) the following packages:

.. code-block:: console

    zlib1g-dev libblas-dev liblapack-dev libopenmpi-dev

Note also that `pt-scotch` relies on `bison` and `flex`.
The developer should provide these tools too.

GEOSX uses extensively the `BLT <https://github.com/LLNL/blt>`_ (a compilation of ``cmake`` functions for high performance computing) provided as a git submodule.
This brings also the `googletest <https://github.com/google/googletest>`_ framework which is used by GEOSX.

If you want further details about building the GEOSX and its TPLs, (which packages to install for example),
consider having a look at the `Dockerfiles in the third party library repository <https://github.com/GEOSX/thirdPartyLibs/tree/master/docker>`_.

.. _Third_party_libraries_build_management_pattern:

Third party libraries build management pattern
==============================================

Each dependency of GEOSX has its own build system.
The `thirdPartyLibs/CMakeLists.txt <https://github.com/GEOSX/thirdPartyLibs/blob/master/CMakeLists.txt>`_ cmake script was written to automate the build of all the external libraries.
In case the developer is no cmake expert (or simply to save time on a repetitive task),
it is convenient to use the `thirdPartyLibs/scripts/config-build.py <https://github.com/GEOSX/thirdPartyLibs/blob/master/scripts/config-build.py>`_ script that will manage the main key parameters of the dependency build.
This python script runs the proper cmake command line; generating the Makefiles accordingly.
Then running the ``make`` command (possibly with a specific target) will compile (and deploy) the libraries.

These two scripts are also used to build the docker images (see :ref:`Continuous_Integration_process`).

The most crucial parameters of the python script are ``--installpath``, ``--buildtype``, ``--hostconfig``.
(Other parameters do exist, check the script).

* ``--installpath`` is the installation directory. It wraps ``CMAKE_INSTALL_PREFIX``.
* ``--buildtype`` is a wrapper to the `CMAKE_BUILD_TYPE <https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html>`_ option.
* The ``--hostconfig`` option requires a cmake file containing some build parameters (compiler locations and flags, etc.).
  You may find some examples in the host-configs folders of the `third party library <https://github.com/GEOSX/thirdPartyLibs/tree/master/host-configs>`_ of from `GEOSX <https://github.com/GEOSX/GEOSX/tree/develop/host-configs>`_

To be more practical, you may need to run the following command line

.. code-block:: console

    python scripts/config-build.py --hostconfig=/path/to/your-platform.cmake --buildtype=Release --installpath=/path/to/install/dir

We do recommend using a *host config cmake file* for fine grained control of the build.
Have a look at some of the `already existing examples <https://github.com/GEOSX/GEOSX/blob/develop/host-configs>`_

Last, note that any extra argument will be transferred directly as a ``cmake`` argument.
For example, use the ``-DNUM_PROC=2`` to compile the TPL using two threads.

If you want to directly write the `cmake` command line, we advise you to dig into the `config-build.py <https://github.com/GEOSX/GEOSX/blob/develop/scripts/config-build.py>`_ python code.

Building GEOSX
==============

The same kind of `thirdPartyLibs/scripts/config-build.py <https://github.com/GEOSX/GEOSX/blob/develop/scripts/config-build.py>`_ (with the same main options) is used to build GEOSX.
In order to further customize your build, you can append any additional variable at the end of your command line.

Here is a non exhaustive list of options you may want to specify.

- ``-DNUM_PROC=4`` will allow you to compile with 4 parallel threads. (In GEOSX: to change this for the third party libraries, please modify in the code).
- ``-DGEOSX_TPL_DIR=/path/to/TPLs`` in case you did not use the default folder while building GEOSX and its third party libraries, you can use this options so GEOSX can find them.
- Some of the third party libraries can be activated/deactivated. Generally, the corresponding option looks like ``ENABLE_VTK``, ``ENABLE_CALIPER``...
- Computational features of GEOSX are activated with the following self-explanatory options: ``ENABLE_CUDA``, ``ENABLE_MPI``, ``ENABLE_OPENMP``.
- Building the documentation is controlled by the ``ENABLE_DOCS`` option.
- ``ENABLE_WARNINGS_AS_ERRORS``: GEOSX considers every warning as an error. When developing, you may face warnings however. You can modify this options (at your own risk) directly in the cmake scripts. Please understand that you won't be able to merge your code like this :)
