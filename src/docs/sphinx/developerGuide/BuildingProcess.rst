Building process
################

Dependencies
============

GEOSX relies on multiple dependencies

The build management process it self uses

- `git <https://git-scm.com/>`_ and its large file extension plugin `git-lfs <https://git-lfs.github.com/>`_ (git 2.20 and git-lfs 2.7.1 are tested and working versions)
- `cmake <https://cmake.org/>`_ version 3.9 (cmake 3.13 is tested and working).
- C/C++/FORTRAN compilers (GEOSX uses c++14).
- `python <https://www.python.org/>`_ (version 2.7 tested and validated, only standard modules are used).

**It is the developer’s responsibility to provide these tools.**
On Debian flavor distributions, consider the following command line.
.. code-block:: console
    apt-get install git git-lfs gcc g++ gfortran python2.7

The GEOSX software makes use of multiple libraries.
**Most of them are mirrored in the `thidPartyLib/tplMirror <https://github.com/GEOSX/thirdPartyLibs/tree/master/tplMirror>`_ folder and will be configured and build by the GEOSX TPL building process.**
You may want to check the `CMakeLists.txt <https://github.com/GEOSX/thirdPartyLibs/blob/master/CMakeLists.txt>`_ that contains the versions of the dependencies.

HPC LLNL's libraries and tool boxes...

- `axom <https://github.com/LLNL/axom>`_
- `caliper <https://github.com/LLNL/Caliper>`_
- `conduit <https://github.com/LLNL/conduit>`_
- `CHAI <https://github.com/LLNL/CHAI>`_
- `RAJA <https://github.com/LLNL/RAJA>`_

... and more mainstream C++ stuff.

- `asmjit <https://github.com/asmjit/asmjit>`_
- `fparser <http://warp.povusers.org/FunctionParser/>`_
- `pugixml <https://pugixml.org/>`_
- `mathpresso <https://github.com/kobalicek/mathpresso>`_

Linear (or more) algebra solvers

- `hypre <https://github.com/hypre-space/hypre>`_
- `petsc <https://www.mcs.anl.gov/petsc/>`_
- `superlu_dist <https://portal.nersc.gov/project/sparse/superlu/>`_
- `trilinos <https://trilinos.github.io/>`_

A graph partitioning tool

- `parmetis <http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview>`_

A large and comple data collection library

- `hdf5 <https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git>`_ or its `clone on github <https://github.com/live-clones/hdf5>`_ and `its project page <https://portal.hdfgroup.org/display/knowledge>`_

A LLNL's visualisation tool

- `silo <https://wci.llnl.gov/simulation/computer-codes/silo>`_

Code formatters & linters

- `astyle <http://astyle.sourceforge.net/>`_
- `uncrustify <http://uncrustify.sourceforge.net/>`_

Some of these libraries rely on blas, lapack & MPI implementations.
As well as the zlib compression library.
**It is the developer’s responsibility to provide them.**

On Debian flavored distribution, consider installing (apt-get install) the

.. code-block:: console
    zlib1g-dev libblas-dev liblapack-dev libopenmpi-dev

packages.

Third party libraries build management pattern
==============================================

Each dependency of GEOSX has its own build system.
The `thirdPartyLibs/CMakeLists.txt <https://github.com/GEOSX/thirdPartyLibs/blob/master/CMakeLists.txt>`_ cmake script was written to automatise the build of all the external libraries.
In case the developer is no cmake expert (or simply to save time on a repetitive task),
it is convenient to use the `thirdPartyLibs/scripts/config-build.py <https://github.com/GEOSX/thirdPartyLibs/blob/master/scripts/config-build.py>`_ script that will manage the main key parameters of the dependency build.
This python script runs the proper cmake command line; generating the Makefiles accordingly.
Then running the `make` command (possibly with a specific target) will compile (and deploy) the libraries.

These two scripts are also used the build the docker images easier without repeating ourselves (see `Continuous Integration process`_).

The most crucial parameters of the python script are `--installpath`, `--buildtype`, `--hostconfig`.
(Other parameters do exist, check the script).
While the first two parameters are obvious, the last one requires some explaination.
The `-—hostconfig` option requires a cmake file containing some build parameters (compiler location, etc.).
You may find some examples in the host-configs folders of the `third party library <https://github.com/GEOSX/thirdPartyLibs/tree/master/host-configs>`_ of from `GEOSX <https://github.com/GEOSX/GEOSX/tree/develop/host-configs>`_

Continuous Integration process
==============================

TODO

Troubleshooting
===============

TODO