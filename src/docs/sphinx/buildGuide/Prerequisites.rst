.. _Prerequisites:

System prerequisites
====================

To configure and build GEOS you will need the following tools available on your system.

List of prerequisites
---------------------

Minimal requirements:

- `CMake <https://cmake.org/>`_ build system generator (3.24+).
- build tools (`GNU make <https://www.gnu.org/software/make/>`_ or `ninja <https://ninja-build.org/>`_ on Linux, XCode on MacOS).
- a C++ compiler with C++20 support. GEOS is tested with `gcc <https://gcc.gnu.org/>`_ 12 and 13, `clang <https://clang.llvm.org/>`_ 14 to 22, and Apple Clang on macOS.
- `python <https://www.python.org/>`_ 3.9-3.11 (versions 3.12+ are untested).
- :code:`zlib`, :code:`blas` and :code:`lapack` libraries
- any compatible MPI runtime and compilers (if building with MPI)

If you want to build from a repository check out (instead of a release tarball):

- `git <https://git-scm.com/>`_ (2.20+ is tested, but most versions should work fine)

If you plan on building bundled third-party library (TPLs) dependencies yourself:

- Compatible C and Fortran compilers
- internet access (the library sources are downloaded during the build)

If you will be running the integrated tests (see :ref:`IntegratedTests`):

- `h5py <https://www.h5py.org/>`_ and `mpi4py <https://pypi.org/project/mpi4py/>`_ python modules

If you are interested in building Doxygen documentation:

- `GNU bison <https://www.gnu.org/software/bison/>`_
- `LaTeX <https://www.latex-project.org/>`_
- `ghostscript <https://www.ghostscript.com/>`_
- `Graphviz <https://graphviz.org/>`_

In order for XML validation to work (executed as an optional build step):

- `xmllint <http://xmlsoft.org/xmllint.html>`_

Installing prerequisites
------------------------

On a local development machine with sudo/root privileges, most of these dependencies can be installed with a system package manager.
For example, on a Debian-based system such as Ubuntu 24.04 (check your package manager for specific package names):

.. code-block:: console

    sudo apt install build-essential git gcc g++ gfortran cmake libopenmpi-dev libblas-dev liblapack-dev zlib1g-dev python3 python3-h5py python3-mpi4py libxml2-utils

On HPC systems it is typical for these tools to be installed by system administrators and provided via `modules <http://modules.sourceforge.net/>`_.
To list available modules, type:

.. code-block:: console

    module avail

Then load the appropriate modules using :code:`module load` command.
Please contact your system administrator if you need help choosing or installing appropriate modules.

On macOS, see :ref:`AppleMacOS`.
