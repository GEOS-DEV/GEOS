Continuous Integration process
==============================

To save building time, the third party libraries (that do not change so often) and GEOSX are build separately.

Everytime a pull is requested in the TPL repository, a docker image is generated and deployed on `dockerhub <https://hub.docker.com/r/geosx/compiler>`_.
The date (`YYYY-MM-DD`) is appended to the tag name so the client code (i.e. GEOSX) can select the version it needs
(the `DOCKER_DATE` env variable is defined in the `GEOSX's .travis.yml <https://github.com/GEOSX/GEOSX/blob/develop/.travis.yml>`_).

For the OSX builds, we build a tarball of the TPLs and save them a remote location.
The client (GEOSX again) will select the version it needs by defining the  `TPL_OSX_TRAVIS_BUILD_NUMBER` environment variable in the `.travis.yml <https://github.com/GEOSX/GEOSX/blob/develop/.travis.yml>`_ file.
An important counterpart to using a tarball and not a docker image is that the tarball does not provide the whole system the precompiled binaries rely on.
Problems may arise since we use the rolling release `Homebrew <https://brew.sh/>`_ (to install open-mpi in particular).
To circumvent this potential issue, the brew version is fixed to a specific commit (see BREW_HASH variable in `third party's .travis.yml <https://github.com/GEOSX/thirdPartyLibs/blob/master/.travis.yml>`_) and stored in a `brew_hash.txt` file at the root folder of the TPLs.
It is therefore possible for GEOSX to build against the same revision of brew packages.

It must be mentionned that one and only one version of the compiled TPL tarball is stored per pull request (older ones are removed automatically).
Therefore, a client building against a work in progress PR may experience a 404 error sooner or later.

It must be noted that there are now two different ways to designate the same version of the TPL.
An effort should be done to make this homogemneous.

Building docker images
----------------------

Our continuous integration process builds the TPL and GEOSX against two operating systems (ubuntu 18.04 and centos 7.5.1804) with respectively two compilers (clang versions 6 and 7, gcc on version 7 and 8).
Two dockerfiles are used, one per operating system, the compiler version being a parameter of the `docker build` stage.
The docker files use `multi-stage builds <https://docs.docker.com/develop/develop-images/multistage-build/>`_ in order to minimise the sizes of the images.

* First stage installs and defines all the elements that are commons to both TPL and GEOSX (for example, MPI and c++ compiler, BLAS, LAPACK, path to the installation directory...).
* As a second stage, we install everything needed to build (`not run`) the TPLs.
  We keep nothing from this second step for GEOSX, except the compiled TPL themselves.
  For example, a fortran compiler is needed by the TPL but not by GEOSX: it shall be installed during this step, so GEOSX won't access a fortran compiler (it does not have to).
* Last stage copies the compiled TPL from second stage and installs the elements only required by GEOSX (there are few).

Docker images `contract`
------------------------

GEOSX will find find a compiled version of the third party libraries.

As part of the contract provided by the TPL, the docker images also defines several environment variables.
The 

.. code-block::

    GEOSX_TPL_DIR

variable contains the absolute path of the installation root directory of the third party libraries.
GEOSX must use it when building.

Other variables are classical absolute path compiler variables.

.. code-block::

    CC
    CXX
    MPICC
    MPICXX

Anf the absolute path the mpirun (or equivalent) command.

.. code-block::

    MPIEXEC

The following openmpi environment variables allow it to work properly in the docker container.
But there should be no reason to access or use them explicitely.

.. code-block::

    OMPI_CC=$CC
    OMPI_CXX=$CXX

Using the TPL precompiled binaries
----------------------------------

For development purposes, you may want to use the publicly available docker images or the OSX tarball instead of compiling them yourself.
While this is surely possible, please note that *this is not supported by the GEOSX team that reserves the right to modify its workflow or delete elements on which you may have build your own workflow*.
