##################################################################################
Getting Started with the GEOSX Solid Mechanics Lagrangian Finite Element Proxy-App
##################################################################################

Summary of Code
=================================
The Solid Mechanics Lagrangian Finite Element Proxy-App is designed to be a framework to explore different data layouts
for the most compute intensive kernel within the GEOSX's physics solver package.
The code includes a single makefile that will build the code in a variety of platforms.
Using macros and pre-processors, the code provides a simple way to alter how data is laid out
in large arrays. The rest of the document will outline how to build
the proxy. A more detailed description of the code is found in :ref:`about-proxy`


Getting Ready
=================================
The proxy is found inside the PhysicsPackage1 folder::

    /PhysicsSolverPackage1/src/SolidMechanicsLagrangianFEM-MiniApp

In particular, ``feature/artv3/mini-app-shapefun`` branch in the PhysicsPackage1 repository holds the
latest version of the mini-app. For instructions on cloning a copy from the Github repository please see::

https://github.com/GEOSX/GEOSX/blob/develop/src/docs/sphinx/getting_started.rst

Compiling and Running Code
=================================
A single makefile simplifies building the mini-app but does require users to provide paths
to the RAJA and CUDA directories::

  RAJA_DIR ?= /g/g17/vargas45/Git-Repo/RAJA/develop/build
  CUDA_DIR ?= /usr/tce/packages/cuda/cuda-9.2.64

The RAJA directory points to the installed directory of RAJA, for instructions on installing RAJA please
see the RAJA Github page::

  https://github.com/llnl/raja

If CUDA is available, the location of the CUDA directory may be found by invoking
the ``which`` command::

  $ which nvcc
  /usr/tce/packages/cuda/cuda-9.2.88/bin/nvcc

In this case the location of the CUDA directory is ::

  /usr/tce/packages/cuda/cuda-9.2.88

The included makefile will compile and build the code on the following machines and compilers:

**Quartz**

* ``GNU-7``
* ``Clang-4.0.0``

**Ray**

* ``xl-beta-2018-06-01``
* ``Clang-coral-2018-05-23``

The default the GPU compiler is assumed to be nvcc and may be paired with the Clang-coral
compiler as the host compiler.

**CORI**

* ``Intel`` - default version

The command ``make`` will build the code and ``make clean`` will remove executables.
Swapping between the compilers and compiler flags within a machine is accomplished by
commenting in/out lines the makefile. Post the build processes two executables will be generated

* ``ArrayOfObjects``

* ``ObjectOfArrays``

The ArrayOfObjects executable will store nodal degrees of freedom and cell centered values
in an array of objects.
The ObjectOfArrays executables will store nodal degrees of freedom and cell centered values as an object of arrays.
Lastly, running the code is accomplished in the following manner::

  ./ArrayOfObjects NoElem Niter
  ./ObjectOfArrays NoElem Niter

where the arguments are

* NoElem - The number of elements in a cartesian dimension

* Niter - The number of times the kernel will be executed.

Post the computation the drivers will print out on the terminal the minimum kernel runtime,
whether the computation was done on the CPU/GPU, and total number of elements
used in the computation. 
