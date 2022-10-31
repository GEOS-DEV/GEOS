/**
 * @file GeosxConfig.hpp
 *
 * GEOSX build configuration file.
 * Contains a CMake-generated list of macros that define a particular build configuration.
 */

#ifndef GEOSX_COMMON_CONFIG_HPP
#define GEOSX_COMMON_CONFIG_HPP

/// Enables floating point exceptions
#define GEOSX_USE_FPE

/// Enables bounds check in LvArray classes (CMake option ARRAY_BOUNDS_CHECK)
#define GEOSX_USE_ARRAY_BOUNDS_CHECK

/// Enables use of Caliper (CMake option ENABLE_CALIPER)
#define GEOSX_USE_CALIPER

/// Enables use of CHAI (CMake option ENABLE_CHAI)
#define GEOSX_USE_CHAI

/// Enables use of Mathpresso library (CMake option ENABLE_MATHPRESSO)
#define GEOSX_USE_MATHPRESSO

/// Enables use of MPI (CMake option ENABLE_MPI)
#define GEOSX_USE_MPI

/// Enables use of OpenMP (CMake option ENABLE_OPENMP)
#define GEOSX_USE_OPENMP

/// Enables use of CUDA (CMake option ENABLE_CUDA)
#define GEOSX_USE_CUDA

/// Enables use of CUDA NVToolsExt (CMake option ENABLE_CUDA_NVTOOLSEXT)
#define GEOSX_USE_CUDA_NVTOOLSEXT

/// Enables use of PVTPackage (CMake option ENABLE_PVTPackage)
#define GEOSX_USE_PVTPackage

/// Enables use of Python (CMake option ENABLE_PYTHON)
#define GEOSX_USE_PYGEOSX

/// Enables use of RAJA (CMake option ENABLE_RAJA)
#define GEOSX_USE_RAJA

/// Enables use of sys/time.h based timers (CMake option ENABLE_TIMERS)
#define GEOSX_USE_TIMERS

/// Enables use of additional debugging interface for TotalView (Cmake option ENABLE_TOTALVIEW_OUTPUT)
#define GEOSX_USE_TOTALVIEW_OUTPUT

/// Enables use of Intel MKL (CMake option ENABLE_MKL)
#define GEOSX_USE_MKL

/// Enables use of Trilinos library (CMake option ENABLE_TRILINOS)
#define GEOSX_USE_TRILINOS

/// Enables use of Hypre library (CMake option ENABLE_HYPRE)
#define GEOSX_USE_HYPRE

/// Macro defined when using cuda in HYPRE  (CMake option ENABLE_HYPRE_CUDA)
#define GEOSX_USE_HYPRE_CUDA

/// Enables use of PETSc library (CMake option ENABLE_PETSC)
#define GEOSX_USE_PETSC

/// Enables use of Scotch library (CMake option ENABLE_SCOTCH)
#define GEOSX_USE_SCOTCH

/// Choice of global linear algebra interface (CMake option GEOSX_LA_INTERFACE)
#define GEOSX_LA_INTERFACE Hypre
/// Macro defined when Trilinos interface is selected
/* #undef GEOSX_LA_INTERFACE_TRILINOS */
/// Macro defined when Hypre interface is selected
#define GEOSX_LA_INTERFACE_HYPRE
/// Macro defined when PETSc interface is selected
/* #undef GEOSX_LA_INTERFACE_PETSC */

/// Platform-dependent mangling of fortran function names (CMake option FORTRAN_MANGLE_NO_UNDERSCORE)
#define FORTRAN_MANGLE_NO_UNDERSCORE

/// USE OF SEPARATION COEFFICIENT IN FRACTURE FLOW
#define GEOSX_USE_SEPARATION_COEFFICIENT

/// CMake option CMAKE_BUILD_TYPE
#define GEOSX_CMAKE_BUILD_TYPE "Release"

/// The type that localIndex will be aliased to.
#define GEOSX_LOCALINDEX_TYPE std::ptrdiff_t

/// An integer flag representing the type that localIndex will be aliased to.
#define GEOSX_LOCALINDEX_TYPE_FLAG 3

/// The type that globalIndex will be aliased to.
#define GEOSX_GLOBALINDEX_TYPE long long int

/// An integer flag representing the type that globalIndex will be aliased to.
#define GEOSX_GLOBALINDEX_TYPE_FLAG 2

/// Version information for HDF5
#define HDF5_VERSION 1.12.1

/// Version information for Conduit
#define Conduit_VERSION 0.8.2

/// Version information for RAJA
#define RAJA_VERSION 2022.3.0

/// Version information for umpire
#define umpire_VERSION 2022.3.0

/// Version information for chai
/* #undef chai_VERSION */

/// Version information for adiak
#define adiak_VERSION ..

/// Version information for caliper
#define caliper_VERSION 2.8.0

/// Version information for Metis
#define METIS_VERSION 5.1.0

/// Version information for ParMetis
#define PARAMETIS_VERSION 4.0.3

/// Version information for scotch 
#define scotch_VERSION 6.0.9

/// Version information for superlu_dist
#define superlu_dist_VERSION 6.3.0

/// Version information for suitesparse
#define suitesparse_VERSION 5.7.9

/// Version information for VTK
#define VTK_VERSION 9.1.0

/// Version information for fmt
#define fmt_VERSION 8.0.1

/// Version information for python
/* #undef Python3_VERSION */

/// Version information for CUDAToolkit
/* #undef CUDAToolkit_VERSION */


#endif  /* GEOSX_CONFIG_HPP */

