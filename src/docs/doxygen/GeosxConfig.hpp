/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GeosxConfig.hpp
 *
 * GEOSX build configuration file.
 * Contains a CMake-generated list of macros that define a particular build configuration.
 */

#ifndef GEOS_COMMON_CONFIG_HPP
#define GEOS_COMMON_CONFIG_HPP

/// Enables floating point exceptions
#define GEOS_USE_FPE

/// Enables bounds check in LvArray classes (CMake option GEOS_ENABLE_BOUNDS_CHECK)
/* #undef GEOS_USE_BOUNDS_CHECK */

/// Enables use of Caliper (CMake option ENABLE_CALIPER)
#define GEOS_USE_CALIPER

/// Enables use of Caliper (CMake option ENABLE_ADIAK)
/* #undef GEOS_USE_ADIAK */

/// Enables use of CHAI (CMake option ENABLE_CHAI)
#define GEOS_USE_CHAI

/// Enables use of Mathpresso library (CMake option ENABLE_MATHPRESSO)
#define GEOS_USE_MATHPRESSO

/// Enables use of MPI (CMake option ENABLE_MPI)
#define GEOS_USE_MPI

/// Enables detection of MPI desynchronization (CMake option ENABLE_MPI_DESYNC_DETECTION)
#define GEOS_USE_MPI_DESYNC_DETECTION

/// Enables use of OpenMP (CMake option ENABLE_OPENMP)
#define GEOS_USE_OPENMP

/// Enables use of CUDA (CMake option ENABLE_CUDA)
/* #undef GEOS_USE_CUDA */

/// Enables use of CUDA NVToolsExt (CMake option ENABLE_CUDA_NVTOOLSEXT)
/* #undef GEOS_USE_CUDA_NVTOOLSEXT */

/// Sets the CUDA stack size in KB (CMake option CUDA_STACK_SIZE)
/* #undef GEOS_CUDA_STACK_SIZE */

/// Enables use of HIP (CMake option ENABLE_HIP)
/* #undef GEOS_USE_HIP */

/// Workaround for FMT compilation issue on some NVCC/PowerPC machines (CMake option ENABLE_FMT_CONST_FORMATTER_WORKAROUND)
/* #undef GEOS_USE_FMT_CONST_FORMATTER_WORKAROUND */

/// Enables use of Python (CMake option ENABLE_PYTHON)
/* #undef GEOS_USE_PYGEOSX */

/// Enables use of RAJA (CMake option ENABLE_RAJA)
#define GEOS_USE_RAJA

/// Enables use of sys/time.h based timers (CMake option ENABLE_TIMERS)
/* #undef GEOS_USE_TIMERS */

/// Enables use of additional debugging interface for TotalView (Cmake option ENABLE_TOTALVIEW_OUTPUT)
/* #undef GEOS_USE_TOTALVIEW_OUTPUT */

/// Enables use of Intel MKL (CMake option ENABLE_MKL)
/* #undef GEOS_USE_MKL */

/// Enables use of Trilinos library (CMake option ENABLE_TRILINOS)
#define GEOS_USE_TRILINOS

/// Enables use of Hypre library (CMake option ENABLE_HYPRE)
#define GEOS_USE_HYPRE

/// Enables use of HYPREDRV library (CMake option ENABLE_HYPREDRV)
/* #undef GEOS_USE_HYPREDRV */

/// Denotes HYPRE using CPU
#define GEOS_USE_HYPRE_CPU 0
/// Denotes HYPRE using CUDA
#define GEOS_USE_HYPRE_CUDA 1
/// Denotes HYPRE using HIP
#define GEOS_USE_HYPRE_HIP 2
/// Macro determining what parellel interface hypre is using
#define GEOS_USE_HYPRE_DEVICE GEOS_USE_HYPRE_CPU

/// Enables use of SuperLU_dist library through HYPRE (CMake option ENABLE_SUPERLU_DIST)
#define GEOS_USE_SUPERLU_DIST

/// Enables use of PETSc library (CMake option ENABLE_PETSC)
/* #undef GEOS_USE_PETSC */

///  Enables use of METIS library (CMake option ENABLE_METIS)
#define GEOS_USE_METIS

///  Enables use of ParMETIS library (CMake option ENABLE_PARMETIS)
#define GEOS_USE_PARMETIS

/// Enables use of Scotch library (CMake option ENABLE_SCOTCH)
#define GEOS_USE_SCOTCH

/// Choice of global linear algebra interface (CMake option GEOS_LA_INTERFACE)
#define GEOS_LA_INTERFACE Hypre
/// Macro defined when Trilinos interface is selected
/* #undef GEOS_LA_INTERFACE_TRILINOS */
/// Macro defined when Hypre interface is selected
#define GEOS_LA_INTERFACE_HYPRE
/// Macro defined when PETSc interface is selected
/* #undef GEOS_LA_INTERFACE_PETSC */

/// Platform-dependent mangling of fortran function names (CMake option FORTRAN_MANGLE_NO_UNDERSCORE)
/* #undef FORTRAN_MANGLE_NO_UNDERSCORE */

/// USE OF SEPARATION COEFFICIENT IN FRACTURE FLOW
/* #undef GEOS_USE_SEPARATION_COEFFICIENT */

/// CMake option CMAKE_BUILD_TYPE
#define GEOS_CMAKE_BUILD_TYPE "Release"

/// The type that localIndex will be aliased to.
#define GEOS_LOCALINDEX_TYPE int

/// An integer flag representing the type that localIndex will be aliased to.
#define GEOS_LOCALINDEX_TYPE_FLAG 0

/// The type that globalIndex will be aliased to.
#define GEOS_GLOBALINDEX_TYPE long long int

/// An integer flag representing the type that globalIndex will be aliased to.
#define GEOS_GLOBALINDEX_TYPE_FLAG 2

/// The default block size for GEOSX on this platform
#define GEOS_BLOCK_SIZE 32

/// Version information for HDF5
#define HDF5_VERSION 1.12.1

/// Version information for Conduit
#define Conduit_VERSION 0.9.5

/// Version information for RAJA
#define RAJA_VERSION 2025.12.0

/// Version information for umpire
#define umpire_VERSION 2025.12.0

/// Version information for chai
#define chai_VERSION 2025.12.0

/// Version information for adiak
#define adiak_VERSION ..

/// Version information for caliper
#define caliper_VERSION ..

/// Version information for Metis
#define metis_VERSION 5.1.0

/// Version information for ParMetis
#define parmetis_VERSION 4.0.0

/// Version information for scotch
/* #undef scotch_VERSION */

/// Version information for superlu_dist
#define superlu_dist_VERSION 6.3.0

/// Version information for suitesparse
#define suitesparse_VERSION 5.7.9

/// Version information for hypre
#define hypre_VERSION 3.1.0

/// Version information for trilinos
#define trilinos_VERSION 16.1.0

/// Version information for petsc
/* #undef petsc_VERSION */

/// Version information for VTK
#define VTK_VERSION 9.4.2

/// Version information for fmt
/* #undef fmt_VERSION */

/// Version information for python
#define Python3_VERSION 3.9.13

/// Version information for CUDAToolkit
/* #undef CUDAToolkit_VERSION */

#if defined(GEOS_USE_CUDA) || defined(GEOS_USE_HIP)
// This needs to be placed into this header in order to appropriately replace
//  the old usage of GEOS_USE_CUDA, since we detect whether it is defined
//  rather than a value, not having it in the *same* header can cauase nebulous
//  compilation problems including the USD of arrays changing depending the scope
#define GEOS_USE_DEVICE
#endif

#endif  /* GEOS_CONFIG_HPP */
