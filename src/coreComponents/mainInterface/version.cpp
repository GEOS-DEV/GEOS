/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "common/GeosxConfig.hpp"
#include "common/Logger.hpp"
#include "mainInterface/GeosxVersion.hpp"

namespace geos
{

std::string getVersion()
{
#if defined(GEOS_GIT_BRANCH) && defined(GEOS_GIT_HASH)
  return GEOS_VERSION_FULL " (" GEOS_GIT_BRANCH ", sha1: " GEOS_GIT_HASH ")";
#else
  return GEOS_VERSION_FULL;
#endif
}

static std::string getCppCompilerIdString()
{
  std::ostringstream oss;

  #if defined(__clang__)
  oss << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#if defined(__apple_build_version__)
  oss << " (apple version " << __apple_build_version__ << ")";
#endif
#if defined(__ibmxl_vrm__)
  oss << "IBMXL "
      << __ibmxl_version__ << "."
      << __ibmxl_release__ << "."
      << __ibmxl_modification__ << "."
      << __ibmxl_ptf_fix_level__;
#endif
#elif defined(__GNUC__)
  oss << "gcc "
      << __GNUC__ << "."
      << __GNUC_MINOR__ << "."
      << __GNUC_PATCHLEVEL__;
#endif
  return oss.str();
}

static std::string getGpuCompilerIdString()
{
  std::ostringstream oss;

#if defined( GEOS_USE_CUDA )
  oss << "  - CUDA compiler version: " << CUDA_VERSION/10/100 << "." << CUDA_VERSION/10%100;
#endif
#if defined( GEOS_USE_HIP )
  oss << "  - ROCm compiler version: " << ROCM_VERSION/100/100 << "." << ROCM_VERSION/100%100;
#endif
  return oss.str();
}

void outputVersionInfo()
{
  // TODO use a table layout here

  logger.rank0log( "GEOS version: ", getVersion() );

  logger.rank0log( "  - c++ compiler: ", getCppCompilerIdString() );

  std::string const gpuCompilerIdString = getGpuCompilerIdString();
  if (!gpuCompilerIdString.empty())
    logger.rank0log( "  - gpu compiler: ", gpuCompilerIdString );

#if defined(_OPENMP)
  logger.rank0log( "  - openmp version: ", _OPENMP );
#endif

#if defined(GEOS_USE_MPI)
  {
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int len;
    MPI_Get_library_version( version, &len );
    logger.rank0log( "  - MPI version: ", version );
  }
#endif

#if defined(HDF5_VERSION)
  logger.rank0log( "  - HDF5 version: ", STRINGIZE( HDF5_VERSION ) );
#endif

#if defined(Conduit_VERSION)
  logger.rank0log( "  - Conduit version: ", STRINGIZE( Conduit_VERSION ) );
#endif

#if defined(VTK_VERSION)
  logger.rank0log( "  - VTK version: ", STRINGIZE( VTK_VERSION ) );
#endif

#if defined(RAJA_VERSION)
  logger.rank0log( "  - RAJA version: ", STRINGIZE( RAJA_VERSION ) );
#endif

#if defined(umpire_VERSION)
  logger.rank0log( "  - umpire version: ", STRINGIZE( umpire_VERSION ) );
#endif

#if defined(chai_VERSION)
  logger.rank0log( "  - chai version: ", STRINGIZE( chai_VERSION ) );
#endif

#if defined(adiak_VERSION)
  logger.rank0log( "  - adiak version: ", STRINGIZE( adiak_VERSION ) );
#endif

#if defined(caliper_VERSION)
  logger.rank0log( "  - caliper version: ", STRINGIZE( caliper_VERSION ) );
#endif

#if defined(metis_VERSION)
  logger.rank0log( "  - METIS version: ", STRINGIZE( metis_VERSION ) );
#endif

#if defined(parmetis_VERSION)
  logger.rank0log( "  - PARAMETIS version: ", STRINGIZE( parmetis_VERSION ) );
#endif

#if defined(scotch_VERSION)
  logger.rank0log( "  - scotch version: ", STRINGIZE( scotch_VERSION ) );
#endif

#if defined(superlu_dist_VERSION)
  logger.rank0log( "  - superlu_dist version: ", STRINGIZE( superlu_dist_VERSION ) );
#endif

#if defined(suitesparse_VERSION)
  logger.rank0log( "  - suitesparse version: ", STRINGIZE( suitesparse_VERSION ) );
#endif

#if defined(hypre_VERSION)
  logger.rank0log( "  - hypre version: ", STRINGIZE( hypre_VERSION ) );
#endif

#if defined(trilinos_VERSION)
  logger.rank0log( "  - trilinos version: ", STRINGIZE( trilinos_VERSION ) );
#endif

#if defined(petsc_VERSION)
  logger.rank0log( "  - petsc version: ", STRINGIZE( petsc_VERSION ) );
#endif

#if defined(Python3_VERSION)
  logger.rank0log( "  - Python3 version: ", STRINGIZE( Python3_VERSION ) );
#endif

#if defined(CUDAToolkit_VERSION)
  logger.rank0log( "  - CUDAToolkit version: ", STRINGIZE( CUDAToolkit_VERSION ) );
#endif

#if \
  defined(GEOS_USE_DEVICE) && \
  defined(GEOS_USE_HYPRE) && \
  ( GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CPU )
  logger.rank0log( "" );
  logger.rank0log( "**************************************************" );
  logger.rank0log( "*                   WARNING!!!                   *" );
  logger.rank0log( "*                                                *" );
  logger.rank0log( "*  GEOS has GPU support enabled, but not HYPRE!  *" );
  logger.rank0log( "**************************************************" );
  logger.rank0log( "" );
#endif

}

}
