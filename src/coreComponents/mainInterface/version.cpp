/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "version.hpp"
#include "mainInterface/GeosxVersion.hpp"

namespace geos
{

string getVersion()
{
#if defined(GEOSX_GIT_BRANCH) && defined(GEOSX_GIT_HASH)
  return GEOSX_VERSION_FULL " (" GEOSX_GIT_BRANCH ", sha1: " GEOSX_GIT_HASH ")";
#else
  return GEOSX_VERSION_FULL;
#endif
}

string getCppCompilerIdString()
{
  std::ostringstream oss;

  #if defined(__clang__)
  oss<<"clang " <<__clang_major__<<"."<<__clang_minor__<<"."<<__clang_patchlevel__;
#if defined(__apple_build_version__)
  oss<<" ( apple version " <<__apple_build_version__<<" )";
#endif
#if defined(__ibmxl_vrm__)
  oss<<"IBMXL "
     <<__ibmxl_version__<<"."
     <<__ibmxl_release__<<"."
     <<__ibmxl_modification__<<"."
     <<__ibmxl_ptf_fix_level__;
#endif
#elif defined(__GNUC__)
  oss<<"gcc "
     <<__GNUC__<<"."
     <<__GNUC_MINOR__<<"."
     <<__GNUC_PATCHLEVEL__;
#endif
  return oss.str();
}

string getGpuCompilerIdString()
{
  std::ostringstream oss;

#if defined( GEOS_USE_CUDA )
  oss<<"  - cuda compiler version: " <<CUDA_VERSION/1000<<"."<<CUDA_VERSION/10%100;
#endif
  return oss.str();
}

void outputVersionInfo()
{

  logger.rank0Log( "GEOSX version: ", getVersion() );

  logger.rank0Log( "  - c++ compiler: ", getCppCompilerIdString() );

  string const gpuCompilerIdString = getGpuCompilerIdString();
  GEOS_LOG_RANK_0_IF( !gpuCompilerIdString.empty(), gpuCompilerIdString );


#if defined(_OPENMP)
  logger.rank0Log( "  - openmp version: ", _OPENMP );
#endif

#if defined(GEOSX_USE_MPI)
  {
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int len;
    MPI_Get_library_version( version, &len );
    logger.rank0Log( "  - MPI version: ", version );
  }
#endif

#if defined(HDF5_VERSION)
  logger.rank0Log( "  - HDF5 version: ", STRINGIZE( HDF5_VERSION ) );
#endif

#if defined(Conduit_VERSION)
  logger.rank0Log( "  - Conduit version: ", STRINGIZE( Conduit_VERSION ) );
#endif

#if defined(VTK_VERSION)
  logger.rank0Log( "  - VTK version: ", STRINGIZE( VTK_VERSION ) );
#endif

#if defined(RAJA_VERSION)
  logger.rank0Log( "  - RAJA version: ", STRINGIZE( RAJA_VERSION ) );
#endif
#if defined(umpire_VERSION)
  logger.rank0Log( "  - umpire version: ", STRINGIZE( umpire_VERSION ) );
#endif
#if defined(chai_VERSION)
  logger.rank0Log( "  - chai version: ", STRINGIZE( chai_VERSION ) );
#endif

#if defined(adiak_VERSION)
  logger.rank0Log( "  -  adiak version: ", STRINGIZE( adiak_VERSION ) );
#endif

#if defined(caliper_VERSION)
  logger.rank0Log( "  - caliper version: ", STRINGIZE( caliper_VERSION ) );
#endif

#if defined(METIS_VERSION)
  logger.rank0Log( "  - METIS version: ", STRINGIZE( METIS_VERSION ) );
#endif

#if defined(PARAMETIS_VERSION)
  logger.rank0Log( "  - PARAMETIS version: ", STRINGIZE( PARAMETIS_VERSION ) );
#endif

#if defined(scotch_VERSION)
  logger.rank0Log( "  - scotch version: ", STRINGIZE( scotch_VERSION ) );
#endif

#if defined(superlu_dist_VERSION)
  logger.rank0Log( "  - superlu_dist version: ", STRINGIZE( superlu_dist_VERSION ) );
#endif

#if defined(suitesparse_VERSION)
  logger.rank0Log( "  - suitesparse version: ", STRINGIZE( suitesparse_VERSION ) );
#endif

#if defined(Python3_VERSION)
  logger.rank0Log( "  - Python3 version: ", STRINGIZE( Python3_VERSION ) );
#endif

#if defined(CUDAToolkit_VERSION)
  logger.rank0Log( "  - CUDAToolkit version: ", STRINGIZE( CUDAToolkit_VERSION ) );
#endif

}

}
