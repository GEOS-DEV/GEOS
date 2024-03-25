/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-LiCense-Identifier: LGPL-2.1-only
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

#include "initializeEnvironment.hpp"

#include "TimingMacros.hpp"
#include "Path.hpp"
#include "LvArray/src/system.hpp"
#include "codingUtilities/TableLayout.hpp"
#include "codingUtilities/TableData.hpp"
#include "codingUtilities/TableFormatter.hpp"
// TPL includes
#include <umpire/ResourceManager.hpp>

#if defined( GEOSX_USE_CALIPER )
#include <caliper/cali-manager.h>
#if defined( GEOSX_USE_ADIAK )
#include <adiak.hpp>
#endif
#endif

// System includes
#include <iomanip>

#if defined( GEOSX_USE_MKL )
#include <mkl.h>
#endif

#if defined( GEOSX_USE_OPENMP )
#include <omp.h>
#endif

#if defined( GEOS_USE_CUDA )
#include <cuda.h>
#endif

#if defined( GEOS_USE_HIP )
#include <hip/hip_runtime.h>
#endif
#include <cfenv>

namespace geos
{


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupLogger()
{
#ifdef GEOSX_USE_MPI
  logger::InitializeLogger( MPI_COMM_GEOSX );
#else
  logger::InitializeLogger();
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeLogger()
{
  logger::FinalizeLogger();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupLvArray()
{
  LvArray::system::setErrorHandler( []()
  {
  #if defined( GEOSX_USE_MPI )
    int mpi = 0;
    MPI_Initialized( &mpi );
    if( mpi )
    {
      MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
    }
  #endif
    std::abort();
  } );

  LvArray::system::setSignalHandling( []( int const signal ) { LvArray::system::stackTraceHandler( signal, true ); } );

#if defined(GEOSX_USE_FPE)
  LvArray::system::setFPE();
#else
  LvArray::system::disableFloatingPointExceptions( FE_ALL_EXCEPT );
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupMKL()
{
#ifdef GEOSX_USE_MKL
  GEOS_LOG_RANK_0( "MKL max threads: " << mkl_get_max_threads() );
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupOpenMP()
{
#ifdef GEOSX_USE_OPENMP
  GEOS_LOG_RANK_0( "Max threads: " << omp_get_max_threads() );
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupMPI( int argc, char * argv[] )
{
  if( !MpiWrapper::initialized() )
  {
    MpiWrapper::init( &argc, &argv );
  }

  MPI_COMM_GEOSX = MpiWrapper::commDup( MPI_COMM_WORLD );

  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
  {
    // Can't use logging macros prior to logger init
    std::cout << "Num ranks: " << MpiWrapper::commSize( MPI_COMM_GEOSX ) << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeMPI()
{
  MpiWrapper::commFree( MPI_COMM_GEOSX );
  MpiWrapper::finalize();
}

#if defined( GEOSX_USE_CALIPER )

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupCaliper( cali::ConfigManager & caliperManager,
                   CommandLineOptions const & commandLineOptions )
{
  caliperManager.add( commandLineOptions.timerOutput.c_str() );
  GEOS_ERROR_IF( caliperManager.error(), "Caliper config error: " << caliperManager.error_msg() );
  caliperManager.start();

#if defined( GEOSX_USE_ADIAK )
#if defined( GEOSX_USE_MPI )
  adiak::init( &MPI_COMM_GEOSX );
#else
  adiak::init( nullptr );
#endif

  GEOS_WARNING_IF( !adiak::uid(), "Error getting the user info." );
  GEOS_WARNING_IF( !adiak::launchdate(), "Error getting the launch date info." );
  GEOS_WARNING_IF( !adiak::cmdline(), "Error getting the command line args." );
  GEOS_WARNING_IF( !adiak::clustername(), "Error getting the clustername." );
  GEOS_WARNING_IF( !adiak::walltime(), "Error getting the walltime." );
  GEOS_WARNING_IF( !adiak::systime(), "Error getting the systime." );
  GEOS_WARNING_IF( !adiak::cputime(), "Error getting the cputime." );

  for( auto & fileName: commandLineOptions.inputFileNames )
  {
    adiak::value( "XML File", splitPath( fileName ).second );
  }
  adiak::value( "Problem name", commandLineOptions.problemName );

  // MPI info
#if defined( GEOSX_USE_MPI )
  adiak::value( "MPI", "On" );
  adiak::value( "mpi ranks", MpiWrapper::commSize() );
#else
  adiak::value( "MPI", "Off" );
  adiak::value( "mpi ranks", 1 );
#endif

  // Build info
#if defined( __clang_version__ )
  adiak::value( "compiler", "clang" );
  adiak::value( "compiler version", adiak::version( "clang" __clang_version__ ) );
#elif defined( __INTEL_COMPILER )
  adiak::value( "compiler", "intel" );
  adiak::value( "compiler version", adiak::version( "intel" __INTEL_COMPILER ) );
#elif defined( __GNUC__ )
  adiak::value( "compiler", "gcc" );
  adiak::value( "compiler version", adiak::version( "gcc" __VERSION__ ) );
#else
  adiak::value( "compiler", "unknown" );
  adiak::value ( "compiler version", "unknown" );
#endif

  adiak::value( "build type", GEOSX_CMAKE_BUILD_TYPE );
  adiak::value( "compilation date", __DATE__ );

  // OpenMP info
#if defined( GEOSX_USE_OPENMP )
  std::int64_t const numThreads = omp_get_max_threads();
  adiak::value( "OpenMP", "On" );
#else
  std::int64_t const numThreads = 1;
  adiak::value( "OpenMP", "Off" );
#endif
  pushStatsIntoAdiak( "numThreads", numThreads );

  // CUDA info
  int cudaRuntimeVersion = 0;
  int cudaDriverVersion = 0;
#if defined( GEOS_USE_CUDA )
  adiak::value( "CUDA", "On" );
  GEOS_ERROR_IF_NE( cudaSuccess, cudaRuntimeGetVersion( &cudaRuntimeVersion ) );
  GEOS_ERROR_IF_NE( cudaSuccess, cudaDriverGetVersion( &cudaDriverVersion ) );
#else
  adiak::value( "CUDA", "Off" );
#endif
  adiak::value( "CUDA runtime version", cudaRuntimeVersion );
  adiak::value( "CUDA driver version", cudaDriverVersion );

  // HIP info
  int hipRuntimeVersion = 0;
  int hipDriverVersion = 0;
#if defined( GESOX_USE_HIP )
  adiak::value( "HIP", "On" )
  GEOSX_ERROR_IF_NE( hipSuccess, hipRuntimeGetVersion( &hipRuntimeVersion ) );
  GEOSX_ERROR_IF_NE( hipSuccess, hipDriverGetVersion( &hipDriverVersion ) );
#else
  adiak::value( "HIP", "Off" );
#endif
  adiak::value( "HIP runtime version", hipRuntimeVersion );
  adiak::value( "HIP driver version", hipDriverVersion );
#endif // defined( GEOSX_USE ADIAK )
}
#endif // defined( GEOSX_USE_CALIPER )

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeCaliper()
{
#if defined( GEOSX_USE_CALIPER )and defined( GEOSX_USE_ADIAK )
  adiak::fini();
#endif
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief For each Umpire::Allocator compute the total high water mark across all ranks
 *        and if using Adiak add statistics about the high water mark.
 */
static void addUmpireHighWaterMarks()
{
  umpire::ResourceManager & rm = umpire::ResourceManager::getInstance();

  // Get a list of all the allocators and sort it so that it's in the same order on each rank.
  std::vector< string > allocatorNames = rm.getAllocatorNames();
  std::sort( allocatorNames.begin(), allocatorNames.end() );

  // If each rank doesn't have the same number of allocators you can't aggregate them.
  std::size_t const numAllocators = allocatorNames.size();
  std::size_t const minNumAllocators = MpiWrapper::min( numAllocators );

  if( numAllocators != minNumAllocators )
  {
    GEOS_WARNING( "Not all ranks have created the same number of umpire allocators, cannot compute high water marks." );
    return;
  }

  // Loop over the allocators.
  unsigned MAX_NAME_LENGTH = 100;

  TableData tableData;
  for( string const & allocatorName : allocatorNames )
  {
    // Skip umpire internal allocators.
    if( allocatorName.rfind( "__umpire_internal", 0 ) == 0 )
      continue;

    GEOS_ERROR_IF_GT( allocatorName.size(), MAX_NAME_LENGTH );
    string allocatorNameFixedSize = allocatorName;
    allocatorNameFixedSize.resize( MAX_NAME_LENGTH, '\0' );
    string allocatorNameMinChars = string( MAX_NAME_LENGTH, '\0' );

    // Make sure that each rank is looking at the same allocator.
    MpiWrapper::allReduce( allocatorNameFixedSize.c_str(), &allocatorNameMinChars.front(), MAX_NAME_LENGTH, MPI_MIN, MPI_COMM_GEOSX );
    GEOS_LOG_RANK_0( " allocatorNameBuffer " <<  allocatorNameFixedSize << "allocatorNameMinCharsBuffer" << allocatorNameMinChars << std::endl );
    if( allocatorNameFixedSize != allocatorNameMinChars )
    {
      GEOS_WARNING( "Not all ranks have an allocator named " << allocatorNameFixedSize << ", cannot compute high water mark." );
      continue;
    }

    // Get the total number of bytes allocated with this allocator across ranks.
    // This is a little redundant since
    std::size_t const mark = rm.getAllocator( allocatorName ).getHighWatermark();

    string const minMarkValue = GEOS_FMT( "{} ({}%)",
                                          LvArray::system::calculateSize( MpiWrapper::min( mark ) ) );
    string const maxMarkValue = GEOS_FMT( "{} ({}%)",
                                          LvArray::system::calculateSize( MpiWrapper::max( mark ) ));
    string const avgMarkValue = GEOS_FMT( "{} ({}%)",
                                          LvArray::system::calculateSize( MpiWrapper::max( mark ) ) );
    string const sumMarkValue = GEOS_FMT( "{} ({}%)",
                                          LvArray::system::calculateSize( MpiWrapper::sum( mark ) ) );

    tableData.addRow( allocatorName,
                      minMarkValue,
                      maxMarkValue,
                      avgMarkValue,
                      sumMarkValue );

    pushStatsIntoAdiak( allocatorName + " sum across ranks", mark );
    pushStatsIntoAdiak( allocatorName + " rank max", mark );
  }

  TableLayout const memoryStatLayout ( {"Umpire Pool", "Min (GB/%)\nover ranks", "Max (GB/%)\nover ranks", "Avg (GB/%)\nover ranks", "Sum (GB/%)\nover ranks" } );
  TableTextFormatter const memoryStatLog( memoryStatLayout );
  GEOS_LOG_RANK_0( memoryStatLog.toString( tableData ));

  // GEOS_LOG_RANK_0( allocatorLog.ToString( tableData )); // prettier-ignore

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupEnvironment( int argc, char * argv[] )
{
  setupMPI( argc, argv );
  setupLogger();
  setupLvArray();
  setupOpenMP();
  setupMKL();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cleanupEnvironment()
{
  LvArray::system::resetSignalHandling();
  finalizeLogger();
  addUmpireHighWaterMarks();
  finalizeCaliper();
  finalizeMPI();
}


} // namespace geos
