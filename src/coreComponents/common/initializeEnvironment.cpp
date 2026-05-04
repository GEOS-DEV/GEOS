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

#include "initializeEnvironment.hpp"

#include "TimingMacros.hpp"
#include "Path.hpp"
#include "LvArray/src/system.hpp"
#include "common/LifoStorageCommon.hpp"
#include "logger/ErrorHandling.hpp"
#include "logger/ExternalErrorHandler.hpp"
#include <umpire/TypedAllocator.hpp>
// TPL includes
#include <umpire/ResourceManager.hpp>
#include <umpire/Allocator.hpp>
#include <umpire/strategy/AllocationStrategy.hpp>
#include "umpire/util/MemoryResourceTraits.hpp"
#include "umpire/util/Platform.hpp"

#if defined( GEOS_USE_CALIPER )
#include <caliper/cali-manager.h>
#if defined( GEOS_USE_ADIAK )
#include <adiak.hpp>
#endif
#endif

// System includes
#include <iomanip>

#if defined( GEOS_USE_MKL )
#include <mkl.h>
#endif

#if defined( GEOS_USE_OPENMP )
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
#ifdef GEOS_USE_MPI
  logger::InitializeLogger( MPI_COMM_GEOS );
#else
  logger::InitializeLogger();
#endif

  { // setup error handling (using LvArray helper system functions)

    ExternalErrorHandler::instance().enableStderrPipeDeviation( true );

    ///// set external error handling behaviour /////
    ExternalErrorHandler::instance().setErrorHandling( []( string_view errorMsg,
                                                           string_view detectionLocation )
    {
      // Filter out INFO level messages from external libraries (e.g., VTK)
      // ( error / signal lambda would calls either an error function or an info function, depending on a filtering function )
      if( ExternalErrorHandler::isNotAnErrorMsg( errorMsg ) )
      {
        // Just print the message without error formatting
        GEOS_LOG( errorMsg );
        return;
      }
      else
      {
        std::string const stackHistory = LvArray::system::stackTrace( true );
        DiagnosticMsg diagnosticMsg;
        ErrorLogger::global().flushErrorMsg( DiagnosticMsgBuilder::init( diagnosticMsg,
                                                                         MsgType::Error, errorMsg,
                                                                         ::geos::logger::internal::g_rank )
                                               .addCallStackInfo( stackHistory )
                                               .addDetectionLocation( detectionLocation )
                                               .getDiagnosticMsg() );

        // we do not terminate the program as 1. the error could be non-fatal, 2. there may be more messages to output.
      }
    } );

    ///// set signal handling behaviour /////
    LvArray::system::setSignalHandling( []( int const signal )
    {
      // Disable signal handling to prevent catching exit signal (infinite loop)
      LvArray::system::setSignalHandling( nullptr );

      // first of all, external error can await to be output, we must output them
      ExternalErrorHandler::instance().flush( "before signal error output" );

      // error message output
      std::string const stackHistory = LvArray::system::stackTrace( true );
      DiagnosticMsg diagnosticMsg;
      ErrorLogger::global().flushErrorMsg( DiagnosticMsgBuilder::init( diagnosticMsg,
                                                                       MsgType::ExternalError, "",
                                                                       ::geos::logger::internal::g_rank )
                                             .addSignal( signal )
                                             .addCallStackInfo( stackHistory )
                                             .getDiagnosticMsg() );

      // call program termination
      LvArray::system::callErrorHandler();
    } );

    ///// set Post-Handled Error behaviour /////
    LvArray::system::setErrorHandler( []()
    {
  #if defined( GEOS_USE_MPI )
      int mpi = 0;
      MPI_Initialized( &mpi );
      if( mpi )
      {
        MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
      }
  #endif
      std::abort();
    } );
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeLogger()
{
  logger::FinalizeLogger();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupLvArray()
{
#if defined(GEOS_USE_FPE)
  LvArray::system::setFPE();
#else
  LvArray::system::disableFloatingPointExceptions( FE_ALL_EXCEPT );
#endif

  /* Disable chai callbacks by default */
  chai::ArrayManager::getInstance()->disableCallbacks();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupMKL()
{
#ifdef GEOS_USE_MKL
  GEOS_LOG_RANK_0( "MKL max threads: " << mkl_get_max_threads() );
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupOpenMP()
{
#ifdef GEOS_USE_OPENMP
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

  MPI_COMM_GEOS = MpiWrapper::commDup( MPI_COMM_WORLD );

  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
  {
    // Can't use logging macros prior to logger init
    std::cout << "Num ranks: " << MpiWrapper::commSize( MPI_COMM_GEOS ) << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeMPI( bool inError )
{
  if( !inError )
  {
    MpiWrapper::finalize();
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupCUDA()
{
#if defined( GEOS_USE_CUDA ) && defined( GEOS_CUDA_STACK_SIZE )
  size_t const stackSize = (GEOS_CUDA_STACK_SIZE) * 1024;
  if( 0 < stackSize )
  {
    cudaError_t status = cudaDeviceSetLimit( cudaLimitStackSize, stackSize );
    GEOS_ERROR_IF( status != cudaSuccess,
                   GEOS_FMT( "Failed to set CUDA stack size. Error {}: {}", status, cudaGetErrorString( status ) ) );
  }
#endif
}

#if defined( GEOS_USE_CALIPER )

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupCaliper( cali::ConfigManager & caliperManager,
                   CommandLineOptions const & commandLineOptions )
{
  caliperManager.add( commandLineOptions.timerOutput.c_str() );
  GEOS_ERROR_IF( caliperManager.error(),
                 GEOS_FMT( "Caliper config error: {}", caliperManager.error_msg() ) );
  caliperManager.start();

#if defined( GEOS_USE_ADIAK )
#if defined( GEOS_USE_MPI )
  adiak::init( &MPI_COMM_GEOS );
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
#if defined( GEOS_USE_MPI )
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

  adiak::value( "build type", GEOS_CMAKE_BUILD_TYPE );
  adiak::value( "compilation date", __DATE__ );

  // OpenMP info
#if defined( GEOS_USE_OPENMP )
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
  GEOS_ERROR_IF_NE( hipSuccess, hipRuntimeGetVersion( &hipRuntimeVersion ) );
  GEOS_ERROR_IF_NE( hipSuccess, hipDriverGetVersion( &hipDriverVersion ) );
#else
  adiak::value( "HIP", "Off" );
#endif
  adiak::value( "HIP runtime version", hipRuntimeVersion );
  adiak::value( "HIP driver version", hipDriverVersion );
#endif // defined( GEOS_USE ADIAK )
}
#endif // defined( GEOS_USE_CALIPER )

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeCaliper()
{
#if defined( GEOS_USE_CALIPER )and defined( GEOS_USE_ADIAK )
  adiak::fini();
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupEnvironment( int argc, char * argv[] )
{
  setupCUDA();
  setupMPI( argc, argv );
  setupLogger();
  setupLvArray();
  setupOpenMP();
  setupMKL();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cleanupEnvironment( bool inError )
{
  LvArray::system::resetSignalHandling();
  finalizeLogger();
  finalizeCaliper();
  finalizeMPI( inError );
}

} // namespace geos
