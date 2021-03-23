/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_MANAGERS_INITIALIZATION_HPP_
#define GEOSX_MANAGERS_INITIALIZATION_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

// TPL includes
#ifdef GEOSX_USE_CALIPER
  #include <adiak.hpp>

//Forward declaration of cali::ConfigManager.
namespace cali
{
class ConfigManager;
}
#endif

namespace geosx
{

/**
 * CommandLineOptions class containing the parsed command line options.
 */
struct CommandLineOptions
{
  /// The path to the input xml.
  string inputFileName;

  /// True iff restarting from the middle of an existing run.
  bool beginFromRestart = false;

  /// The path to the restart file, if specified.
  string restartFileName;

  /// The number of partitions in the x direction.
  integer xPartitionsOverride;

  /// The number of partitions in the y direction.
  integer yPartitionsOverride;

  /// The number of partitions in the z direction.
  integer zPartitionsOverride;

  /// True if using the partition override.
  integer overridePartitionNumbers = false;

  /// True if processing mpi communications in any order.
  /// But leads to non-reproducible results.
  integer useNonblockingMPI = false;

  /// True iff supress the use of pinned memory buffers
  /// ( if available ) for MPI communication.
  /// Generally only used by the integration tests.
  integer suppressPinned = false;

  /// The name of the schema.
  string schemaName;

  /// The name of the problem being run.
  string problemName;

  /// The directory to put all output.
  string outputDirectory = ".";

  /// The string used to initialize caliper.
  string timerOutput = "";

  /// Suppress logging of host-device data migration.
  integer suppressMoveLogging = false;
};

/**
 * @brief Parse the command line options and populate @p commandLineOptions with the results.
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @return The command line options.
 */
std::unique_ptr< CommandLineOptions > parseCommandLineOptions( int argc, char * * argv );

/**
 * @brief Perform the basic GEOSX initialization and optionally parse the command line input.
 * @param [in] argc The number of command line arguments.
 * @param [in,out] argv The command line arguments.
 * @param [in] parseCommandLine True iff the command line options should be parsed.
 * @return The command line options, if @c parseCommandLine is @c false then the returned value
 *   is default constructed (empty).
 */
std::unique_ptr< CommandLineOptions > basicSetup( int argc, char * argv[], bool const parseCommandLine=false );

/**
 * @brief Perform the basic GEOSX cleanup.
 */
void basicCleanup();

/**
 * @brief Initialize the logger.
 */
void setupLogger();

/**
 * @brief Finalize the logger.
 */
void finalizeLogger();

/**
 * @brief Finalize Caliper and Adiak.
 */
void finalizeCaliper();

/**
 * @brief Setup the LvArray library. This initializes signal handling
 *        and the floating point environment.
 */
void setupLvArray();

/**
 * @brief Setup MKL if in use.
 */
void setupMKL();

/**
 * @brief Setup OpenMP.
 */
void setupOpenMP();

/**
 * @brief Setup MPI.
 * @param [in] argc the number of command line arguments.
 * @param [in,out] argv the command line arguments.
 */
void setupMPI( int argc, char * argv[] );

/**
 * @brief Finalize MPI.
 */
void finalizeMPI();

#if defined( GEOSX_USE_CALIPER )

/**
 * @brief Setup Caliper and Adiak.
 * @param caliperManager The Caliper ConfigManager to initialize.
 * @param commandLineOptions The command line options.
 */
void setupCaliper( cali::ConfigManager & caliperManager,
                   CommandLineOptions const & commandLineOptions );


#endif

/**
 * @brief Compute the sum, mean, min, and max of @p value across ranks and push
 *        them into Adiak using @p name.
 * @tparam T The type of @p value.
 * @param name The name to use when adding the stats to Adiak.
 * @param value The value to compute the statistics of.
 */
template< typename T >
void pushStatsIntoAdiak( string const & name, T const value )
{
#if defined( GEOSX_USE_CALIPER ) && !defined(__APPLE__)
  // Apple clang doesn't like adiak.
  T const total = MpiWrapper::sum( value );
  adiak::value( name + " sum", total );
  adiak::value( name + " mean", double( total ) / MpiWrapper::commSize() );
  adiak::value( name + " min", MpiWrapper::min( value ) );
  adiak::value( name + " max", MpiWrapper::max( value ) );
#else
  GEOSX_UNUSED_VAR( name );
  GEOSX_UNUSED_VAR( value );
#endif
}

} // namespace geosx

#endif // GEOSX_MANAGERS_INITIALIZATION_HPP_
