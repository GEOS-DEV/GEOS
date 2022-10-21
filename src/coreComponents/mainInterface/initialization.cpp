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

#include "initialization.hpp"

#include "common/DataTypes.hpp"
#include "common/Path.hpp"
#include "LvArray/src/system.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mainInterface/GeosxVersion.hpp"

// TPL includes
#include <optionparser.h>
#include <hdf5.h>
#include <vtkVersion.h>
#include <umpire/Umpire.hpp>
#include <conduit/conduit.hpp>
namespace geosx
{

/**
 * @class Arg a class inheriting from option::Arg that can parse a command line argument.
 */
struct Arg : public option::Arg
{
  /**
   * @brief Parse an unknown option. Unknown options aren't supported so this throws an error.
   * @param option the option to parse.
   * @return option::ARG_ILLEGAL.
   */
  static option::ArgStatus unknown( option::Option const & option, bool )
  {
    GEOSX_LOG_RANK( "Unknown option: " << option.name );
    return option::ARG_ILLEGAL;
  }

  /**
   * @brief Parse a non-empty string option.
   * @param option the option to parse.
   * @return option::ARK_OK if the parse was successful, option::ARG_ILLEGAL otherwise.
   */
  static option::ArgStatus nonEmpty( const option::Option & option, bool )
  {
    if((option.arg != nullptr) && (option.arg[0] != 0))
    {
      return option::ARG_OK;
    }

    GEOSX_LOG_RANK( "Error: " << option.name << " requires a non-empty argument!" );
    return option::ARG_ILLEGAL;
  }

  /**
   * @brief Parse a numeric string option.
   * @param option the option to parse.
   * @return option::ARK_OK if the parse was successful, option::ARG_ILLEGAL otherwise.
   */
  static option::ArgStatus numeric( const option::Option & option, bool )
  {
    char * endptr = nullptr;
    if((option.arg != nullptr) && strtol( option.arg, &endptr, 10 )) {}
    if((endptr != option.arg) && (*endptr == 0))
    {
      return option::ARG_OK;
    }

    GEOSX_LOG_RANK( "Error: " << option.name << " requires a long-int argument!" );
    return option::ARG_ILLEGAL;
  }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unique_ptr< CommandLineOptions > parseCommandLineOptions( int argc, char * * argv )
{
  std::unique_ptr< CommandLineOptions > commandLineOptions = std::make_unique< CommandLineOptions >();

  // Set the options structs and parse
  enum optionIndex
  {
    UNKNOWN,
    HELP,
    INPUT,
    RESTART,
    XPAR,
    YPAR,
    ZPAR,
    SCHEMA,
    NONBLOCKING_MPI,
    SUPPRESS_PINNED,
    PROBLEMNAME,
    OUTPUTDIR,
    TIMERS,
    TRACE_DATA_MIGRATION,
    PAUSE_FOR,
  };

  const option::Descriptor usage[] =
  {
    { UNKNOWN, 0, "", "", Arg::unknown, "USAGE: geosx -i input.xml [options]\n\nOptions:" },
    { HELP, 0, "?", "help", Arg::None, "\t-?, --help" },
    { INPUT, 0, "i", "input", Arg::nonEmpty, "\t-i, --input, \t Input xml filename (required)" },
    { RESTART, 0, "r", "restart", Arg::nonEmpty, "\t-r, --restart, \t Target restart filename" },
    { XPAR, 0, "x", "xpartitions", Arg::numeric, "\t-x, --x-partitions, \t Number of partitions in the x-direction" },
    { YPAR, 0, "y", "ypartitions", Arg::numeric, "\t-y, --y-partitions, \t Number of partitions in the y-direction" },
    { ZPAR, 0, "z", "zpartitions", Arg::numeric, "\t-z, --z-partitions, \t Number of partitions in the z-direction" },
    { SCHEMA, 0, "s", "schema", Arg::nonEmpty, "\t-s, --schema, \t Name of the output schema" },
    { NONBLOCKING_MPI, 0, "b", "use-nonblocking", Arg::None, "\t-b, --use-nonblocking, \t Use non-blocking MPI communication" },
    { PROBLEMNAME, 0, "n", "name", Arg::nonEmpty, "\t-n, --name, \t Name of the problem, used for output" },
    { SUPPRESS_PINNED, 0, "s", "suppress-pinned", Arg::None, "\t-s, --suppress-pinned, \t Suppress usage of pinned memory for MPI communication buffers" },
    { OUTPUTDIR, 0, "o", "output", Arg::nonEmpty, "\t-o, --output, \t Directory to put the output files" },
    { TIMERS, 0, "t", "timers", Arg::nonEmpty, "\t-t, --timers, \t String specifying the type of timer output" },
    { TRACE_DATA_MIGRATION, 0, "", "trace-data-migration", Arg::None, "\t--trace-data-migration, \t Trace host-device data migration" },
    { PAUSE_FOR, 0, "", "pause-for", Arg::numeric, "\t--pause-for, \t Pause geosx for a given number of seconds before starting execution" },
    { 0, 0, nullptr, nullptr, nullptr, nullptr }
  };

  argc -= ( argc > 0 );
  argv += ( argc > 0 );
  option::Stats stats( usage, argc, argv );
  option::Option options[ 100 ];//stats.options_max];
  option::Option buffer[ 100 ];//stats.buffer_max];
  option::Parser parse( usage, argc, argv, options, buffer );

  // Handle special cases
  bool const noXML = options[INPUT].count() == 0 && options[SCHEMA].count() == 0;
  if( parse.error() || options[HELP] || (argc == 0) || noXML )
  {
    int columns = getenv( "COLUMNS" ) ? atoi( getenv( "COLUMNS" )) : 120;
    option::printUsage( fwrite, stdout, usage, columns );

    if( options[HELP] )
    {
      throw NotAnError();
    }

    GEOSX_THROW( "Bad command line arguments.", InputError );
  }

  // Iterate over the remaining inputs
  for( int ii=0; ii<parse.optionsCount(); ++ii )
  {
    option::Option & opt = buffer[ii];
    switch( opt.index() )
    {
      case UNKNOWN:
      {}
      break;
      case HELP:
      {}
      break;
      case INPUT:
      {
        commandLineOptions->inputFileNames.emplace_back( opt.arg );
      }
      break;
      case RESTART:
      {
        commandLineOptions->restartFileName = opt.arg;
        commandLineOptions->beginFromRestart = true;
      }
      break;
      case XPAR:
      {
        commandLineOptions->xPartitionsOverride = std::stoi( opt.arg );
        commandLineOptions->overridePartitionNumbers = 1;
      }
      break;
      case YPAR:
      {
        commandLineOptions->yPartitionsOverride = std::stoi( opt.arg );
        commandLineOptions->overridePartitionNumbers = 1;
      }
      break;
      case ZPAR:
      {
        commandLineOptions->zPartitionsOverride = std::stoi( opt.arg );
        commandLineOptions->overridePartitionNumbers = 1;
      }
      break;
      case NONBLOCKING_MPI:
      {
        commandLineOptions->useNonblockingMPI = true;
      }
      break;
      case SUPPRESS_PINNED:
      {
        commandLineOptions->suppressPinned = true;
      }
      break;
      case SCHEMA:
      {
        commandLineOptions->schemaName = opt.arg;
      }
      break;
      case PROBLEMNAME:
      {
        commandLineOptions->problemName = opt.arg;
      }
      break;
      case OUTPUTDIR:
      {
        commandLineOptions->outputDirectory = opt.arg;
      }
      break;
      case TIMERS:
      {
        commandLineOptions->timerOutput = opt.arg;
      }
      break;
      case TRACE_DATA_MIGRATION:
      {
        commandLineOptions->traceDataMigration = true;
      }
      break;
      case PAUSE_FOR:
      {
        // we should store this in commandLineOptions and sleep in main
        integer const duration = std::stoi( opt.arg );
        GEOSX_LOG_RANK_0( "Paused for " << duration << " s" );
        std::this_thread::sleep_for( std::chrono::seconds( duration ) );
      }
      break;
    }
  }

  if( commandLineOptions->problemName.empty() && options[INPUT].count() > 0 )
  {
    string & inputFileName = commandLineOptions->inputFileNames[0];
    if( inputFileName.length() > 4 && inputFileName.substr( inputFileName.length() - 4, 4 ) == ".xml" )
    {
      string::size_type start = inputFileName.find_last_of( '/' ) + 1;
      if( start >= inputFileName.length())
      {
        start = 0;
      }
      commandLineOptions->problemName.assign( inputFileName, start, inputFileName.length() - 4 - start );
    }
    else
    {
      commandLineOptions->problemName.assign( inputFileName );
    }
  }

  return commandLineOptions;
}


void outputVersionInfo()
{
#if defined(__clang__)
  GEOSX_LOG_RANK_0( "  - clang version: " << 
                    __clang_major__<<"."<<
                    __clang_minor__<<"."<<
                    __clang_patchlevel__ );
#if defined(__apple_build_version__)
  GEOSX_LOG_RANK_0( "    apple clang version: " << 
                    __apple_build_version__);
#endif
#if defined(__ibmxl_vrm__)
  GEOSX_LOG_RANK_0( "  - IBM compiler version: " <<
                    __ibmxl_version__<<"."<<
                    __ibmxl_release__<<"."<<
                    __ibmxl_modification__<<"."<<
                    __ibmxl_ptf_fix_level__ );
#endif
#elif defined(__GNUC__)
  GEOSX_LOG_RANK_0( "  - gcc version: " << 
                    __GNUC__<<"."<<
                    __GNUC_MINOR__<<"."<<
                    __GNUC_PATCHLEVEL__ );
#endif

#if defined( GEOSX_USE_CUDA )
  GEOSX_LOG_RANK_0( "  - cuda version: " <<
                    CUDA_VERSION/1000<<"."<<
                    CUDA_VERSION/10%100 );
#endif

  {
#if defined(_OPENMP)
    GEOSX_LOG_RANK_0( "  - openmp version: "<<_OPENMP )
#endif
  }


  {
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int len;

    MPI_Get_library_version( version, &len );
    GEOSX_LOG_RANK_0( "  - MPI version: " << version );

  }

  {
    unsigned majnum, minnum, relnum;
    H5get_libversion( &majnum, &minnum, &relnum );
    GEOSX_LOG_RANK_0( "  - HDF5 version: " << majnum<<"."<<minnum<<"."<<relnum );
  }

  {
    static const char * vtkVersionString = vtkVersion::GetVTKVersion();
    GEOSX_LOG_RANK_0( "  - VTK version: " << vtkVersionString );

  }

  {
    GEOSX_LOG_RANK_0( "  - RAJA version: " << 
                      RAJA_VERSION_MAJOR<<"."<<
                      RAJA_VERSION_MINOR<<"."<<
                      RAJA_VERSION_PATCHLEVEL );

    GEOSX_LOG_RANK_0( "  - UMPIRE version: " << 
                      umpire::get_major_version()<<"."<<
                      umpire::get_minor_version()<<"."<<
                      umpire::get_patch_version() );
  }

  {
    GEOSX_LOG_RANK_0( "  - Conduit Version: "<<CONDUIT_VERSION );
  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unique_ptr< CommandLineOptions > basicSetup( int argc, char * argv[], bool const parseCommandLine )
{
  setupEnvironment( argc, argv );

  GEOSX_LOG_RANK_0( "GEOSX version: " << getVersion() );
  outputVersionInfo();

  setupLAI();

  if( parseCommandLine )
  {
    return parseCommandLineOptions( argc, argv );
  }
  else
  {
    return std::make_unique< CommandLineOptions >();
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void basicCleanup()
{
  finalizeLAI();
  cleanupEnvironment();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string getVersion()
{
#if defined(GEOSX_GIT_BRANCH) && defined(GEOSX_GIT_HASH)
  return GEOSX_VERSION_FULL " (" GEOSX_GIT_BRANCH ", sha1: " GEOSX_GIT_HASH ")";
#else
  return GEOSX_VERSION_FULL;
#endif
}


} // namespace geosx
