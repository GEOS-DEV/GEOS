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

// forcefully enable asserts macros for this unit test
#define GEOS_ASSERT_ENABLED

#include "common/logger/ErrorHandling.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/DataContext.hpp"
#include "common/initializeEnvironment.hpp"

#include <gtest/gtest.h>
#include <filesystem>

using namespace geos;
using namespace dataRepository;

namespace fs = std::filesystem;

// declare a constant which value is the source file line (to predict the error file output).
#define GET_LINE( lineVar ) static size_t constexpr lineVar = __LINE__

// various dummy test values and contexts
double testMinPrecision = 1e-6;
double testMaxPrecision = 1e-3;
int testValue = 5;
DataFileContext const context = DataFileContext( "Base Test Class", "/path/to/file.xml", 23 );
DataFileContext const additionalContext = DataFileContext( "Additional Test Class", "/path/to/file.xml", 32 );
DataFileContext const importantAdditionalContext = DataFileContext( "Important Additional Test Class", "/path/to/file.xml", 64 );

/**
 * @brief begin a test with a local logger
 * @param errorLogger local error logger instance
 * @param filename output error filename
 */
void beginLocalLoggerTest( ErrorLogger & errorLogger, string_view filename )
{
  errorLogger.enableFileOutput( true );
  errorLogger.setOutputFilename( filename );
  errorLogger.createFile();
}

/**
 * @brief end the local logger test by reading the logger file output, comparing it to a reference, and removing it.
 * @param errorLogger local error logger instance
 * @param expectedFileBits reference file parts that must be in the logger file output
 */
void endLocalLoggerTest( ErrorLogger & errorLogger,
                         std::vector< string > expectedFileBits )
{
  auto const readFile = [] ( string_view filename ) {
    if( !fs::exists( filename ))
      throw std::runtime_error( "File not found: " + std::string( filename ) );

    std::ifstream file{ std::string( filename ) };
    if( !file.is_open())
      throw std::runtime_error( "Failed to open file: " + std::string( filename ) );

    std::stringstream buffer;
    std::string line;
    while( std::getline( file, line ))
      buffer << line << '\n';

    return buffer.str();
  };

  string_view filename = errorLogger.getOutputFilename();
  string fileContent = readFile( filename );
  bool testFailed = false;
  for( size_t i = 0; i < expectedFileBits.size(); ++i )
  {
    bool const foundFileBit = fileContent.find( expectedFileBits[i] ) != string::npos;
    EXPECT_TRUE( foundFileBit ) << "Expected bit not found (no." << i << "):\n"
                                << "-----------------------\n"
                                << expectedFileBits[i] << '\n'
                                << "-----------------------\n";
    testFailed |= !foundFileBit;
  }
  EXPECT_FALSE( testFailed ) << "Generated error file content:\n"
                             << "-----------------------\n"
                             << fileContent << '\n'
                             << "-----------------------\n";

  if( fs::exists( filename ) )
    fs::remove( filename );
}

TEST( ErrorHandling, testYamlFileWarningOutput )
{
  ErrorLogger g_errorLogger;  // Local overriding of global 'g_errorLogger' (to contain test macros effects to local scope)
  beginLocalLoggerTest( g_errorLogger, "warningTestOutput.yaml" );

  GET_LINE( line1 ); GEOS_WARNING( "Conflicting pressure boundary conditions" );

  GET_LINE( line2 ); GEOS_WARNING_IF_GT_MSG( testValue, testMaxPrecision, "Pressure value is too high." );

  GET_LINE( line3 ); GEOS_WARNING_CTX_IF( testValue == 5,
                                          GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.",
                                                    context.toString(), testMinPrecision, testMaxPrecision, testMinPrecision ),
                                          context, additionalContext );

  endLocalLoggerTest( g_errorLogger, {
    R"(errors:)",

    GEOS_FMT(
      R"(- type: Warning
    rank: 0
    message: >-
      Conflicting pressure boundary conditions
    sourceLocation:
      file: {}
      line: {})",
      __FILE__, line1 ),

    GEOS_FMT(
      R"(- type: Warning
    rank: 0
    message: >-
      Pressure value is too high.
    cause: >-
      Expected: testValue <= testMaxPrecision
      * testValue = 5
      * testMaxPrecision = 0.001
    sourceLocation:
      file: {}
      line: {})",
      __FILE__, line2 ),

    GEOS_FMT(
      R"(- type: Warning
    rank: 0
    message: >-
      Base Test Class (file.xml, l.23): option should be between 1e-06 and 0.001. A value of 1e-06 will be used.
    contexts:
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 32
    cause: >-
      Warning cause: testValue == 5
    sourceLocation:
      file: {}
      line: {})",
      __FILE__, line3 ),
  } );
}

TEST( ErrorHandling, testYamlFileExceptionOutput )
{
  ErrorLogger g_errorLogger;  // Local overriding of global 'g_errorLogger' (to contain test macros effects to local scope)
  beginLocalLoggerTest( g_errorLogger, "exceptionTestOutput.yaml" );
  size_t line1;

  // Stacked exception test (contexts must appear sorted by priority)
  try
  {
    line1 = __LINE__; GEOS_THROW_CTX_IF( testValue == 5,
                                         "Group " << context.toString() << " has no wrapper named" << std::endl,
                                         std::domain_error,
                                         context.getContextInfo().setPriority( 1 ) );
  }
  catch( std::domain_error const & ex )
  {
    string const errorMsg = "Table input error.\n";
    g_errorLogger.currentErrorMsg()
      .addToMsg( errorMsg )
      .addContextInfo( additionalContext.getContextInfo() )
      .addContextInfo( importantAdditionalContext.getContextInfo().setPriority( 2 ) );
  }
  g_errorLogger.flushErrorMsg( g_errorLogger.currentErrorMsg() );

  endLocalLoggerTest( g_errorLogger, {
    R"(errors:)",

    GEOS_FMT(
      R"(- type: Exception
    rank: 0
    message: >-
      Table input error.
      Group Base Test Class (file.xml, l.23) has no wrapper named
    contexts:
      - priority: 2
        inputFile: /path/to/file.xml
        inputLine: 64
      - priority: 1
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 32
    sourceLocation:
      file: {}
      line: {}
    sourceCallStack:)",
      __FILE__, line1 ),
    "- frame0: ",
    "- frame1: ",
    "- frame2: "
  } );
}

TEST( ErrorHandling, testYamlFileErrorOutput )
{
  ErrorLogger g_errorLogger; // Local overriding of global 'g_errorLogger' (to contain test macros effects to local scope)
  beginLocalLoggerTest( g_errorLogger, "errorTestOutput.yaml" );

  GET_LINE( line1 ); EXPECT_EXIT( GEOS_ERROR_IF_GT_MSG( testValue, testMaxPrecision,
                                                        GEOS_FMT( "{}: option should be lower than {}.",
                                                                  context.toString(), testMaxPrecision ),
                                                        context,
                                                        additionalContext,
                                                        importantAdditionalContext.getContextInfo().setPriority( 2 ) ),
                                  ::testing::ExitedWithCode( 1 ),
                                  ".*" );

  endLocalLoggerTest( g_errorLogger, {
    R"(errors:)",

    GEOS_FMT(
      R"(- type: Error
    rank: 0
    message: >-
      Base Test Class (file.xml, l.23): option should be lower than 0.001.
    contexts:
      - priority: 2
        inputFile: /path/to/file.xml
        inputLine: 64
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 32
    cause: >-
      Expected: testValue <= testMaxPrecision
      * testValue = 5
      * testMaxPrecision = 0.001
    sourceLocation:
      file: {}
      line: {}
    sourceCallStack:)",
      __FILE__, line1 ),
    "- frame0: ",
    "- frame1: ",
    "- frame2: "
  } );
}

#ifdef GEOS_ASSERT_ENABLED
TEST( ErrorHandling, testYamlFileAssertOutput )
{
  ErrorLogger g_errorLogger; // Local overriding of global 'g_errorLogger' (to contain test macros effects to local scope)
  beginLocalLoggerTest( g_errorLogger, "assertTestOutput.yaml" );

  GET_LINE( line1 ); EXPECT_EXIT( GEOS_ASSERT_MSG( testValue > testMinPrecision && testValue < testMaxPrecision,
                                                   GEOS_FMT( "{}: value should be between {} and {}, but is {}.",
                                                             context.toString(), testMinPrecision, testMaxPrecision, testValue ),
                                                   context,
                                                   additionalContext ),
                                  ::testing::ExitedWithCode( 1 ),
                                  ".*" );

  endLocalLoggerTest( g_errorLogger, {
    R"(errors:)",

    GEOS_FMT(
      R"(- type: Error
    rank: 0
    message: >-
      Base Test Class (file.xml, l.23): value should be between 1e-06 and 0.001, but is 5.
    contexts:
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 32
    cause: >-
      Expected: testValue > testMinPrecision && testValue < testMaxPrecision
    sourceLocation:
      file: {}
      line: {}
    sourceCallStack:)",
      __FILE__, line1 ),
    "- frame0: ",
    "- frame1: ",
    "- frame2: "
  } );
}
#endif

int main( int ac, char * av[] )
{
  ::testing::GTEST_FLAG( death_test_style ) = "threadsafe";
  ::testing::InitGoogleTest( &ac, av );
  geos::setupEnvironment( ac, av );
  int const result = RUN_ALL_TESTS();
  geos::cleanupEnvironment();
  return result;
}
