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
#include "LvArray/src/system.hpp"
#include "common/logger/ExternalErrorHandler.hpp"
#include "gtest/gtest.h"
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

// redeging logger instance to test macros with a local instance (to prevent any side effect)
#undef GEOS_GLOBAL_LOGGER
#define GEOS_GLOBAL_LOGGER testErrorLogger

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
                         stdVector< string > expectedFileBits )
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
  ErrorLogger testErrorLogger;

  beginLocalLoggerTest( testErrorLogger, "warningTestOutput.yaml" );

  GET_LINE( line1 ); GEOS_WARNING( "Conflicting pressure boundary conditions" );
  GET_LINE( line2 ); GEOS_WARNING_IF_GT_MSG( testValue, testMaxPrecision, "Pressure value is too high." );

  string const warningMsg = GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.",
                                      context.toString(), testMinPrecision, testMaxPrecision, testMinPrecision );
  GET_LINE( line3 ); GEOS_WARNING_IF( testValue == 5, warningMsg, context, additionalContext );
  endLocalLoggerTest( testErrorLogger, {
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
        description: Base Test Class (file.xml, l.23)
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        description: Additional Test Class (file.xml, l.32)
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
  ErrorLogger testErrorLogger;
  string const file = "exceptionTestOutput.yaml";
  beginLocalLoggerTest( testErrorLogger, file );
  size_t line1;

  // Stacked exception test (contexts must appear sorted by priority)
  try
  {
    line1 = __LINE__; GEOS_THROW_IF( testValue == 5, "Empty Group", geos::DomainError, context.getContextInfo().setPriority( 3 ) );
  }
  catch( geos::DomainError const & ex )
  {
    string const errorMsg = "Table input error.\n";
    testErrorLogger.modifyCurrentExceptionMessage()
      .addToMsg( errorMsg )
      .addContextInfo( additionalContext.getContextInfo() )
      .addContextInfo( importantAdditionalContext.getContextInfo().setPriority( 2 ) );
    string const whatExpected = GEOS_FMT( "***** Exception\n"
                                          "***** LOCATION: {}:{}\n"
                                          "***** Error cause: testValue == 5\n"
                                          "***** Rank 0\n"
                                          "***** Message from {}:\n"
                                          "Empty Group",
                                          testErrorLogger.getCurrentExceptionMsg().m_file, line1,
                                          testErrorLogger.getCurrentExceptionMsg().m_contextsInfo.front().m_formattedContext );
    GEOS_ERROR_IF_EQ_MSG( string( ex.what()).find( whatExpected ), string::npos,
                          "The error message was not containing the expected sequence.\n" <<
                          "  Error message :\n" << ex.what() <<
                          "  expected sequence :\n" << whatExpected );
  }

  testErrorLogger.flushCurrentExceptionMessage();
  endLocalLoggerTest( testErrorLogger, {
    R"(errors:)",
    GEOS_FMT(
      R"(- type: Exception
    rank: 0
    message: >-
      Table input error.
      Empty Group
    contexts:
      - priority: 3
        description: Base Test Class (file.xml, l.23)
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 2
        description: Important Additional Test Class (file.xml, l.64)
        inputFile: /path/to/file.xml
        inputLine: 64
      - priority: 0
        description: Additional Test Class (file.xml, l.32)
        inputFile: /path/to/file.xml
        inputLine: 32
    cause: >-
      Error cause: testValue == 5
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
  ErrorLogger testErrorLogger;

  beginLocalLoggerTest( testErrorLogger, "errorTestOutput.yaml" );

  EXPECT_EXIT( GEOS_ERROR_IF_GT_MSG( testValue, testMaxPrecision,
                                     GEOS_FMT( "{}: option should be lower than {}.",
                                               context.toString(), testMaxPrecision ),
                                     context,
                                     additionalContext,
                                     importantAdditionalContext.getContextInfo().setPriority( 2 ) ),
               ::testing::ExitedWithCode( 1 ),
               ".*" );

  endLocalLoggerTest( testErrorLogger, {
    R"(errors:)",

    // we won't test the line index for this test as it cannot be a one-liner.
    R"(- type: Error
    rank: 0
    message: >-
      Base Test Class (file.xml, l.23): option should be lower than 0.001.
    contexts:
      - priority: 2
        description: Important Additional Test Class (file.xml, l.64)
        inputFile: /path/to/file.xml
        inputLine: 64
      - priority: 0
        description: Base Test Class (file.xml, l.23)
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        description: Additional Test Class (file.xml, l.32)
        inputFile: /path/to/file.xml
        inputLine: 32
    cause: >-
      Expected: testValue <= testMaxPrecision
      * testValue = 5
      * testMaxPrecision = 0.001
    sourceLocation:)",
    "  file: ",
    "  line: ",
    "sourceCallStack:",
    "- frame0: ",
    "- frame1: ",
    "- frame2: "
  } );
}


TEST( ErrorHandling, testLogFileExceptionOutput )
{
  ErrorLogger testErrorLogger;

  size_t line1;
  // test that the first context with default priority comes after the most important one, but before
  // a second context with default priotity
  try
  {
    line1 = __LINE__; GEOS_THROW_IF( testValue == 5, "Empty Group", geos::DomainError, context );
  }
  catch( geos::DomainError const & ex )
  {
    string const errorMsg = "Table input error.\n";
    testErrorLogger.modifyCurrentExceptionMessage()
      .addToMsg( errorMsg )
      .addContextInfo( additionalContext.getContextInfo() )
      .addContextInfo( importantAdditionalContext.getContextInfo().setPriority( 2 ));
    string const streamExpected = GEOS_FMT(
      "***** Exception\n"
      "***** LOCATION: {}:{}\n"
      "***** {}\n"
      "***** Rank 0\n"
      "***** Message from {}:\n"
      "{}\n"
      "***** Additional contexts:\n"
      "***** - {}\n"
      "***** - {}\n",
      testErrorLogger.getCurrentExceptionMsg().m_file, line1,
      testErrorLogger.getCurrentExceptionMsg().m_cause,
      importantAdditionalContext.toString(),
      testErrorLogger.getCurrentExceptionMsg().m_msg,
      context.toString(),
      additionalContext.toString(),
      testErrorLogger.getCurrentExceptionMsg().m_sourceCallStack );
    std::ostringstream oss;
    ErrorLogger::formatMsgForLog( testErrorLogger.getCurrentExceptionMsg(), oss );
    GEOS_ERROR_IF_EQ_MSG( oss.str().find( streamExpected ), string::npos,
                          "The error message was not containing the expected sequence.\n" <<
                          "The error message was not containing the expected sequence.\n" <<
                          "  Error message :\n" <<oss.str() <<
                          "  expected sequence :\n" << streamExpected );
  }
}

TEST( ErrorHandling, testStdException )
{
  ErrorLogger testErrorLogger;

  try
  {

    throw std::invalid_argument( "received negative value" );
  }
  catch( std::exception & e )
  {

    testErrorLogger.initCurrentExceptionMessage( MsgType::Exception, e.what(),
                                                 ::geos::logger::internal::g_rank )
      .addCallStackInfo( LvArray::system::stackTrace( true ) );

    std::ostringstream oss;
    ErrorLogger::formatMsgForLog( testErrorLogger.getCurrentExceptionMsg(), oss );
    string const streamExpected = GEOS_FMT(
      "***** Exception\n"
      "***** Rank 0\n"
      "***** Message :\n"
      "{}\n",
      testErrorLogger.getCurrentExceptionMsg().m_msg );
    GEOS_ERROR_IF_EQ_MSG( oss.str().find( streamExpected ), string::npos,
                          "The error message was not containing the expected sequence.\n" <<
                          "The error message was not containing the expected sequence.\n" <<
                          "  Error message :\n" << oss.str() <<
                          "  expected sequence :\n" << streamExpected );
  }
}

TEST( ErrorHandling, testYamlFileAssertOutput )
{
  ErrorLogger testErrorLogger;

  beginLocalLoggerTest( testErrorLogger, "assertTestOutput.yaml" );

  EXPECT_EXIT( GEOS_ASSERT_MSG( testValue > testMinPrecision && testValue < testMaxPrecision,
                                GEOS_FMT( "{}: value should be between {} and {}, but is {}.",
                                          context.toString(), testMinPrecision, testMaxPrecision, testValue ),
                                context,
                                additionalContext ),
               ::testing::ExitedWithCode( 1 ),
               ".*" );

  endLocalLoggerTest( testErrorLogger, {
    R"(errors:)",

    // we won't test the line index for this test as it cannot be a one-liner.
    R"(- type: Error
    rank: 0
    message: >-
      Base Test Class (file.xml, l.23): value should be between 1e-06 and 0.001, but is 5.
    contexts:
      - priority: 0
        description: Base Test Class (file.xml, l.23)
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        description: Additional Test Class (file.xml, l.32)
        inputFile: /path/to/file.xml
        inputLine: 32
    cause: >-
      Expected: testValue > testMinPrecision && testValue < testMaxPrecision
    sourceLocation:)",
    "  file: ",
    "  line: ",
    "sourceCallStack:",
    "- frame0: ",
    "- frame1: ",
    "- frame2: "
  } );
}


TEST( ErrorHandling, VerifySignalHandlerLogs )
{
  ErrorLogger testErrorLogger;

  beginLocalLoggerTest( testErrorLogger, "errors.yaml" );

  setupLogger();

  bool signalHappened = false;
  LvArray::system::setErrorHandler( [&]()
  {
    endLocalLoggerTest( testErrorLogger, {
      R"(- type: ExternalError
    rank: 0
    message: >-
      Signal encountered (no. 2): Interrupt
    contexts:
      - priority: 0
        description: Signal (detected from Signal Handler)
        detectionLocation: Signal handler
        signal: 2
    sourceCallStack:)",
      "- frame0: ",
      "- frame1: ",
      "- frame2: "
    } );

    signalHappened = true;
  } );

  ErrorLogger::global().enableFileOutput( true );

  raise( SIGINT );

  EXPECT_TRUE( signalHappened );

  setupLogger();
}

int main( int ac, char * av[] )
{
  ::testing::GTEST_FLAG( death_test_style ) = "threadsafe";
  ::testing::InitGoogleTest( &ac, av );
  geos::setupEnvironment( ac, av );
  int const result = RUN_ALL_TESTS();
  geos::cleanupEnvironment( );
  return result;
}
