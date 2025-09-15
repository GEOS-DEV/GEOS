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
#include "common/logger/ErrorHandling.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/DataContext.hpp"
#include "common/initializeEnvironment.hpp"

#include <gtest/gtest.h>
#include <filesystem>

using namespace geos;
using namespace dataRepository;

static constexpr std::string_view filename = "errorsOutput.yaml";

namespace fs = std::filesystem;

void removeFile()
{
  if( fs::exists( filename ) )
  {
    fs::remove( filename );
  }
}


std::string readFile( std::optional< size_t > startLine = std::nullopt,
                      std::optional< size_t > endLine = std::nullopt )
{
  if( !fs::exists( filename ))
  {
    throw std::runtime_error( "File not found: " + std::string( filename ) );
  }

  std::ifstream file{ std::string( filename ) };
  if( !file.is_open())
  {
    throw std::runtime_error( "Failed to open file: " + std::string( filename ) );
  }

  std::stringstream buffer;
  std::string line;
  size_t currentLine = 0;

  while( std::getline( file, line ))
  {
    if((!startLine || currentLine >= *startLine) &&
       (!endLine || currentLine < *endLine))
    {
      buffer << line << '\n';
    }
    currentLine++;
  }

  return buffer.str();
}

// separated file bits, which allow us to ignore the absolute path of the workspace ("file:" attribute).
// note: "line:" attribute of "sourceLocation:" need to be manually updated if test code changes (to
//       verify they are correctly reported)
static constexpr std::array< std::string_view, 12 > expectedFileBits = {
  R"(errors:)",

  R"(- type: Warning
    rank: 0
    message: >-
      Conflicting pressure boundary conditions
    sourceLocation:
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 233
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start)",

  R"(- type: Warning
    rank: 0
    message: >-
      Pressure value is too small.
    sourceLocation:
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 235
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start)",

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
    sourceLocation:
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 236
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start)",

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
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 244
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start)",

  R"(- type: Error
    rank: 0
    message: >-
      Base Test Class (file.xml, l.23): option should be between 1e-06 and 0.001. A value of 1e-06 will be used.
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
    sourceLocation:
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 259
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start)"
};

static constexpr std::string_view exceptionFormat =
  R"(    message: >-
      Table input error.
      Group Base Test Class (file.xml, l.23) has no wrapper named
    contexts:
      - priority: 2
        inputFile: /path/to/file.xml
        inputLine: 23
)";

double testMinPrecision = 1e-6;
double testMaxPrecision = 1e-3;
int testValue = 5;

TEST( ErrorHandling, testYamlFileOutputFormat )
{
  { // building the error yaml test case file
    // Local overriding of global 'g_errorLogger' (to contain test macros effects to local scope)
    ErrorLogger g_errorLogger;

    g_errorLogger.enableFileOutput( true );
    g_errorLogger.setOutputFilename( filename );

    DataFileContext const context = DataFileContext( "Base Test Class", "/path/to/file.xml", 23 );
    DataFileContext const additionalContext = DataFileContext( "Additional Test Class", "/path/to/file.xml", 32 );
    DataFileContext const importantAdditionalContext = DataFileContext( "Additional Test Class", "/path/to/file.xml", 64 );

    if( g_errorLogger.isOutputFileEnabled() )
    {
      g_errorLogger.createFile();
    }

    // Warning tests
    GEOS_WARNING( "Conflicting pressure boundary conditions" );

    GEOS_WARNING_IF( testValue == 5, "Pressure value is too small." );
    GEOS_WARNING_CTX_IF( testValue == 5,
                         GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.",
                                   context.toString(), testMinPrecision, testMaxPrecision, testMinPrecision ),
                         context, additionalContext );

    // Stacked exception test (contexts must appear sorted by priority)
    try
    {
      GEOS_THROW_CTX_IF( testValue == 5,
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

    EXPECT_EXIT( GEOS_ERROR_IF( testValue == 5,
                                GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.",
                                          context.toString(), testMinPrecision, testMaxPrecision, testMinPrecision ),
                                context,
                                additionalContext,
                                importantAdditionalContext.getContextInfo().setPriority( 2 ) ),
                 ::testing::ExitedWithCode( 1 ),
                 ".*" );
  }


  { // read back yaml file and check its content
    std::string fileContent = readFile();
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

    removeFile();
  }
}

int main( int ac, char * av[] )
{
  ::testing::GTEST_FLAG( death_test_style ) = "threadsafe";
  ::testing::InitGoogleTest( &ac, av );
  geos::setupEnvironment( ac, av );
  int const result = RUN_ALL_TESTS();
  geos::cleanupEnvironment();
  return result;
}
