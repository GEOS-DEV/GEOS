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

static constexpr std::array< std::string_view, 5 > expectedFileBits = {
  R"(errors: 

  - type: Warning
    rank: 0
    message: >-
      Conflicting pressure boundary conditions
    sourceLocation:
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 194
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start 

  - type: Warning
    rank: 0
    message: >-
      Pressure value is too small.
    sourceLocation:
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 196
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start 

  - type: Warning
    rank: 0
    message: >-
      Base Test Class (file.xml, l.23): option should be between 1e-06 and 0.001. A value of 1e-06 will be used.
    contexts:
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 12
    sourceLocation:
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 197
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start 

  - type: Exception
    rank: 0
    message: >-
      Table input error.
      Group Base Test Class (file.xml, l.23) has no wrapper named
    contexts:
      - priority: 2
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 23
      - priority: 0
        inputFile: /path/to/file.xml
        inputLine: 12
    sourceLocation:
      file: )",
  R"(src/coreComponents/dataRepository/unitTests/testErrorHandling.cpp
      line: 203
    sourceCallStack:
      - frame0:  void testing::internal::HandleExceptionsInMethodIfSupported<testing::Test, void>(testing::Test*, void (testing::Test::*)(), char const*) 
      - frame1:  testing::Test::Run() 
      - frame2:  testing::TestInfo::Run() 
      - frame3:  testing::TestSuite::Run() 
      - frame4:  testing::internal::UnitTestImpl::RunAllTests() 
      - frame5:  testing::UnitTest::Run() 
      - frame6:  main 
      - frame7:  __libc_start_main 
      - frame8:  _start 
)" };

static constexpr std::string_view exceptionFormat =
  R"(    message: >-
      Table input error.
      Group Base Test Class (file.xml, l.23) has no wrapper named
    contexts:
      - priority: 2
        inputFile: /path/to/file.xml
        inputLine: 23
)";

TEST( ErrorHandling, testYaml )
{
  g_errorLogger.setOutputFilename( filename );
  g_errorLogger.enableFileOutput( true );
  double minPrecision = 1e-6;
  double maxPrecision = 1e-3;
  int x = 5;

  DataFileContext const context = DataFileContext( "Base Test Class", "/path/to/file.xml", 23 );
  DataFileContext const additionalContext = DataFileContext( "Additional Test Class", "/path/to/file.xml", 12 );

  if( g_errorLogger.isOutputFileEnabled() )
  {
    g_errorLogger.createFile();
  }

  GEOS_WARNING( "Conflicting pressure boundary conditions" );

  GEOS_WARNING_IF( x == 5, "Pressure value is too small." );
  GEOS_WARNING_CTX_IF( x == 5,
                       GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.",
                                 context.toString(), minPrecision, maxPrecision, minPrecision ),
                       context, additionalContext );
  try
  {
    GEOS_THROW_CTX_IF( x == 5,
                       "Group " << context.toString() << " has no wrapper named" << std::endl,
                       std::domain_error,
                       context, additionalContext );
  }
  catch( std::domain_error const & ex )
  {
    string const errorMsg = "Table input error.\n";
    g_errorLogger.currentErrorMsg()
      .addToMsg( errorMsg )
      .addContextInfo( context.getContextInfo().setPriority( 2 ) );
  }

  if( g_errorLogger.isOutputFileEnabled() )
  {
    g_errorLogger.flushErrorMsg( g_errorLogger.currentErrorMsg() );
  }

  std::string fileContent = readFile();

  for( size_t i = 0; i < expectedFileBits.size(); ++i )
  {
    auto it = fileContent.find( expectedFileBits[i] );
    EXPECT_NE( it, std::string::npos ) << "Expected bit not found: " << expectedFileBits[i];
  }

  std::string additionalExceptionInformation = readFile( 65, 72 );
  EXPECT_EQ( additionalExceptionInformation, exceptionFormat );

  EXPECT_EXIT( GEOS_ERROR_CTX_IF( x == 5,
                                  GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.PID  {}",
                                            context.toString(), minPrecision, maxPrecision, minPrecision, getpid() ),
                                  context, additionalContext ),
               ::testing::ExitedWithCode( 1 ),
               ".*" );

  g_errorLogger = ErrorLogger{};
  removeFile();
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
