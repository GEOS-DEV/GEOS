/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "gtest/gtest.h"

#include "managers/ProblemManager.hpp"

using namespace geosx;
using namespace dataRepository;
namespace
{
int global_argc;
char * * global_argv;
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

#ifdef GEOSX_USE_MPI
  int rank;
  MPI_Init( &argc, &argv );
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
  logger::InitializeLogger( MPI_COMM_GEOSX );
#else
  logger::InitializeLogger();
#endif

  global_argc = argc;
  global_argv = new char *[global_argc];
  for( int i = 0; i < argc; ++i )
  {
    global_argv[i] = argv[i];
    std::cout<<argv[i]<<std::endl;
  }

  int const result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

#ifdef GEOSX_USE_MPI
  MPI_Finalize();
#endif

  return result;
}


TEST( testXML, testXML )
{
  ProblemManager problemManager( "Problem", nullptr );

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( global_argc, global_argv );
  // {
  //   dataRepository::Group * commandLine = problemManager.GetGroup<Group>(std::string("commandLine"));
  //   Wrapper<std::string>::rtype  inputFileName = commandLine->getData<std::string>(std::string("inputFileName"));
  //   inputFileName = "../../src/components/core/tests/xmlTests/basic_input.xml";
  // }
  problemManager.ParseInputFile();
}
