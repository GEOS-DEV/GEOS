/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#endif

#include "managers/ProblemManager.hpp"

using namespace geosx;
using namespace dataRepository;
namespace
{
int global_argc;
char** global_argv;
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>(global_argc)];
  for( int i=0 ; i<argc ; ++i )
  {
    global_argv[i] = argv[i];
    std::cout<<argv[i]<<std::endl;
  }

  return RUN_ALL_TESTS();
}

TEST(testXML,testXML)
{
  ProblemManager problemManager("ProblemManager",nullptr);

  problemManager.BuildDataStructure(nullptr);
  problemManager.SetDocumentationNodes( &problemManager );

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( global_argc, global_argv );
  // {
  //   dataRepository::ManagedGroup * commandLine = problemManager.GetGroup<ManagedGroup>(std::string("commandLine"));
  //   ViewWrapper<std::string>::rtype  inputFileName = commandLine->getData<std::string>(std::string("inputFileName"));
  //   inputFileName = "../../src/components/core/tests/xmlTests/basic_input.xml";
  // }
  problemManager.ParseInputFile();

}
