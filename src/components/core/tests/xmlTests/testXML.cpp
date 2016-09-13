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

int global_argc;
char** global_argv;

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  global_argc = argc;
  global_argv = new char*[argc];
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

  problemManager.Registration(nullptr);

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( global_argc, global_argv );
  problemManager.ParseInputFile();

}
