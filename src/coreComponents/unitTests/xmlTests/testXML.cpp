/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
#define __null nullptr
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

#ifdef GEOSX_USE_MPI
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
#endif

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>(global_argc)];
  for( int i=0 ; i<argc ; ++i )
  {
    global_argv[i] = argv[i];
    std::cout<<argv[i]<<std::endl;
  }

  int const result = RUN_ALL_TESTS();

#ifdef GEOSX_USE_MPI
  MPI_Finalize();
#endif

  return result;
}


TEST(testXML,testXML)
{
  ProblemManager problemManager("ProblemManager",nullptr);

  problemManager.BuildDataStructure(nullptr);
  problemManager.SetDocumentationNodes();

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( global_argc, global_argv );
  // {
  //   dataRepository::ManagedGroup * commandLine = problemManager.GetGroup<ManagedGroup>(std::string("commandLine"));
  //   ViewWrapper<std::string>::rtype  inputFileName = commandLine->getData<std::string>(std::string("inputFileName"));
  //   inputFileName = "../../src/components/core/tests/xmlTests/basic_input.xml";
  // }
  problemManager.ParseInputFile();

}
