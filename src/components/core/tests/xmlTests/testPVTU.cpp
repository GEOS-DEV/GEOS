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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include "gtest/gtest.h"
#if USE_MPI
#include <mpi.h>
#endif

#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "fileIO/pvtu/PvtuFile.hpp"
using namespace geosx;
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


TEST(testPVTU,testPVTU)
{
#if USE_MPI
    MPI_Init(0, nullptr);
#endif
    PvtuFile vtupFile;
    vtupFile.Load("4layers.pvtu");
#if USE_MPI
    MPI_Finalize();
#endif
}
