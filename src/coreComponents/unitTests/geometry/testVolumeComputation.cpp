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

#include "meshUtilities/ComputationalGeometry.hpp"

using namespace geosx;
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
    R1Tensor hex[8] = {R1Tensor(3), R1Tensor(3), R1Tensor(3), R1Tensor(3), 
                       R1Tensor(3), R1Tensor(3), R1Tensor(3), R1Tensor(3)};
    hex[0][0] = 0.; hex[0][1] = 0. ; hex[0][2]= 0;
    hex[1][0] = 1.; hex[1][1] = 0. ; hex[0][2]= 0;
    hex[2][0] = 0.; hex[3][1] = 1. ; hex[0][2]= 0;
    hex[3][0] = 1.; hex[2][1] = 1. ; hex[0][2]= 0;
    hex[4][0] = 0.; hex[4][1] = 0. ; hex[0][2]= 1;
    hex[5][0] = 1.; hex[5][1] = 0. ; hex[0][2]= 1;
    hex[6][0] = 0.; hex[7][1] = 1. ; hex[0][2]= 1;
    hex[7][0] = 1.; hex[6][1] = 1. ; hex[0][2]= 1;

    if( geosx::computationalGeometry::HexVolume( hex ) != 1.) {
        std::cout << geosx::computationalGeometry::HexVolume( hex ) << std::endl;
        GEOS_ERROR("Problem detected in the volume computation of hex");
    }


}
