/*
 * Logger.cpp
 *
 *  Created on: Aug 31, 2017
 *      Author: settgast
 */

#include "Logger.hpp"
#include <mpi.h>
#include <stdlib.h>
#include "stackTrace.hpp"
#include<iostream>

namespace geosx
{

void geos_abort( std::string message )
{
  std::cerr<<message<<std::endl;
  cxx_utilities::handler1(EXIT_FAILURE);
#if USE_MPI == 1
  int mpi = 0;
  MPI_Initialized( &mpi );
  if ( mpi )
  {
    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
  }
  else
#endif
  {
    exit( EXIT_FAILURE );
  }
}

}
