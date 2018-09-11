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
 * Logger.cpp
 *
 *  Created on: Aug 31, 2017
 *      Author: settgast
 */

#include "Logger.hpp"
#include "DataTypes.hpp"
#include <mpi.h>
#include <stdlib.h>
#include "stackTrace.hpp"
#include<iostream>

namespace geosx
{

namespace logger
{

int rank = 0;

void geos_abort()
{
  cxx_utilities::handler1(EXIT_FAILURE);
#ifdef GEOSX_USE_MPI
  int mpi = 0;
  MPI_Initialized( &mpi );
  if ( mpi )
  {
    MPI_Abort( MPI_COMM_GEOSX, EXIT_FAILURE );
  }
  else
#endif
  {
    exit( EXIT_FAILURE );
  }
}

} /* namespace logger */

} /* namespace geosx */
