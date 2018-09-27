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

#ifdef GEOSX_USE_ATK

#include "slic/GenericOutputStream.hpp"

#ifdef GEOSX_USE_MPI
#include "slic/LumberjackStream.hpp"
#endif

#endif

namespace geosx
{

namespace logger
{

int rank = 0;


void InitializeLogger()
{
#ifdef GEOSX_USE_ATK
  axom::slic::initialize();
  axom::slic::setLoggingMsgLevel( axom::slic::message::Debug );

#ifdef GEOSX_USE_MPI
  std::string format =  std::string( "***********************************\n" )+
                        std::string( "MESSAGE=<MESSAGE>\n" ) +
                        std::string( "\t<TIMESTAMP>\n\n" ) +
                        std::string( "\tLEVEL=<LEVEL>\n" ) +
                        std::string( "\tRANKS=<RANK>\n") +
                        std::string( "\tFILE=<FILE>\n" ) +
                        std::string( "\tLINE=<LINE>\n" ) +
                        std::string( "***********************************\n" );

  const int ranks_limit = 5;
  axom::slic::LumberjackStream * const stream = new axom::slic::LumberjackStream(&std::cout, MPI_COMM_GEOSX, ranks_limit, format);
#else /* #ifdef GEOSX_USE_MPI */
  std::string format =  std::string( "***********************************\n" )+
                        std::string( "MESSAGE=<MESSAGE>\n" ) +
                        std::string( "\t<TIMESTAMP>\n\n" ) +
                        std::string( "\tLEVEL=<LEVEL>\n" ) +
                        std::string( "\tFILE=<FILE>\n" ) +
                        std::string( "\tLINE=<LINE>\n" ) +
                        std::string( "***********************************\n" );

  axom::slic::GenericOutputStream * const stream = new axom::slic::GenericOutputStream(&std::cout, format );
#endif /* #ifdef GEOSX_USE_MPI */
  axom::slic::addStreamToAllMsgLevels( stream );
#endif /* #ifdef GEOSX_USE_ATK */
}

void FinalizeLogger()
{
#ifdef GEOSX_USE_ATK
  axom::slic::flushStreams();
  axom::slic::finalize();
#endif

}

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
