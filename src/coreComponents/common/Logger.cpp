/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/**
 * @file Logger.cpp
 */

// Source includes
#include "common/Logger.hpp"

namespace geosx
{

namespace logger
{

namespace internal
{

int rank = 0;

int n_ranks = 1;

std::ostream * rankStream = nullptr;

#ifdef USE_MPI
MPI_Comm comm;
#endif

} // namespace internal

#ifdef USE_MPI

void InitializeLogger( MPI_Comm mpi_comm, const std::string & rankOutputDir )
{
  internal::comm = mpi_comm;
  MPI_Comm_rank( mpi_comm, &internal::rank );
  MPI_Comm_size( mpi_comm, &internal::n_ranks );

  if( rankOutputDir != "" )
  {
    if( internal::rank != 0 )
    {
      MPI_Barrier( mpi_comm );
    }
    else
    {
      std::string cmd = "mkdir -p " + rankOutputDir;
      int ret = std::system( cmd.c_str());
      if( ret != 0 )
      {
        GEOSX_ERROR( "Failed to initialize Logger: command '" << cmd << "' exited with code " << std::to_string( ret ));
      }
      MPI_Barrier( mpi_comm );
    }

    std::string outputFilePath = rankOutputDir + "/rank_" + std::to_string( internal::rank ) + ".out";
    internal::rankStream = new std::ofstream( outputFilePath );
  }
  else
  {
    internal::rankStream = &std::cout;
  }
}

#endif

void InitializeLogger( const std::string & rankOutputDir )
{
  if( rankOutputDir != "" )
  {
    std::string cmd = "mkdir -p " + rankOutputDir;
    int ret = std::system( cmd.c_str());
    if( ret != 0 )
    {
      GEOSX_LOG( "Failed to initialize Logger: command '" << cmd << "' exited with code " << std::to_string( ret ));
      abort();
    }

    std::string outputFilePath = rankOutputDir + "/rank_" + std::to_string( internal::rank ) + ".out";
    internal::rankStream = new std::ofstream( outputFilePath );
  }
  else
  {
    internal::rankStream = &std::cout;
  }
}

void FinalizeLogger()
{
  if( internal::rankStream != &std::cout )
  {
    delete internal::rankStream;
  }

  internal::rankStream = nullptr;
}

} // namespace logger

} // namespace geosx
