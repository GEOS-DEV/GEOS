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
std::string rankString = "0";

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

  internal::rankString = std::to_string( internal::rank );

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

    std::string outputFilePath = rankOutputDir + "/rank_" + internal::rankString + ".out";
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

    std::string outputFilePath = rankOutputDir + "/rank_" + internal::rankString + ".out";
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
