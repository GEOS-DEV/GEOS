/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Logger.cpp
 */

// Source includes
#include "Logger.hpp"
#include "common/MpiWrapper.hpp"
#include "common/Path.hpp"

namespace geos
{

namespace logger
{

namespace internal
{

int g_rank = 0;
int g_n_ranks = 1;
std::string g_rankString = "?";
std::ostream * g_rankStream = nullptr;

} // namespace internal

#ifdef GEOS_USE_MPI

void InitializeLogger( MPI_Comm mpi_comm, const std::string & rankOutputDir )
{
  MPI_Comm_rank( mpi_comm, &internal::g_rank );
  MPI_Comm_size( mpi_comm, &internal::g_n_ranks );

  internal::g_rankString = std::to_string( internal::g_rank );

  if( rankOutputDir != "" )
  {
    if( internal::g_rank == 0 )
    {
      makeDirsForPath( rankOutputDir );
    }

    MpiWrapper::barrier( mpi_comm );
    std::string outputFilePath = rankOutputDir + "/rank_" + internal::g_rankString + ".out";
    internal::g_rankStream = new std::ofstream( outputFilePath );
  }
  else
  {
    internal::g_rankStream = &std::cout;
  }
}

#endif

void InitializeLogger( const std::string & rankOutputDir )
{
  if( rankOutputDir != "" )
  {
    makeDirsForPath( rankOutputDir );

    std::string outputFilePath = rankOutputDir + "/rank_" + internal::g_rankString + ".out";
    internal::g_rankStream = new std::ofstream( outputFilePath );
  }
  else
  {
    internal::g_rankStream = &std::cout;
  }
}

void FinalizeLogger()
{
  if( internal::g_rankStream != &std::cout )
  {
    delete internal::g_rankStream;
  }

  internal::g_rankStream = nullptr;
}

} // namespace logger

} // namespace geos
