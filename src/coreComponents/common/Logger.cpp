/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
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
#include "Path.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geos
{

/**
 * @brief Insert an exception message in another one.
 * @param originalMsg original exception message (i.e. thrown from LVARRAY_THROW or GEOSX_THROW)
 * @param msgToInsert message to insert at the top of the originalMsg
 */
std::string InsertExMsg( std::string const & originalMsg, std::string const & msgToInsert )
{
  std::string newMsg( originalMsg );

  size_t insertPos = 0;
  // for readability purposes, we try to insert the message after the "***** Rank N: " or after "***** " instead of at the top.
  static auto constexpr rankLogStart =  "***** Rank ";
  static auto constexpr rankLogEnd =  ": ";
  static auto constexpr simpleLogStart =  "***** ";
  if( ( insertPos = newMsg.find( rankLogStart ) ) != std::string::npos )
  {
    insertPos = newMsg.find( rankLogEnd, insertPos + stringutilities::cstrlen( rankLogStart ) )
                + stringutilities::cstrlen( rankLogEnd );
  }
  else if( ( insertPos = newMsg.find_last_of( simpleLogStart ) ) != std::string::npos )
  {
    insertPos += stringutilities::cstrlen( simpleLogStart );
  }
  newMsg.insert( insertPos, msgToInsert );
  return newMsg;
}

InputError::InputError( std::exception const & subException, std::string const & msgToInsert ):
  std::runtime_error( InsertExMsg( subException.what(), msgToInsert ) )
{}

namespace logger
{

namespace internal
{

int rank = 0;
std::string rankString = "0";

int n_ranks = 1;

std::ostream * rankStream = nullptr;

#ifdef GEOSX_USE_MPI
MPI_Comm comm;
#endif

} // namespace internal

#ifdef GEOSX_USE_MPI

void InitializeLogger( MPI_Comm mpi_comm, const std::string & rankOutputDir )
{
  internal::comm = mpi_comm;
  MPI_Comm_rank( mpi_comm, &internal::rank );
  MPI_Comm_size( mpi_comm, &internal::n_ranks );

  internal::rankString = std::to_string( internal::rank );

  if( rankOutputDir != "" )
  {
    if( internal::rank == 0 )
    {
      makeDirsForPath( rankOutputDir );
    }

    MPI_Barrier( mpi_comm );
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
    makeDirsForPath( rankOutputDir );

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

} // namespace geos
