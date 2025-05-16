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
 * @file PartitionerBase.cpp
 */


#include "PartitionerBase.hpp"

#include "common/MpiWrapper.hpp"


namespace geos
{


void PartitionerBase::buildNeighbors()
{
  m_neighbors.clear();

  for( int const rank : m_neighborsRank )
  {
    m_neighbors.push_back( NeighborCommunicator( rank ));
  }

  m_numFirstOrderNeighbors = getNumNeighbors();
}


void PartitionerBase::appendNeighborsOfNeighbors()
{
  for( int const neighborRank : buildNeighborsOfNeighborsRank() )
  {
    m_neighborsRank.emplace_back( neighborRank );
    m_neighbors.emplace_back( neighborRank );
  }
}

std::set< int > PartitionerBase::buildNeighborsOfNeighborsRank()
{

  // Send this list of neighbors to all neighbors.
  int neighborsTag = 54;
  stdVector< MPI_Request > requests( getNumNeighbors(), MPI_REQUEST_NULL );

  for( int i = 0; i < getNumNeighbors(); ++i )
  {
    MpiWrapper::iSend( getNeighborsRank().data(), getNumNeighbors(), getNeighborsRank()[i],
                       neighborsTag, MPI_COMM_GEOS, &requests[ i ] );
  }

  // This set will contain the second (neighbor of) neighbors ranks.
  std::set< int > neighborsOfNeighborsList;
  for( int i = 0; i < getNumNeighbors(); ++i )
  {
    array1d< int > neighborsOfOneNeighbor;
    MpiWrapper::recv( neighborsOfOneNeighbor, getNeighborsRank()[ i ], neighborsTag, MPI_COMM_GEOS, MPI_STATUS_IGNORE );

    // Insert the neighbors of the current neighbor into the set of second neighbors.
    neighborsOfNeighborsList.insert( neighborsOfOneNeighbor.begin(), neighborsOfOneNeighbor.end() );
  }

  // Ensure all send operations are complete
  MpiWrapper::waitAll( requests.size(), requests.data(), MPI_STATUSES_IGNORE );


  // Remove yourself and all the first neighbors from the second neighbors.
  neighborsOfNeighborsList.erase( MpiWrapper::commRank() );
  for( int const & neighborRank : getNeighborsRank() )
  {
    neighborsOfNeighborsList.erase( neighborRank );
  }

  return neighborsOfNeighborsList;
}


} // namespace geos
