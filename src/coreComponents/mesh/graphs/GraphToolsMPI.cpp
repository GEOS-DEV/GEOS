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
 * @file GraphToolsMPI.cpp
 *
 */

#include "GraphToolsMPI.hpp"
#include <numeric>


namespace geos
{
namespace graph
{

std::pair< std::vector< camp::idx_t >, std::vector< camp::idx_t > >
scatterGraphData( const std::vector< camp::idx_t > & xadj,
                  const std::vector< camp::idx_t > & adjncy,
                  MPI_Comm comm )
{
  //int const rank = MpiWrapper::commRank( comm );
  //int const size = MpiWrapper::commSize( comm );
  int rank, size;
  MPI_Comm_rank( comm, &rank );
  MPI_Comm_size( comm, &size );

  // Check that the number of ranks is the same as the size of xadj-1
  if( rank == 0 )
  {
    //GEOS_ASSERT_EQ_MSG( size, static_cast< int >(xadj.size()) - 1, "Number of ranks does not match the size of xadj-1" );
  }

  std::vector< int > sendCounts;
  std::vector< int > displacements;
  std::vector< camp::idx_t > xadjToScatter;

  if( rank == 0 )
  {
    sendCounts.resize( size );
    displacements.resize( size );

    // Calculate send counts and displacements for MPI_Scatterv
    for( int i = 0; i < size; ++i )
    {
      sendCounts[i] = static_cast< int >(xadj[i + 1] - xadj[i]);
      displacements[i] = static_cast< int >(xadj[i]);
    }

    xadjToScatter.resize( 2*size );
    for( int i = 0; i < size; ++i )
    {
      xadjToScatter[2 * i] = xadj[i];
      xadjToScatter[2 * i + 1] = xadj[i + 1];
    }
  }

  std::vector< camp::idx_t > localXadj( 2 ); // Each rank will have two elements in localXadj: xadj[i] and xadj[i+1]
  MPI_Scatter( xadjToScatter.data(), 2, MPI_LONG, localXadj.data(), 2, MPI_LONG, 0, comm );

  int localSize;
  MPI_Scatter( sendCounts.data(), 1, MPI_INT, &localSize, 1, MPI_INT, 0, comm );

  std::vector< camp::idx_t > localAdjncy( localSize );
  MPI_Scatterv( adjncy.data(), sendCounts.data(), displacements.data(), MPI_LONG,
                localAdjncy.data(), localSize, MPI_LONG, 0, comm );

  return {localXadj, localAdjncy};
}



std::pair< std::vector< camp::idx_t >, std::vector< camp::idx_t > >
gatherGraphData( const std::vector< camp::idx_t > & localXadj,
                 const std::vector< camp::idx_t > & localAdjncy,
                 MPI_Comm comm )
{
  //int const rank = MpiWrapper::commRank( comm );
  //int const size = MpiWrapper::commSize( comm );
  int rank, size;
  MPI_Comm_rank( comm, &rank );
  MPI_Comm_size( comm, &size );

  // Determine the number of elements to send
  int nCounts = (rank < size - 1) ? static_cast< int >(localXadj.size() - 1) : static_cast< int >(localXadj.size());

  // Gather the counts of elements from each rank
  std::vector< int > recvCounts( size );
  MPI_Gather( &nCounts, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, comm );

  // Calculate displacements for the gathered data
  std::vector< int > displacements( size );
  if( rank == 0 )
  {
    displacements[0] = 0;
    for( int i = 1; i < size; ++i )
    {
      displacements[i] = displacements[i - 1] + recvCounts[i - 1];
    }
  }

  // Resize the xadj vector on rank 0 to hold the gathered data
  std::vector< camp::idx_t > xadj;
  if( rank == 0 )
  {
    int totalSize = std::accumulate( recvCounts.begin(), recvCounts.end(), 0 );
    xadj.resize( totalSize );
  }

  // Gather the xadj data from all ranks
  MPI_Gatherv( localXadj.data(), nCounts, MPI_LONG,
               xadj.data(), recvCounts.data(), displacements.data(), MPI_LONG, 0, comm );


  // Gather the counts of elements from each rank for localAdjncy
  int localAdjncySize = static_cast< int >(localAdjncy.size());
  MPI_Gather( &localAdjncySize, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, comm );

  // Calculate displacements for the gathered data
  if( rank == 0 )
  {
    displacements[0] = 0;
    for( int i = 1; i < size; ++i )
    {
      displacements[i] = displacements[i - 1] + recvCounts[i - 1];
    }
  }

  // Resize the adjncy vector on rank 0 to hold the gathered data
  std::vector< camp::idx_t > adjncy;
  if( rank == 0 )
  {
    int totalSize = std::accumulate( recvCounts.begin(), recvCounts.end(), 0 );
    adjncy.resize( totalSize );
  }

  // Gather the adjncy data from all ranks
  MPI_Gatherv( localAdjncy.data(), localAdjncySize, MPI_LONG,
               adjncy.data(), recvCounts.data(), displacements.data(), MPI_LONG, 0, comm );

  return {xadj, adjncy};
}



std::vector< camp::idx_t > createXadjFromAdjncy( const std::vector< camp::idx_t > & localAdjncy, MPI_Comm comm )
{
  int rank, size;
  MPI_Comm_rank( comm, &rank );
  MPI_Comm_size( comm, &size );

  // Gather the counts of elements from each rank for localAdjncy
  int localAdjncySize = static_cast< int >(localAdjncy.size());
  std::vector< int > adjncyCounts( size );
  MPI_Allgather( &localAdjncySize, 1, MPI_INT, adjncyCounts.data(), 1, MPI_INT, comm );

  // Calculate xadj
  std::vector< camp::idx_t > xadj( size + 1, 0 );
  std::partial_sum( adjncyCounts.begin(), adjncyCounts.end(), xadj.begin() + 1 );

  // Prepare data to scatter
  std::vector< camp::idx_t > xadjToScatter( 2 * size );
  for( int i = 0; i < size; ++i )
  {
    xadjToScatter[2 * i] = xadj[i];
    xadjToScatter[2 * i + 1] = xadj[i + 1];
  }

  // Scatter the xadj data
  std::vector< camp::idx_t > localXadj( 2 );
  MPI_Scatter( xadjToScatter.data(), 2, MPI_LONG, localXadj.data(), 2, MPI_LONG, 0, comm );

  return localXadj;
}



std::vector< int > createVertexGlobalID( const std::vector< camp::idx_t > & localXadj, MPI_Comm comm )
{
  int rank, size;
  MPI_Comm_rank( comm, &rank );
  MPI_Comm_size( comm, &size );

  int localNumVertices = static_cast< int >(localXadj.size()) - 1;

  std::vector< int > allNumVertices( size );
  MPI_Allgather( &localNumVertices, 1, MPI_INT, allNumVertices.data(), 1, MPI_INT, comm );

  // Compute displacement for the current rank
  std::vector< int > displacements( size, 0 );
  std::partial_sum( allNumVertices.begin(), allNumVertices.end() - 1, displacements.begin() + 1 );
  int displacement = displacements[rank];

  std::vector< int > globalID( localNumVertices );
  std::iota( globalID.begin(), globalID.end(), displacement );

  return globalID;
}



} // namespace graph
} // namespace geos
