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
 * @file DistributedMatrixStatistics.cpp
 */

#include "DistributedMatrixStatistics.hpp"

#include "common/MpiWrapper.hpp"
#include "codingUtilities/Utilities.hpp"

#include <algorithm>

namespace geos
{

namespace
{

enum Metric : localIndex
{
  firstRow,
  numRows,
  numNonzeros,
  numOffRankNonzeros,
  numHaloReceiveDofs,
  numHaloSendDofs,
  numReceiveNeighbors,
  numSendNeighbors,
  maxRowNonzeros,
  numEmptyRows,
  numMetrics
};

int findOwner( globalIndex const column,
               arrayView1d< globalIndex const > const & ownership,
               int const numRanks )
{
  int lowerRank = 0;
  int upperRank = numRanks;
  while( lowerRank < upperRank )
  {
    int const middleRank = lowerRank + ( upperRank - lowerRank ) / 2;
    globalIndex const rowEnd = ownership[2 * middleRank + 1];
    if( column < rowEnd )
    {
      upperRank = middleRank;
    }
    else
    {
      lowerRank = middleRank + 1;
    }
  }

  if( lowerRank == numRanks || column < ownership[2 * lowerRank] )
  {
    return -1;
  }
  return lowerRank;
}

} // namespace

DistributedMatrixStatistics
computeDistributedMatrixStatistics( CRSMatrixView< real64 const, globalIndex const > const & localMatrix,
                                    globalIndex const firstLocalRow,
                                    MPI_Comm const comm )
{
  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );
  globalIndex const localRows = LvArray::integerConversion< globalIndex >( localMatrix.numRows() );
  globalIndex const localRowEnd = firstLocalRow + localRows;

  array1d< globalIndex > localOwnership( 2 );
  localOwnership[0] = firstLocalRow;
  localOwnership[1] = localRowEnd;
  array1d< globalIndex > ownership;
  MpiWrapper::allGather( localOwnership.toViewConst(), ownership, comm );

  for( int owner = 1; owner < numRanks; ++owner )
  {
    GEOS_ERROR_IF_NE_MSG( ownership[2 * owner], ownership[2 * ( owner - 1 ) + 1],
                          "Distributed matrix row ownership must be contiguous" );
  }

  // Matrix assembly may leave the CRS buffers in the backend's memory space.
  // Diagnostics are explicitly opt-in, so migrate the pattern to the host for inspection.
  localMatrix.move( hostMemorySpace, false );

  globalIndex localNonzeros = 0;
  globalIndex localOffRankNonzeros = 0;
  globalIndex localMaxRowNonzeros = 0;
  globalIndex localEmptyRows = 0;
  std::vector< globalIndex > remoteColumns;

  for( localIndex localRow = 0; localRow < localMatrix.numRows(); ++localRow )
  {
    localIndex const rowNonzeros = localMatrix.numNonZeros( localRow );
    localNonzeros += LvArray::integerConversion< globalIndex >( rowNonzeros );
    localMaxRowNonzeros = std::max( localMaxRowNonzeros,
                                    LvArray::integerConversion< globalIndex >( rowNonzeros ) );
    localEmptyRows += rowNonzeros == 0;

    arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( localRow );
    for( localIndex entry = 0; entry < rowNonzeros; ++entry )
    {
      globalIndex const column = columns[entry];
      if( column < firstLocalRow || column >= localRowEnd )
      {
        ++localOffRankNonzeros;
        remoteColumns.push_back( column );
      }
    }
  }

  std::sort( remoteColumns.begin(), remoteColumns.end() );
  remoteColumns.erase( std::unique( remoteColumns.begin(), remoteColumns.end() ), remoteColumns.end() );

  array1d< globalIndex > receiveCountsByOwner;
  receiveCountsByOwner.resizeDefault( numRanks, 0 );
  for( globalIndex const column : remoteColumns )
  {
    int const owner = findOwner( column, ownership.toViewConst(), numRanks );
    GEOS_ERROR_IF_LT_MSG( owner, 0, GEOS_FMT( "Matrix column {} has no owning rank", column ) );
    GEOS_ERROR_IF_EQ_MSG( owner, rank,
                          GEOS_FMT( "Column {} was classified as off-rank but is locally owned", column ) );
    ++receiveCountsByOwner[owner];
  }

  array1d< globalIndex > sendCountsByRequester( numRanks );
#ifdef GEOS_USE_MPI
  MPI_CHECK_ERROR( MPI_Alltoall( receiveCountsByOwner.data(),
                                 1,
                                 internal::getMpiType< globalIndex >(),
                                 sendCountsByRequester.data(),
                                 1,
                                 internal::getMpiType< globalIndex >(),
                                 comm ) );
#else
  for( int peer = 0; peer < numRanks; ++peer )
  {
    sendCountsByRequester[peer] = receiveCountsByOwner[peer];
  }
#endif

  globalIndex localHaloSendDofs = 0;
  globalIndex localReceiveNeighbors = 0;
  globalIndex localSendNeighbors = 0;
  for( int peer = 0; peer < numRanks; ++peer )
  {
    localReceiveNeighbors += receiveCountsByOwner[peer] > 0;
    globalIndex const sendCount = sendCountsByRequester[peer];
    localHaloSendDofs += sendCount;
    localSendNeighbors += sendCount > 0;
  }

  array1d< globalIndex > localMetrics( static_cast< localIndex >( numMetrics ) );
  localMetrics[firstRow] = firstLocalRow;
  localMetrics[numRows] = localRows;
  localMetrics[numNonzeros] = localNonzeros;
  localMetrics[numOffRankNonzeros] = localOffRankNonzeros;
  localMetrics[numHaloReceiveDofs] = LvArray::integerConversion< globalIndex >( remoteColumns.size() );
  localMetrics[numHaloSendDofs] = localHaloSendDofs;
  localMetrics[numReceiveNeighbors] = localReceiveNeighbors;
  localMetrics[numSendNeighbors] = localSendNeighbors;
  localMetrics[maxRowNonzeros] = localMaxRowNonzeros;
  localMetrics[numEmptyRows] = localEmptyRows;

  array1d< globalIndex > allMetrics;
  MpiWrapper::allGather( localMetrics.toViewConst(), allMetrics, comm );

  DistributedMatrixStatistics result;
  result.ranks.resize( numRanks );
  for( int metricRank = 0; metricRank < numRanks; ++metricRank )
  {
    globalIndex const * const metrics = allMetrics.data() + metricRank * numMetrics;
    MatrixRankStatistics & rankStatistics = result.ranks[metricRank];
    rankStatistics.firstRow = metrics[firstRow];
    rankStatistics.numRows = metrics[numRows];
    rankStatistics.numNonzeros = metrics[numNonzeros];
    rankStatistics.numOffRankNonzeros = metrics[numOffRankNonzeros];
    rankStatistics.numHaloReceiveDofs = metrics[numHaloReceiveDofs];
    rankStatistics.numHaloSendDofs = metrics[numHaloSendDofs];
    rankStatistics.numReceiveNeighbors = metrics[numReceiveNeighbors];
    rankStatistics.numSendNeighbors = metrics[numSendNeighbors];
    rankStatistics.maxRowNonzeros = metrics[maxRowNonzeros];
    rankStatistics.numEmptyRows = metrics[numEmptyRows];
  }

  return result;
}

} // namespace geos
