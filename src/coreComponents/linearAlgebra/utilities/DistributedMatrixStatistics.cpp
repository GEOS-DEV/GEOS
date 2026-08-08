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
#include <cmath>
#include <fstream>
#include <vector>

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

void writeDistributedMatrixMarket( CRSMatrixView< real64 const, globalIndex const > const & localMatrix,
                                   globalIndex const firstLocalRow,
                                   globalIndex const numGlobalRows,
                                   globalIndex const numGlobalColumns,
                                   string const & filenamePrefix,
                                   MPI_Comm const comm )
{
  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );
  globalIndex const numLocalRows = LvArray::integerConversion< globalIndex >( localMatrix.numRows() );

  localMatrix.move( hostMemorySpace, false );
  globalIndex localNonzeros = 0;
  for( localIndex row = 0; row < localMatrix.numRows(); ++row )
  {
    localNonzeros += LvArray::integerConversion< globalIndex >( localMatrix.numNonZeros( row ) );
  }

  string const pieceFilename = GEOS_FMT( "{}_rank_{:06}.mtx", filenamePrefix, rank );
  std::vector< char > outputBuffer( 1024 * 1024 );
  std::ofstream output;
  output.rdbuf()->pubsetbuf( outputBuffer.data(), outputBuffer.size() );
  output.open( pieceFilename );
  GEOS_ERROR_IF( !output, GEOS_FMT( "Unable to open distributed matrix piece {}", pieceFilename ) );
  output << "%%MatrixMarket matrix coordinate real general\n";
  output << GEOS_FMT( "% distributed rank {} of {}; zero-based owned rows [{}, {})\n",
                      rank, numRanks, firstLocalRow, firstLocalRow + numLocalRows );
  output << GEOS_FMT( "{} {} {}\n", numGlobalRows, numGlobalColumns, localNonzeros );

  char line[96];
  globalIndex const maximumDimension = std::max( { numGlobalRows, numGlobalColumns, globalIndex( 1 ) } );
  int const width = static_cast< int >( std::log10( maximumDimension ) ) + 1;
  for( localIndex localRow = 0; localRow < localMatrix.numRows(); ++localRow )
  {
    arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( localRow );
    arraySlice1d< real64 const > const values = localMatrix.getEntries( localRow );
    for( localIndex entry = 0; entry < columns.size(); ++entry )
    {
      GEOS_FMT_TO( line, sizeof( line ), "{1:>{0}} {2:>{0}} {3:>24.16e}\n", width,
                   firstLocalRow + localRow + 1, columns[entry] + 1, values[entry] );
      output << line;
    }
  }
  output.close();

  array1d< globalIndex > localMetadata( 3 );
  localMetadata[0] = firstLocalRow;
  localMetadata[1] = numLocalRows;
  localMetadata[2] = localNonzeros;
  array1d< globalIndex > allMetadata;
  MpiWrapper::allGather( localMetadata.toViewConst(), allMetadata, comm );
  MpiWrapper::barrier( comm );

  if( rank == 0 )
  {
    string const manifestFilename = filenamePrefix + "_manifest.json";
    std::ofstream manifest( manifestFilename );
    GEOS_ERROR_IF( !manifest,
                   GEOS_FMT( "Unable to open distributed matrix manifest {}", manifestFilename ) );
    manifest << "{\n"
             << "  \"format\": \"distributed-matrix-market-v1\",\n"
             << "  \"global_rows\": " << numGlobalRows << ",\n"
             << "  \"global_columns\": " << numGlobalColumns << ",\n"
             << "  \"num_ranks\": " << numRanks << ",\n"
             << "  \"parts\": [\n";
    for( int partRank = 0; partRank < numRanks; ++partRank )
    {
      globalIndex const partFirstRow = allMetadata[3 * partRank];
      globalIndex const partNumRows = allMetadata[3 * partRank + 1];
      globalIndex const partNonzeros = allMetadata[3 * partRank + 2];
      manifest << "    { \"rank\": " << partRank
               << ", \"first_row\": " << partFirstRow
               << ", \"last_row_exclusive\": " << partFirstRow + partNumRows
               << ", \"nonzeros\": " << partNonzeros
               << ", \"file\": \"" << GEOS_FMT( "{}_rank_{:06}.mtx", filenamePrefix, partRank )
               << "\" }" << ( partRank + 1 == numRanks ? "\n" : ",\n" );
    }
    manifest << "  ]\n}\n";
  }
  MpiWrapper::barrier( comm );
}

} // namespace geos
