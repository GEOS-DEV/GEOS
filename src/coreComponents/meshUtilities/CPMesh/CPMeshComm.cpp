/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CPMeshBuilder.cpp
 */

#include "CPMeshBuilder.hpp"

#include "meshUtilities/CPMesh/CPMeshData.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

namespace geosx
{

namespace CPMesh
{

CPMeshComm::CPMeshComm( string const & name )
  : m_meshName( name ) {}

void CPMeshComm::setupMPIPartition( CPMeshData & cPMeshData ) const
{
  localIndex const myRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  localIndex const size = MpiWrapper::commSize( MPI_COMM_GEOSX );

  if( size == 1 ) // serial run
  {
    cPMeshData.definePartitionBoundaries( 0, 0, cPMeshData.nX()-1, cPMeshData.nY()-1 );
    cPMeshData.definePartitionOverlaps( 0, 0, 0, 0 );

  }
  else // parallel run
  {
    broadcastDomainBoundaries( myRank, cPMeshData );
    scatterPartitionBoundaries( myRank, size, cPMeshData );
  }
  fillOwnershipArrays( cPMeshData );
}

void CPMeshComm::broadcastDomainBoundaries( localIndex const myRank,
                                            CPMeshData & cPMeshData ) const
{
  localIndex nCellsInEachDirection[3];
  localIndex const mainRank = 0;
  if( myRank == mainRank )
  {
    nCellsInEachDirection[0] = cPMeshData.nX();
    nCellsInEachDirection[1] = cPMeshData.nY();
    nCellsInEachDirection[2] = cPMeshData.nZ();
  }
  localIndex const sendRecvSize = 3;
  MpiWrapper::bcast( nCellsInEachDirection, sendRecvSize, mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();

  if( myRank != mainRank )
  {
    cPMeshData.defineDomainBoundaries( nCellsInEachDirection[0],
                                       nCellsInEachDirection[1],
                                       nCellsInEachDirection[2] );
  }
}

void CPMeshComm::scatterPartitionBoundaries( localIndex const myRank,
                                             localIndex const size,
                                             CPMeshData & cPMeshData ) const
{
  localIndex const mainRank = 0;
  localIndex const sizePerRank = 4;
  array1d< localIndex > boundaries( sizePerRank*size );
  array1d< localIndex > overlaps( sizePerRank*size );
  if( myRank == mainRank )
  {
    // TODO: do something smarter here
    // Definitely no what we ultimately want, but makes debugging easier
    // Should have here partitioning in X and Y, taking into account active cells in each partition
    localIndex const spacing = floor( static_cast< real64 >( cPMeshData.nX() )
                                      / static_cast< real64 >( size ) );
    for( localIndex iRank = 0; iRank < size; ++iRank )
    {
      boundaries( iRank*sizePerRank )   = iRank * spacing;
      boundaries( iRank*sizePerRank+1 ) = 0;
      boundaries( iRank*sizePerRank+2 ) = ( iRank < size-1 ) ? (iRank+1) * spacing - 1 : cPMeshData.nZ()-1;
      boundaries( iRank*sizePerRank+3 ) = cPMeshData.nY()-1;

      overlaps( iRank*sizePerRank ) = ( boundaries( iRank*sizePerRank ) == 0 ) ? 0 : 1;
      overlaps( iRank*sizePerRank+1 ) = ( boundaries( iRank*sizePerRank+1 ) == 0 ) ? 0 : 1;
      overlaps( iRank*sizePerRank+2 ) = ( boundaries( iRank*sizePerRank+2 ) == cPMeshData.nX()-1 ) ? 0 : 1;
      overlaps( iRank*sizePerRank+3 ) = ( boundaries( iRank*sizePerRank+3 ) == cPMeshData.nY()-1 ) ? 0 : 1;
    }
  }

  array1d< localIndex > myBoundaries( sizePerRank );
  array1d< localIndex > myOverlaps( sizePerRank );
  // TODO: move to mpiCommunications/MpiWrapper.hpp
  MPI_Scatter( boundaries.data(),
               sizePerRank,
               MpiWrapper::getMpiType< localIndex >(),
               myBoundaries.data(),
               sizePerRank,
               MpiWrapper::getMpiType< localIndex >(),
               mainRank,
               MPI_COMM_GEOSX );
  MpiWrapper::barrier();
  // define the iMin, iMax, jMin, jMax for each rank
  cPMeshData.definePartitionBoundaries( myBoundaries( 0 ), myBoundaries( 1 ),
                                        myBoundaries( 2 ), myBoundaries( 3 ) );

  // TODO: move to mpiCommunications/MpiWrapper.hpp
  MPI_Scatter( overlaps.data(),
               sizePerRank,
               MpiWrapper::getMpiType< localIndex >(),
               myOverlaps.data(),
               sizePerRank,
               MpiWrapper::getMpiType< localIndex >(),
               mainRank,
               MPI_COMM_GEOSX );
  MpiWrapper::barrier();
  // define the iMinOverlap, iMaxOverlap, jMinOverlap, jMaxOverlap for each rank
  cPMeshData.definePartitionOverlaps( myOverlaps( 0 ), myOverlaps( 1 ),
                                      myOverlaps( 2 ), myOverlaps( 3 ) );

  // std::cout << "myRank = " << myRank << std::endl;
  // std::cout << "iMin = " << myBoundaries( 0 ) << " iMax = " << myBoundaries( 2 )
  //           << " jMin = " << myBoundaries( 1 ) << " jMax = " << myBoundaries( 3 ) << std::endl
  //      << "iMinOverlap = " << myOverlaps( 0 ) << " iMaxOverlap = " << myOverlaps( 2 )
  //           << " jMinOverlap = " << myOverlaps( 1 ) << " jMaxOverlap = " << myOverlaps( 3 )
  //           << std::endl;
}

void CPMeshComm::fillOwnershipArrays( CPMeshData & cPMeshData ) const
{
  localIndex const nXLocal = cPMeshData.nXLocal();
  localIndex const nYLocal = cPMeshData.nYLocal();
  localIndex const nZLocal = cPMeshData.nZLocal();

  localIndex const iMinLocal = cPMeshData.iMinLocal();
  localIndex const jMinLocal = cPMeshData.jMinLocal();

  localIndex const iMinOwned = cPMeshData.iMinOwned();
  localIndex const jMinOwned = cPMeshData.jMinOwned();
  localIndex const iMaxOwned = cPMeshData.iMaxOwned();
  localIndex const jMaxOwned = cPMeshData.jMaxOwned();

  array1d< bool > & localCellIsOwned = cPMeshData.localCellIsOwned();
  localCellIsOwned.resize( nXLocal*nYLocal*nZLocal );
  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        localIndex const iLocalCell = k*nXLocal*nYLocal + j*nXLocal + i;
        bool const isOwned = ( i+iMinLocal >= iMinOwned && i+iMinLocal < iMaxOwned
                               && j+jMinLocal >= jMinOwned && j+jMinLocal < jMaxOwned );
        localCellIsOwned( iLocalCell ) = isOwned;
      }
    }
  }
}

localIndex CPMeshComm::gatherLocalCountForEachRank( localIndex const countOnMyRank ) const
{
  localIndex const myRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  localIndex const size = MpiWrapper::commSize( MPI_COMM_GEOSX );
  if( size == 1 )
  {
    return 0;
  }

  array1d< localIndex > countAll( size );
  MpiWrapper::allGather( countOnMyRank, countAll, MPI_COMM_GEOSX );

  localIndex offset = 0;
  for( localIndex iRank = 0; iRank < myRank; ++iRank )
  {
    offset += countAll( iRank );
  }
  return offset;
}

void CPMeshComm::synchronizeBoundaryVertices( CPMeshData & m_cPMeshData ) const
{
  GEOSX_UNUSED_VAR( m_cPMeshData );
}

REGISTER_CATALOG_ENTRY( CPMeshComm, CPMeshComm, string const & )

} // namespace CPMesh

} // end namespace geosx
