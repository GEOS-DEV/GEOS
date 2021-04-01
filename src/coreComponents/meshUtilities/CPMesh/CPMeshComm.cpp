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
 * @file CPMeshComm.cpp
 */

#include "CPMeshComm.hpp"

#include "mpiCommunications/MpiWrapper.hpp"

namespace geosx
{

namespace CPMesh
{

CPMeshComm::CPMeshComm( string const & name )
  :
  m_mainRank( 0 ),
  m_myRank( -1 ),
  m_size( 0 ),
  m_meshName( name )
{}

void CPMeshComm::setupMPIPartition( CPMeshDimensions & meshDims )
{
  m_myRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  m_size = MpiWrapper::commSize( MPI_COMM_GEOSX );

  if( m_size == 1 ) // serial run
  {
    meshDims.definePartitionBoundaries( 0, 0, meshDims.nX()-1, meshDims.nY()-1 );
    meshDims.definePartitionOverlaps( 0, 0, 0, 0 );
  }
  else // parallel run
  {
    broadcastDomainDimensions( meshDims );
    scatterPartitionBoundaries( meshDims );
  }
}

void CPMeshComm::broadcastDomainDimensions( CPMeshDimensions & meshDims ) const
{
  localIndex nCellsInEachDirection[3];

  // the main rank has read the top of the file and knows the domain dimensions
  if( m_myRank == m_mainRank )
  {
    nCellsInEachDirection[0] = meshDims.nX();
    nCellsInEachDirection[1] = meshDims.nY();
    nCellsInEachDirection[2] = meshDims.nZ();
  }
  localIndex const sendRecvSize = 3;
  MpiWrapper::bcast( nCellsInEachDirection, sendRecvSize, m_mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();

  // the other ranks have received the data and can set these values in CPMeshData
  if( m_myRank != m_mainRank )
  {
    meshDims.defineDomainDimensions( nCellsInEachDirection[0],
                                     nCellsInEachDirection[1],
                                     nCellsInEachDirection[2] );
  }
}

void CPMeshComm::scatterPartitionBoundaries( CPMeshDimensions & meshDims )
{
  localIndex const sizePerRank = 4;
  array1d< localIndex > boundaries( sizePerRank*m_size );
  array1d< localIndex > overlaps( sizePerRank*m_size );
  array1d< localIndex > neighbors( sizePerRank*m_size );

  if( m_myRank == m_mainRank )
  {
    // Temporary code to facilitate debugging: 1D partitioning along the x-direction.
    // Ultimately here we will partition in both X and Y, taking into account active cells in each partition.
    // Once I am sure that everything works on large test case, I will rewrite this entirely
    localIndex const spacing = floor( static_cast< real64 >( meshDims.nX() )
                                      / static_cast< real64 >( m_size ) );
    for( localIndex iRank = 0; iRank < m_size; ++iRank )
    {
      boundaries( iRank*sizePerRank )   = iRank * spacing;
      boundaries( iRank*sizePerRank+1 ) = 0;
      boundaries( iRank*sizePerRank+2 ) = ( iRank < m_size-1 ) ? ( iRank+1 ) * spacing - 1 : meshDims.nX()-1;
      boundaries( iRank*sizePerRank+3 ) = meshDims.nY()-1;

      overlaps( iRank*sizePerRank ) = ( boundaries( iRank*sizePerRank ) == 0 ) ? 0 : 1;
      overlaps( iRank*sizePerRank+1 ) = ( boundaries( iRank*sizePerRank+1 ) == 0 ) ? 0 : 1;
      overlaps( iRank*sizePerRank+2 ) = ( boundaries( iRank*sizePerRank+2 ) == meshDims.nX()-1 ) ? 0 : 1;
      overlaps( iRank*sizePerRank+3 ) = ( boundaries( iRank*sizePerRank+3 ) == meshDims.nY()-1 ) ? 0 : 1;

      neighbors( iRank*sizePerRank )   = ( iRank > 0 ) ? iRank - 1 : -1;
      neighbors( iRank*sizePerRank+1 ) = -1;
      neighbors( iRank*sizePerRank+2 ) = ( iRank < m_size-1 ) ? iRank + 1 : -1;
      neighbors( iRank*sizePerRank+3 ) = -1;
    }
  }

  // Step 1: communicate partition boundaries and set them
  array1d< localIndex > myBoundaries( sizePerRank );
  MpiWrapper::scatter( boundaries.toViewConst(), myBoundaries.toView(), m_mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();
  meshDims.definePartitionBoundaries( myBoundaries( 0 ), myBoundaries( 1 ),
                                      myBoundaries( 2 ), myBoundaries( 3 ) );

  // Step 2: communicate overlap between partitions and set them
  array1d< localIndex > myOverlaps( sizePerRank );
  MpiWrapper::scatter( overlaps.toViewConst(), myOverlaps.toView(), m_mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();
  meshDims.definePartitionOverlaps( myOverlaps( 0 ), myOverlaps( 1 ),
                                    myOverlaps( 2 ), myOverlaps( 3 ) );

  // Step 3: communicate neighbor lists
  array1d< localIndex > myNeighbors( sizePerRank );
  MpiWrapper::scatter( neighbors.toViewConst(), myNeighbors.toView(), m_mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();
  for( localIndex i = 0; i < sizePerRank; ++i )
  {
    if( myNeighbors( i ) != -1 )
    {
      m_neighborsList.insert( myNeighbors( i ) );
    }
  }
}

REGISTER_CATALOG_ENTRY( CPMeshComm, CPMeshComm, string const & )

} // namespace CPMesh

} // end namespace geosx
