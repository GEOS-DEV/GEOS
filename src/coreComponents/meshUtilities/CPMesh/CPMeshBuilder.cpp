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

#include "codingUtilities/Utilities.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "meshUtilities/CPMesh/utilities/GeometryUtilities.hpp"
#include "meshUtilities/CPMesh/utilities/OutputUtilities.hpp"

namespace geosx
{

namespace CPMesh
{

using namespace GeometryUtilities;
using namespace OutputUtilities;

CPMeshBuilder::CPMeshBuilder( string const & name )
  :
  m_cPMeshParser( name ),
  m_cPMeshData( name ),
  m_cPMeshComm( name ),
  m_meshName( name )
{}

void CPMeshBuilder::buildMesh( Path const & filePath )
{
  localIndex const iRank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  if( iRank == 0 )
  {
    localIndex nX = 0;
    localIndex nY = 0;
    localIndex nZ = 0;

    // read SPECGRID/DIMENS to get the number of cells in each direction
    m_cPMeshParser.readNumberOfCells( filePath, nX, nY, nZ );
    m_cPMeshData.defineDomainBoundaries( nX, nY, nZ );

    // TODO: read ACTNUM for the inertial partitioning
  }

  // rank 0 partitions the mesh and sends partitions to other ranks
  // rank 0 sets up the communication pattern
  m_cPMeshComm.setupMPIPartition( m_cPMeshData );

  // each rank reads its part of the mesh
  m_cPMeshParser.readMesh( filePath, m_cPMeshData );

  // post-process the mesh to make it conforming
  postProcessMesh();
}

void CPMeshBuilder::postProcessMesh()
{
  // for each active cell, compute the position of the eight CPG vertices
  // at this point, CPG vertices have not been filtered yet
  buildHexahedralMesh();

  // eliminate duplicates from cPVertices
  filterVertices();

  // communicate with other partitions to assign global (active) cell numbers
  assignGlobalActiveCellIndices();

  // communicate with other partitions to assign global vertex numbers
  assignGlobalVertexIndices();

  // match faces, deal with the non-conforming case
  buildFaces();

  // for now, do some debugging
  outputDebugVTKFile( m_cPMeshData );
}

void CPMeshBuilder::buildHexahedralMesh()
{
  localIndex const nXLocal = m_cPMeshData.nXLocal();
  localIndex const nYLocal = m_cPMeshData.nYLocal();
  localIndex const nZLocal = m_cPMeshData.nZLocal();
  localIndex const nLocalCells = nXLocal*nYLocal*nZLocal;

  array1d< localIndex > const & actnum = m_cPMeshData.actnum();
  array1d< real64 > const & coord = m_cPMeshData.coord();
  array1d< real64 > const & zcorn = m_cPMeshData.zcorn();

  array2d< real64 > & cPVertices = m_cPMeshData.cPVertices();
  cPVertices.resizeDimension< 0, 1 >( 8*nLocalCells, 3 );

  localIndex & nOwnedActiveCells = m_cPMeshData.nOwnedActiveCells();
  localIndex & nLocalActiveCells = m_cPMeshData.nLocalActiveCells();
  array1d< bool > const & localCellIsOwned = m_cPMeshData.localCellIsOwned();
  array1d< localIndex > & localActiveCellToLocalCell = m_cPMeshData.localActiveCellToLocalCell();
  array1d< localIndex > & localCellToLocalCPVertices = m_cPMeshData.localCellToLocalCPVertices();
  localActiveCellToLocalCell.reserve( nLocalCells );
  localCellToLocalCPVertices.resize( nLocalCells );

  array1d< real64 > xPos( 8 );
  array1d< real64 > yPos( 8 );
  array1d< real64 > zPos( 8 );

  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        bool activeAndValid = false;

        // compute explicit local index using the ijk structure
        localIndex const iLocalCell = k*nXLocal*nYLocal + j*nXLocal + i;

        if( actnum( iLocalCell ) == 1 )
        {
          // compute the positions of the eight vertices
          activeAndValid = processActiveHexahedron( i, j, k,
                                                    nXLocal, nYLocal,
                                                    coord, zcorn,
                                                    xPos, yPos, zPos );
        }

        if( activeAndValid )
        {
          // save the local indices
          if( localCellIsOwned( iLocalCell ) )
          {
            nOwnedActiveCells++;
          }
          nLocalActiveCells++;
          localActiveCellToLocalCell.emplace_back( iLocalCell );

          // construct the map from local cell to first CP vertex of the cell
          localIndex const iFirstVertexLocal = 8 * iLocalCell;
          localCellToLocalCPVertices( iLocalCell ) = iFirstVertexLocal;

          // save the position of the eight vertices
          localIndex const order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
          for( localIndex pos = 0; pos < 8; ++pos )
          {
            cPVertices( iFirstVertexLocal + pos, 0 ) = xPos( order[pos] );
            cPVertices( iFirstVertexLocal + pos, 1 ) = yPos( order[pos] );
            cPVertices( iFirstVertexLocal + pos, 2 ) = zPos( order[pos] );
          }
        }
      }
    }
  }
}

bool CPMeshBuilder::processActiveHexahedron( localIndex const i, localIndex const j, localIndex const k,
                                             localIndex const nXLocal, localIndex const nYLocal,
                                             array1d< real64 > const & coord,
                                             array1d< real64 > const & zcorn,
                                             array1d< real64 > & xPos,
                                             array1d< real64 > & yPos,
                                             array1d< real64 > & zPos ) const
{
  localIndex const i_xmLocal = k*8*nXLocal*nYLocal + j*4*nXLocal + 2*i;
  localIndex const i_xpLocal = i_xmLocal + 2*nXLocal;

  zPos( 0 ) = zcorn( i_xmLocal );
  zPos( 1 ) = zcorn( i_xmLocal + 1 );
  zPos( 2 ) = zcorn( i_xpLocal );
  zPos( 3 ) = zcorn( i_xpLocal + 1 );
  zPos( 4 ) = zcorn( i_xmLocal + 4*nXLocal*nYLocal );
  zPos( 5 ) = zcorn( i_xmLocal + 4*nXLocal*nYLocal + 1 );
  zPos( 6 ) = zcorn( i_xpLocal + 4*nXLocal*nYLocal );
  zPos( 7 ) = zcorn( i_xpLocal + 4*nXLocal*nYLocal + 1 );

  bool const hexaIsFlat = isZero( zPos( 0 )-zPos( 4 )
                                  + zPos( 1 )-zPos( 5 )
                                  + zPos( 2 )-zPos( 6 )
                                  + zPos( 3 )-zPos( 7 ) );

  if( hexaIsFlat )
  {
    return false;
  }

  // the following code comes from PAMELA
  // TODO: if possible, remove code duplication

  // first pillar
  localIndex const indFirstPillarLocal = j*(nXLocal+1)+i;
  localIndex const ip0 = 6*indFirstPillarLocal - 1;
  real64 const denomFirstPillar = coord( ip0 + 6 ) - coord( ip0 + 3 );

  real64 const slopePos0 = isZero( denomFirstPillar )
    ? 1.0
    : ( zPos( 0 )-coord( ip0 + 3 ) ) / denomFirstPillar;
  xPos( 0 ) = slopePos0 * (coord( ip0 + 4 ) - coord( ip0 + 1 )) + coord( ip0 + 1 );
  yPos( 0 ) = slopePos0 * (coord( ip0 + 5 ) - coord( ip0 + 2 )) + coord( ip0 + 2 );

  real64 const slopePos4 = isZero( denomFirstPillar )
    ? 1.0
    : ( zPos( 4 )-coord( ip0 + 3 ) ) / denomFirstPillar;
  xPos( 4 ) = slopePos4 * (coord( ip0 + 4 ) - coord( ip0 + 1 )) + coord( ip0 + 1 );
  yPos( 4 ) = slopePos4 * (coord( ip0 + 5 ) - coord( ip0 + 2 )) + coord( ip0 + 2 );

  // second pillar
  localIndex const indSecondPillarLocal = (nXLocal+1)*j + i+1;
  localIndex const ip1 = 6*indSecondPillarLocal - 1;
  real64 const denomSecondPillar = coord( ip1 + 6 ) - coord( ip1 + 3 );

  real64 const slopePos1 = isZero( denomSecondPillar )
    ? 1.0
    : ( zPos( 1 )-coord( ip1 + 3 ) ) / denomSecondPillar;
  xPos( 1 ) = slopePos1 * (coord( ip1 + 4 ) - coord( ip1 + 1 )) + coord( ip1 + 1 );
  yPos( 1 ) = slopePos1 * (coord( ip1 + 5 ) - coord( ip1 + 2 )) + coord( ip1 + 2 );

  real64 const slopePos5 = isZero( denomSecondPillar )
    ? 1.0
    : ( zPos( 5 ) - coord( ip1 + 3 ) ) / denomSecondPillar;
  xPos( 5 ) = slopePos5 * (coord( ip1 + 4 ) - coord( ip1 + 1 )) + coord( ip1 + 1 );
  yPos( 5 ) = slopePos5 * (coord( ip1 + 5 ) - coord( ip1 + 2 )) + coord( ip1 + 2 );

  // third pillar
  localIndex const indThirdPillarLocal = (nXLocal+1)*(j+1) + i;
  localIndex const ip2 = 6*indThirdPillarLocal - 1;
  real64 const denomThirdPillar = coord( ip2 + 6 ) - coord( ip2 + 3 );

  real64 const slopePos2 = isZero( denomThirdPillar )
    ? 1.0
    : ( zPos( 2 )-coord( ip2 + 3 ) ) / denomThirdPillar;
  xPos( 2 ) = slopePos2 * (coord( ip2 + 4 ) - coord( ip2 + 1 )) + coord( ip2 + 1 );
  yPos( 2 ) = slopePos2 * (coord( ip2 + 5 ) - coord( ip2 + 2 )) + coord( ip2 + 2 );

  real64 const slopePos6 = isZero( denomThirdPillar )
    ? 1.0
    : ( zPos( 6 )-coord( ip2 + 3 ) ) / denomThirdPillar;
  xPos( 6 ) = slopePos6 * (coord( ip2 + 4 ) - coord( ip2 + 1 )) + coord( ip2 + 1 );
  yPos( 6 ) = slopePos6 * (coord( ip2 + 5 ) - coord( ip2 + 2 )) + coord( ip2 + 2 );

  // fourth pillar
  localIndex const indFourthPillarLocal = (nXLocal+1)*(j+1) + i+1;
  localIndex const ip3 = 6*indFourthPillarLocal - 1;
  real64 const denomFourthPillar = coord( ip3 + 6 ) - coord( ip3 + 3 );

  real64 const slopePos3 = isZero( denomFourthPillar )
    ? 1.0
    : ( zPos( 3 )-coord( ip3 + 3 ) ) / denomFourthPillar;
  xPos( 3 ) = slopePos3 * (coord( ip3 + 4 ) - coord( ip3 + 1 )) + coord( ip3 + 1 );
  yPos( 3 ) = slopePos3 * (coord( ip3 + 5 ) - coord( ip3 + 2 )) + coord( ip3 + 2 );

  real64 const slopePos7 = isZero( denomFourthPillar )
    ? 1.0
    : ( zPos( 7 )-coord( ip3 + 3 ) ) / denomFourthPillar;
  xPos( 7 ) = slopePos7 * (coord( ip3 + 4 ) - coord( ip3 + 1 )) + coord( ip3 + 1 );
  yPos( 7 ) = slopePos7 * (coord( ip3 + 5 ) - coord( ip3 + 2 )) + coord( ip3 + 2 );

  return true;
}

void CPMeshBuilder::filterVertices()
{
  array2d< real64 > const & cPVertices = m_cPMeshData.cPVertices();
  array2d< real64 > & vertices = m_cPMeshData.vertices();
  vertices.reserve( 3*cPVertices.size() );

  array1d< bool > const & cPVertexIsOwned = m_cPMeshData.cPVertexIsOwned();
  array1d< localIndex > & localCPVertexToLocalVertex = m_cPMeshData.localCPVertexToLocalVertex();
  localCPVertexToLocalVertex.resize( cPVertices.size() );

  // First step: filter the unique vertices using a set

  std::set< Vertex, CompareVertices > uniqueVerticesHelper;
  std::set< Vertex, CompareVertices >::iterator it;

  localIndex & nOwnedVertices = m_cPMeshData.nOwnedVertices();
  // loop over of the CP vertices (for now, including those of inactive cells)
  for( localIndex iCPVertex = 0; iCPVertex < cPVertices.size( 0 ); ++iCPVertex )
  {
    Vertex v( cPVertices( iCPVertex, 0 ),
              cPVertices( iCPVertex, 1 ),
              cPVertices( iCPVertex, 2 ) );

    // check if this vertex has already been found and inserted in the set
    it = uniqueVerticesHelper.find( v );

    // if already found, copy the local index in map and move on
    if( it != uniqueVerticesHelper.end() )
    {
      Vertex const & existingVertex = *it;
      localCPVertexToLocalVertex( iCPVertex ) = existingVertex.m_localVertex;
    }
    // if not found yet, insert into the set
    else
    {
      v.m_localVertex = uniqueVerticesHelper.size();
      localCPVertexToLocalVertex( iCPVertex ) = v.m_localVertex;
      // if this vertex is owned by my rank,
      // increment the counter (later used to assign global indices)
      if( cPVertexIsOwned( iCPVertex ) )
      {
        nOwnedVertices++;
      }
      uniqueVerticesHelper.insert( v );
    }
  }

  // Second step: move the data from the set to a array2d, because we don't want the set anymore
  vertices.resizeDimension< 0, 1 >( uniqueVerticesHelper.size(), 3 );
  localIndex & nLocalVertices = m_cPMeshData.nLocalVertices();
  nLocalVertices = uniqueVerticesHelper.size();

  std::for_each( uniqueVerticesHelper.begin(), uniqueVerticesHelper.end(), [&vertices]( Vertex const & v )
  {
    vertices( v.m_localVertex, 0 ) = v.m_x;
    vertices( v.m_localVertex, 1 ) = v.m_y;
    vertices( v.m_localVertex, 2 ) = v.m_z;
  } );
}

void CPMeshBuilder::assignGlobalActiveCellIndices()
{
  // First step: gather all the nLocalActiveCells from all the ranks
  localIndex const nOwnedActiveCells = m_cPMeshData.nOwnedActiveCells();
  localIndex const offset = m_cPMeshComm.gatherLocalCountForEachRank( nOwnedActiveCells );

  // Second step: using the offset for my rank, construct the local to global map
  localIndex const nLocalActiveCells = m_cPMeshData.nLocalActiveCells();
  array1d< bool > const & localCellIsOwned = m_cPMeshData.localCellIsOwned();
  array1d< localIndex > const & localActiveCellToLocalCell = m_cPMeshData.localActiveCellToLocalCell();
  array1d< globalIndex > & localActiveCellToGlobalActiveCell = m_cPMeshData.localActiveCellToGlobalActiveCell();
  localActiveCellToGlobalActiveCell.resize( nOwnedActiveCells );

  localIndex counter = 0;
  for( localIndex iActiveCell = 0; iActiveCell < nLocalActiveCells; ++iActiveCell )
  {
    localIndex const iLocalCell = localActiveCellToLocalCell( iActiveCell );
    if( localCellIsOwned( iLocalCell ) )
    {
      localActiveCellToGlobalActiveCell( iActiveCell ) = offset + counter;
      counter++;
    }
    else
    {
      localActiveCellToGlobalActiveCell( iActiveCell ) = -1;
    }
  }

  GEOSX_THROW_IF( counter != nOwnedActiveCells,
                  "Partitioning error for the cells",
                  InputError );

  // for now, we don't care about global indices for the overlaps...
}

void CPMeshBuilder::assignGlobalVertexIndices()
{
  // First step: gather all the nOwnedVertices from all the ranks
  localIndex & nOwnedVertices = m_cPMeshData.nOwnedVertices();
  localIndex const offset = m_cPMeshComm.gatherLocalCountForEachRank( nOwnedVertices );

  // Second step: using the offset for my rank, construct the local to global map ** for the owned vertices **
  localIndex const nLocalVertices = m_cPMeshData.nLocalVertices();
  array1d< bool > const & vertexIsOwned = m_cPMeshData.vertexIsOwned();
  array1d< globalIndex > & localVertexToGlobalVertex = m_cPMeshData.localVertexToGlobalVertex();
  localVertexToGlobalVertex.resize( nLocalVertices );

  localIndex counter = 0;
  for( localIndex iVertex = 0; iVertex < nLocalVertices; ++iVertex )
  {
    if( vertexIsOwned( iVertex ) )
    {
      localVertexToGlobalVertex( iVertex ) = offset + counter;
      counter++;
    }
    else
    {
      localVertexToGlobalVertex( iVertex ) = -1; // to make sure it crashes if I mess up
    }
  }
  GEOSX_THROW_IF( counter != nOwnedVertices,
                  "Partitioning error for the vertices",
                  InputError );

  // Third step: synchronize vertices shared between partitions
  m_cPMeshComm.synchronizeBoundaryVertices( m_cPMeshData );
}

void CPMeshBuilder::buildFaces()
{
  // TODO
}

REGISTER_CATALOG_ENTRY( CPMeshBuilder, CPMeshBuilder, string const & )

} // namespace CPMesh

} // end namespace geosx
