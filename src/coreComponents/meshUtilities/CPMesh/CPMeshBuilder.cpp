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
  m_meshParser( name ),
  m_meshComm( name ),
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
    m_meshParser.readNumberOfCells( filePath, nX, nY, nZ );
    m_meshDims.defineDomainDimensions( nX, nY, nZ );

    // TODO: read ACTNUM for the inertial partitioning
  }

  // rank 0 partitions the mesh and sends partitions to other ranks
  // rank 0 sets up the communication pattern
  m_meshComm.setupMPIPartition( m_meshDims );

  // each rank reads its part of the mesh
  m_meshParser.readMesh( filePath, m_meshDims );

  // post-process the mesh to make it conforming
  postProcessMesh();
}

void CPMeshBuilder::postProcessMesh()
{
  // for each active cell, compute the position of the eight CPG vertices
  // at this point, CPG vertices have not been filtered yet
  buildCornerPointCells();

  // eliminate duplicates from cPVertices
  filterVertices();

  // match faces, deal with the non-conforming case
  buildFaces();

  // for now, do some debugging
  outputDebugVTKFile( m_meshVertices, m_meshFaces, m_meshCells );
}

void CPMeshBuilder::buildCornerPointCells()
{
  // local and global dimensions
  localIndex const nX = m_meshDims.nX();
  localIndex const nY = m_meshDims.nY();

  localIndex const iMinLocal = m_meshDims.iMinLocal();
  localIndex const jMinLocal = m_meshDims.jMinLocal();
  localIndex const kMinLocal = 0;

  localIndex const iMinOverlap = m_meshDims.iMinOverlap();
  localIndex const jMinOverlap = m_meshDims.jMinOverlap();
  localIndex const iMaxOverlap = m_meshDims.iMaxOverlap();
  localIndex const jMaxOverlap = m_meshDims.jMaxOverlap();

  localIndex const nXLocal = m_meshDims.nXLocal();
  localIndex const nYLocal = m_meshDims.nYLocal();
  localIndex const nZLocal = m_meshDims.nZLocal();
  localIndex const nLocalCells = nXLocal*nYLocal*nZLocal;
  localIndex constexpr nCPVerticesPerCell = 8;
  localIndex const nCPVertices = nCPVerticesPerCell*nLocalCells;

  // corner-point mesh file information
  array1d< localIndex > const & actnum = m_meshParser.actnum();
  array1d< real64 > const & coord = m_meshParser.coord();
  array1d< real64 > const & zcorn = m_meshParser.zcorn();

  // vertex maps
  array2d< real64 > & cPVertices = m_meshVertices.m_cPVertices;
  array1d< globalIndex > & cPVertexToGlobalCPVertex = m_meshVertices.m_cPVertexToGlobalCPVertex;
  array1d< bool > & cPVertexIsInsidePartition = m_meshVertices.m_cPVertexIsInsidePartition;
  cPVertices.resizeDimension< 0, 1 >( nCPVertices, 3 );
  cPVertexToGlobalCPVertex.resize( nCPVertices );
  cPVertexIsInsidePartition.resize( nCPVertices );

  // cell maps
  array1d< localIndex > & activeCellInsidePartitionToActiveCell = m_meshCells.m_activeCellInsidePartitionToActiveCell;
  array1d< localIndex > & activeCellToCell = m_meshCells.m_activeCellToCell;
  array1d< globalIndex > & activeCellToGlobalCell = m_meshCells.m_activeCellToGlobalCell;
  array1d< localIndex > & cellToCPVertices = m_meshCells.m_cellToCPVertices;
  activeCellInsidePartitionToActiveCell.reserve( nLocalCells );
  activeCellToCell.reserve( nLocalCells );
  activeCellToGlobalCell.reserve( nLocalCells );
  cellToCPVertices.resize( nLocalCells );

  array1d< real64 > xPos( nCPVerticesPerCell );
  array1d< real64 > yPos( nCPVerticesPerCell );
  array1d< real64 > zPos( nCPVerticesPerCell );

  array1d< bool > cPVertexIsInside( nCPVerticesPerCell );

  // loop over all cells in the MPI domain, including overlaps
  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        bool activeAndValid = false;

        // compute explicit local and global indices using the ijk structure
        localIndex const iLocalCell = k*nXLocal*nYLocal + j*nXLocal + i;
        globalIndex const iGlobalCell = (k+kMinLocal)*nX*nY + (j+jMinLocal)*nX + (i+iMinLocal);

        if( actnum( iLocalCell ) == 1 )
        {
          // compute the positions of the eight vertices
          activeAndValid = processActiveHexahedron( i, j, k,
                                                    nXLocal, nYLocal,
                                                    iMinOverlap, iMaxOverlap,
                                                    jMinOverlap, jMaxOverlap,
                                                    coord, zcorn,
                                                    xPos, yPos, zPos,
                                                    cPVertexIsInside );
        }

        if( activeAndValid )
        {

          // cell maps

          // assign local and global indices
          activeCellToCell.emplace_back( iLocalCell );
          activeCellToGlobalCell.emplace_back( iGlobalCell );
          // decide flag specifying whether the cell is in the overlap or not
          bool const isInside = !( (i == 0 && iMinOverlap == 1) || (j == 0 && jMinOverlap == 1) ||
                                   (i+1 == nXLocal && iMaxOverlap == 1) || (j+1 == nYLocal && jMaxOverlap == 1) );
          if( isInside )
          {
            activeCellInsidePartitionToActiveCell.emplace_back( activeCellToCell.size()-1 );
          }

          // vertex maps

          // construct the map from local cell to first CP vertex of the cell
          localIndex const iFirstVertexLocal = nCPVerticesPerCell * iLocalCell;
          globalIndex const iFirstVertexGlobal = nCPVerticesPerCell * iGlobalCell;
          cellToCPVertices( iLocalCell ) = iFirstVertexLocal;

          // save the position of the eight vertices, and assign global CP vertex indices
          localIndex const order[nCPVerticesPerCell] = { 0, 1, 3, 2, 4, 5, 7, 6 };
          for( localIndex pos = 0; pos < nCPVerticesPerCell; ++pos )
          {
            cPVertices( iFirstVertexLocal + pos, 0 ) = xPos( order[pos] );
            cPVertices( iFirstVertexLocal + pos, 1 ) = yPos( order[pos] );
            cPVertices( iFirstVertexLocal + pos, 2 ) = zPos( order[pos] );
            cPVertexToGlobalCPVertex( iFirstVertexLocal + pos ) = iFirstVertexGlobal + order[pos];
            cPVertexIsInsidePartition( iFirstVertexLocal + pos ) = cPVertexIsInside( order[pos] );
          }
        }
      }
    }
  }
}

bool CPMeshBuilder::processActiveHexahedron( localIndex const i, localIndex const j, localIndex const k,
                                             localIndex const nXLocal, localIndex const nYLocal,
                                             localIndex const iMinOverlap, localIndex const iMaxOverlap,
                                             localIndex const jMinOverlap, localIndex const jMaxOverlap,
                                             array1d< real64 > const & coord,
                                             array1d< real64 > const & zcorn,
                                             array1d< real64 > & xPos,
                                             array1d< real64 > & yPos,
                                             array1d< real64 > & zPos,
                                             array1d< bool > & cPVertexIsInside ) const
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
  bool const firstPillarIsInside = !( (i == 0 && iMinOverlap == 1) || (j == 0 && jMinOverlap == 1) );

  real64 const slopePos0 = isZero( denomFirstPillar )
    ? 1.0
    : ( zPos( 0 )-coord( ip0 + 3 ) ) / denomFirstPillar;
  xPos( 0 ) = slopePos0 * (coord( ip0 + 4 ) - coord( ip0 + 1 )) + coord( ip0 + 1 );
  yPos( 0 ) = slopePos0 * (coord( ip0 + 5 ) - coord( ip0 + 2 )) + coord( ip0 + 2 );
  cPVertexIsInside( 0 ) = firstPillarIsInside;

  real64 const slopePos4 = isZero( denomFirstPillar )
    ? 1.0
    : ( zPos( 4 )-coord( ip0 + 3 ) ) / denomFirstPillar;
  xPos( 4 ) = slopePos4 * (coord( ip0 + 4 ) - coord( ip0 + 1 )) + coord( ip0 + 1 );
  yPos( 4 ) = slopePos4 * (coord( ip0 + 5 ) - coord( ip0 + 2 )) + coord( ip0 + 2 );
  cPVertexIsInside( 4 ) = firstPillarIsInside;

  // second pillar
  localIndex const indSecondPillarLocal = (nXLocal+1)*j + i+1;
  localIndex const ip1 = 6*indSecondPillarLocal - 1;
  real64 const denomSecondPillar = coord( ip1 + 6 ) - coord( ip1 + 3 );
  bool const secondPillarIsInside = !( (i+1 == nXLocal && iMaxOverlap == 1) || (j == 0 && jMinOverlap == 1) );

  real64 const slopePos1 = isZero( denomSecondPillar )
    ? 1.0
    : ( zPos( 1 )-coord( ip1 + 3 ) ) / denomSecondPillar;
  xPos( 1 ) = slopePos1 * (coord( ip1 + 4 ) - coord( ip1 + 1 )) + coord( ip1 + 1 );
  yPos( 1 ) = slopePos1 * (coord( ip1 + 5 ) - coord( ip1 + 2 )) + coord( ip1 + 2 );
  cPVertexIsInside( 1 ) = secondPillarIsInside;

  real64 const slopePos5 = isZero( denomSecondPillar )
    ? 1.0
    : ( zPos( 5 )-coord( ip1 + 3 ) ) / denomSecondPillar;
  xPos( 5 ) = slopePos5 * (coord( ip1 + 4 ) - coord( ip1 + 1 )) + coord( ip1 + 1 );
  yPos( 5 ) = slopePos5 * (coord( ip1 + 5 ) - coord( ip1 + 2 )) + coord( ip1 + 2 );
  cPVertexIsInside( 5 ) = secondPillarIsInside;

  // third pillar
  localIndex const indThirdPillarLocal = (nXLocal+1)*(j+1) + i;
  localIndex const ip2 = 6*indThirdPillarLocal - 1;
  real64 const denomThirdPillar = coord( ip2 + 6 ) - coord( ip2 + 3 );
  bool const thirdPillarIsInside = !( (i == 0 && iMinOverlap == 1) || (j+1 == nYLocal && jMaxOverlap == 1) );

  real64 const slopePos2 = isZero( denomThirdPillar )
    ? 1.0
    : ( zPos( 2 )-coord( ip2 + 3 ) ) / denomThirdPillar;
  xPos( 2 ) = slopePos2 * (coord( ip2 + 4 ) - coord( ip2 + 1 )) + coord( ip2 + 1 );
  yPos( 2 ) = slopePos2 * (coord( ip2 + 5 ) - coord( ip2 + 2 )) + coord( ip2 + 2 );
  cPVertexIsInside( 2 ) = thirdPillarIsInside;

  real64 const slopePos6 = isZero( denomThirdPillar )
    ? 1.0
    : ( zPos( 6 )-coord( ip2 + 3 ) ) / denomThirdPillar;
  xPos( 6 ) = slopePos6 * (coord( ip2 + 4 ) - coord( ip2 + 1 )) + coord( ip2 + 1 );
  yPos( 6 ) = slopePos6 * (coord( ip2 + 5 ) - coord( ip2 + 2 )) + coord( ip2 + 2 );
  cPVertexIsInside( 6 ) = thirdPillarIsInside;

  // fourth pillar
  localIndex const indFourthPillarLocal = (nXLocal+1)*(j+1) + i+1;
  localIndex const ip3 = 6*indFourthPillarLocal - 1;
  real64 const denomFourthPillar = coord( ip3 + 6 ) - coord( ip3 + 3 );
  bool const fourthPillarIsInside = !( (i+1 == nXLocal && iMaxOverlap == 1) || (j+1 == nYLocal && jMaxOverlap == 1) );

  real64 const slopePos3 = isZero( denomFourthPillar )
    ? 1.0
    : ( zPos( 3 )-coord( ip3 + 3 ) ) / denomFourthPillar;
  xPos( 3 ) = slopePos3 * (coord( ip3 + 4 ) - coord( ip3 + 1 )) + coord( ip3 + 1 );
  yPos( 3 ) = slopePos3 * (coord( ip3 + 5 ) - coord( ip3 + 2 )) + coord( ip3 + 2 );
  cPVertexIsInside( 3 ) = fourthPillarIsInside;

  real64 const slopePos7 = isZero( denomFourthPillar )
    ? 1.0
    : ( zPos( 7 )-coord( ip3 + 3 ) ) / denomFourthPillar;
  xPos( 7 ) = slopePos7 * (coord( ip3 + 4 ) - coord( ip3 + 1 )) + coord( ip3 + 1 );
  yPos( 7 ) = slopePos7 * (coord( ip3 + 5 ) - coord( ip3 + 2 )) + coord( ip3 + 2 );
  cPVertexIsInside( 7 ) = fourthPillarIsInside;

  return true;
}


void CPMeshBuilder::filterVertices()
{
  array2d< real64 > const & cPVertices = m_meshVertices.m_cPVertices;
  array1d< bool > const & cPVertexIsInsidePartition = m_meshVertices.m_cPVertexIsInsidePartition;
  array1d< globalIndex > const & cPVertexToGlobalCPVertex = m_meshVertices.m_cPVertexToGlobalCPVertex;
  array1d< localIndex > & cPVertexToVertex = m_meshVertices.m_cPVertexToVertex;
  cPVertexToVertex.resize( cPVertices.size() );

  // First step: filter the unique vertices using a set

  std::set< Vertex, CompareVertices > uniqueVerticesHelper;
  std::set< Vertex, CompareVertices >::iterator it;

  // loop over of the CP vertices (for now, including those of inactive cells)
  for( localIndex iCPVertex = 0; iCPVertex < cPVertices.size( 0 ); ++iCPVertex )
  {
    // make sure cPVertices in the **exterior** pillars of the overlap are skipped
    if( cPVertexIsInsidePartition( iCPVertex ) )
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
        cPVertexToVertex( iCPVertex ) = existingVertex.m_localIndex;
        if( existingVertex.m_globalIndex > cPVertexToGlobalCPVertex( iCPVertex ) )
        {
          // make sure that the global index does not depend on the order in which cPVertices are processed
          // unfortunately I cannot directly modify existingVertex.m_globalIndex (it is const)
          // TODO: find out is there is a better way to do that
          v.m_localIndex = existingVertex.m_localIndex;
          v.m_globalIndex = cPVertexToGlobalCPVertex( iCPVertex );
          uniqueVerticesHelper.erase( existingVertex );
          uniqueVerticesHelper.insert( v );
        }
      }
      // if not found yet, insert into the set
      else
      {
        v.m_localIndex = uniqueVerticesHelper.size();
        v.m_globalIndex = cPVertexToGlobalCPVertex( iCPVertex );
        cPVertexToVertex( iCPVertex ) = v.m_localIndex;
        uniqueVerticesHelper.insert( v );
      }
    }
  }

  // Second step: move the data from the set to a array2d, because we don't want the set anymore
  // TODO: check if we really need this second step
  array2d< real64 > & vertices = m_meshVertices.m_vertices;
  array1d< globalIndex > & vertexToGlobalVertex = m_meshVertices.m_vertexToGlobalVertex;
  vertices.resizeDimension< 0, 1 >( uniqueVerticesHelper.size(), 3 );
  vertexToGlobalVertex.resize( uniqueVerticesHelper.size() );
  std::for_each( uniqueVerticesHelper.begin(), uniqueVerticesHelper.end(),
                 [&vertices, &vertexToGlobalVertex]( Vertex const & v )
  {
    vertices( v.m_localIndex, 0 ) = v.m_x;
    vertices( v.m_localIndex, 1 ) = v.m_y;
    vertices( v.m_localIndex, 2 ) = v.m_z;
    vertexToGlobalVertex( v.m_localIndex ) = v.m_globalIndex;
  } );
}


void CPMeshBuilder::buildFaces()
{
  // TODO
}

REGISTER_CATALOG_ENTRY( CPMeshBuilder, CPMeshBuilder, string const & )

} // namespace CPMesh

} // end namespace geosx
