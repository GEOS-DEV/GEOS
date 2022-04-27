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
 * @file CornerPointMeshBuilder.cpp
 */

#include "CornerPointMeshBuilder.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/MpiWrapper.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/generators/cornerPointMesh/utilities/GeometryUtilities.hpp"
#include "mesh/generators/cornerPointMesh/utilities/OutputUtilities.hpp"
#include <chrono>

namespace geosx
{

namespace cornerPointMesh
{

using namespace geometryUtilities;
using namespace outputUtilities;

CornerPointMeshBuilder::CornerPointMeshBuilder( string const & name )
  :
  m_parser( name ),
  m_partition( name ),
  m_meshName( name ),
  m_offset( 1.0 )
{}

void CornerPointMeshBuilder::buildMesh( Path const & filePath )
{
  localIndex const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  if( rank == 0 )
  {
    localIndex nX = 0;
    localIndex nY = 0;
    localIndex nZ = 0;

    // read SPECGRID/DIMENS to get the number of cells in each direction
    m_parser.readNumberOfCells( filePath, nX, nY, nZ );
    m_dims.defineDomainDimensions( nX, nY, nZ );

    // TODO: read ACTNUM for the inertial partitioning
  }

  // rank 0 partitions the mesh and sends partitions to other ranks
  // rank 0 sets up the communication pattern
  m_partition.setupMPIPartition( m_dims );

  // each rank reads its part of the mesh
  m_parser.readMesh( filePath, m_dims );

  // post-process the mesh to make it conforming
  postProcessMesh();
}

void CornerPointMeshBuilder::postProcessMesh()
{
  // for each active cell, compute the position of the eight "corner-point" vertices
  // at this point, "corner-point" vertices have not been filtered yet
  buildCornerPointCells();

  // eliminate duplicates from cpVertices
  filterVertices();

  // match faces, deal with the non-conforming case
  buildFaces();

  // construct map linking a region to its cells
  formRegions();

  // copy over the unique vertices (some new vertices may have been added while building faces)
  updateVerticesAfterFaceBuilding();

  // output VTK for testing purposes
  int const numRanks = MpiWrapper::commSize();
  if( numRanks == 1 )
  {
    outputDebugVTKFileWithFaces( m_vertices, m_faces, m_cells );
  }
}

void CornerPointMeshBuilder::updateVerticesAfterFaceBuilding()
{
  std::set< Vertex > const & uniqueVerticesHelper = m_vertices.m_uniqueVerticesHelper;
  array2d< real64 > & vertexPositions = m_vertices.m_vertexPositions;

  vertexPositions.resizeDimension< 0, 1 >( uniqueVerticesHelper.size(), 3 );
  std::for_each( uniqueVerticesHelper.begin(), uniqueVerticesHelper.end(),
                 [&vertexPositions]( geometryUtilities::Vertex const & v )
  {
    vertexPositions( v.m_localIndex, 0 ) = v.m_x;
    vertexPositions( v.m_localIndex, 1 ) = v.m_y;
    vertexPositions( v.m_localIndex, 2 ) = v.m_z;
  } );

}


namespace internal
{

bool computeCellCoordinates( localIndex const i, localIndex const j, localIndex const k,
                             localIndex const nXLocal, localIndex const nYLocal,
                             localIndex const iMinOverlap, localIndex const iMaxOverlap,
                             localIndex const jMinOverlap, localIndex const jMaxOverlap,
                             array1d< real64 > const & coord,
                             array1d< real64 > const & zcorn,
                             array1d< real64 > & xPos,
                             array1d< real64 > & yPos,
                             array1d< real64 > & zPos,
                             array1d< bool > & cpVertexIsInside )
{
  localIndex const iXmLocal = k*8*nXLocal*nYLocal + j*4*nXLocal + 2*i;
  localIndex const iXpLocal = iXmLocal + 2*nXLocal;

  zPos( 0 ) = zcorn( iXmLocal );
  zPos( 1 ) = zcorn( iXmLocal + 1 );
  zPos( 2 ) = zcorn( iXpLocal );
  zPos( 3 ) = zcorn( iXpLocal + 1 );
  zPos( 4 ) = zcorn( iXmLocal + 4*nXLocal*nYLocal );
  zPos( 5 ) = zcorn( iXmLocal + 4*nXLocal*nYLocal + 1 );
  zPos( 6 ) = zcorn( iXpLocal + 4*nXLocal*nYLocal );
  zPos( 7 ) = zcorn( iXpLocal + 4*nXLocal*nYLocal + 1 );

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
  localIndex const iFirstPillarLocal = j*(nXLocal+1)+i;
  localIndex const ip0 = 6*iFirstPillarLocal - 1;
  real64 const denomFirstPillar = coord( ip0 + 6 ) - coord( ip0 + 3 );
  bool const firstPillarIsInside = !( (i == 0 && iMinOverlap == 1) || (j == 0 && jMinOverlap == 1) );

  real64 const slopePos0 = isZero( denomFirstPillar )
    ? 1.0
    : ( zPos( 0 )-coord( ip0 + 3 ) ) / denomFirstPillar;
  xPos( 0 ) = slopePos0 * (coord( ip0 + 4 ) - coord( ip0 + 1 )) + coord( ip0 + 1 );
  yPos( 0 ) = slopePos0 * (coord( ip0 + 5 ) - coord( ip0 + 2 )) + coord( ip0 + 2 );
  cpVertexIsInside( 0 ) = firstPillarIsInside;

  real64 const slopePos4 = isZero( denomFirstPillar )
    ? 1.0
    : ( zPos( 4 )-coord( ip0 + 3 ) ) / denomFirstPillar;
  xPos( 4 ) = slopePos4 * (coord( ip0 + 4 ) - coord( ip0 + 1 )) + coord( ip0 + 1 );
  yPos( 4 ) = slopePos4 * (coord( ip0 + 5 ) - coord( ip0 + 2 )) + coord( ip0 + 2 );
  cpVertexIsInside( 4 ) = firstPillarIsInside;

  // second pillar
  localIndex const iSecondPillarLocal = (nXLocal+1)*j + i+1;
  localIndex const ip1 = 6*iSecondPillarLocal - 1;
  real64 const denomSecondPillar = coord( ip1 + 6 ) - coord( ip1 + 3 );
  bool const secondPillarIsInside = !( (i+1 == nXLocal && iMaxOverlap == 1) || (j == 0 && jMinOverlap == 1) );

  real64 const slopePos1 = isZero( denomSecondPillar )
    ? 1.0
    : ( zPos( 1 )-coord( ip1 + 3 ) ) / denomSecondPillar;
  xPos( 1 ) = slopePos1 * (coord( ip1 + 4 ) - coord( ip1 + 1 )) + coord( ip1 + 1 );
  yPos( 1 ) = slopePos1 * (coord( ip1 + 5 ) - coord( ip1 + 2 )) + coord( ip1 + 2 );
  cpVertexIsInside( 1 ) = secondPillarIsInside;

  real64 const slopePos5 = isZero( denomSecondPillar )
    ? 1.0
    : ( zPos( 5 )-coord( ip1 + 3 ) ) / denomSecondPillar;
  xPos( 5 ) = slopePos5 * (coord( ip1 + 4 ) - coord( ip1 + 1 )) + coord( ip1 + 1 );
  yPos( 5 ) = slopePos5 * (coord( ip1 + 5 ) - coord( ip1 + 2 )) + coord( ip1 + 2 );
  cpVertexIsInside( 5 ) = secondPillarIsInside;

  // third pillar
  localIndex const iThirdPillarLocal = (nXLocal+1)*(j+1) + i;
  localIndex const ip2 = 6*iThirdPillarLocal - 1;
  real64 const denomThirdPillar = coord( ip2 + 6 ) - coord( ip2 + 3 );
  bool const thirdPillarIsInside = !( (i == 0 && iMinOverlap == 1) || (j+1 == nYLocal && jMaxOverlap == 1) );

  real64 const slopePos2 = isZero( denomThirdPillar )
    ? 1.0
    : ( zPos( 2 )-coord( ip2 + 3 ) ) / denomThirdPillar;
  xPos( 2 ) = slopePos2 * (coord( ip2 + 4 ) - coord( ip2 + 1 )) + coord( ip2 + 1 );
  yPos( 2 ) = slopePos2 * (coord( ip2 + 5 ) - coord( ip2 + 2 )) + coord( ip2 + 2 );
  cpVertexIsInside( 2 ) = thirdPillarIsInside;

  real64 const slopePos6 = isZero( denomThirdPillar )
    ? 1.0
    : ( zPos( 6 )-coord( ip2 + 3 ) ) / denomThirdPillar;
  xPos( 6 ) = slopePos6 * (coord( ip2 + 4 ) - coord( ip2 + 1 )) + coord( ip2 + 1 );
  yPos( 6 ) = slopePos6 * (coord( ip2 + 5 ) - coord( ip2 + 2 )) + coord( ip2 + 2 );
  cpVertexIsInside( 6 ) = thirdPillarIsInside;

  // fourth pillar
  localIndex const iFourthPillarLocal = (nXLocal+1)*(j+1) + i+1;
  localIndex const ip3 = 6*iFourthPillarLocal - 1;
  real64 const denomFourthPillar = coord( ip3 + 6 ) - coord( ip3 + 3 );
  bool const fourthPillarIsInside = !( (i+1 == nXLocal && iMaxOverlap == 1) || (j+1 == nYLocal && jMaxOverlap == 1) );

  real64 const slopePos3 = isZero( denomFourthPillar )
    ? 1.0
    : ( zPos( 3 )-coord( ip3 + 3 ) ) / denomFourthPillar;
  xPos( 3 ) = slopePos3 * (coord( ip3 + 4 ) - coord( ip3 + 1 )) + coord( ip3 + 1 );
  yPos( 3 ) = slopePos3 * (coord( ip3 + 5 ) - coord( ip3 + 2 )) + coord( ip3 + 2 );
  cpVertexIsInside( 3 ) = fourthPillarIsInside;

  real64 const slopePos7 = isZero( denomFourthPillar )
    ? 1.0
    : ( zPos( 7 )-coord( ip3 + 3 ) ) / denomFourthPillar;
  xPos( 7 ) = slopePos7 * (coord( ip3 + 4 ) - coord( ip3 + 1 )) + coord( ip3 + 1 );
  yPos( 7 ) = slopePos7 * (coord( ip3 + 5 ) - coord( ip3 + 2 )) + coord( ip3 + 2 );
  cpVertexIsInside( 7 ) = fourthPillarIsInside;

  return true;
}

void populateReverseCellMaps( localIndex const nCells,
                              localIndex const nActiveCells,
                              array1d< localIndex > const & activeCellToCell,
                              array1d< localIndex > const & ownedActiveCellToActiveCell,
                              array1d< localIndex > & cellToActiveCell,
                              array1d< localIndex > & activeCellToOwnedActiveCell )
{
  cellToActiveCell.resize( nCells );
  activeCellToOwnedActiveCell.resize( nActiveCells );
  cellToActiveCell.setValues< serialPolicy >( -1 );
  activeCellToOwnedActiveCell.setValues< serialPolicy >( -1 );

  for( localIndex iActiveCell = 0; iActiveCell < activeCellToCell.size(); ++iActiveCell )
  {
    cellToActiveCell( activeCellToCell( iActiveCell ) ) = iActiveCell;
  }
  for( localIndex iOwnedActiveCell = 0; iOwnedActiveCell < ownedActiveCellToActiveCell.size(); ++iOwnedActiveCell )
  {
    activeCellToOwnedActiveCell( ownedActiveCellToActiveCell( iOwnedActiveCell ) ) = iOwnedActiveCell;
  }
}

bool computeFaceGeometry( std::vector< Face > const & adjacentFaces,
                          std::vector< Vertex > & newVertices )
{
  //  Performs 3D geometrical calculation for two intersected faces, both of which are represented by multiple points

  bool isValid;
  newVertices = adjacentFaces[0].findIntersectionPoints( adjacentFaces[1], isValid );
  return isValid;
}

bool checkFaceOverlap( array2d< real64 > const & vertexPositions,
                       localIndex const (&nextfaceVertices)[ 4 ],
                       localIndex const (&prevfaceVertices)[ 4 ],
                       std::vector< Face > & adjacentFaces )
{
  // TODO: this is super ugly. at some corner cases, the algorithm might fail
  // edges of next cell
  Vertex const pointNext0 = { vertexPositions( nextfaceVertices[0], 0 ),
                              vertexPositions( nextfaceVertices[0], 1 ),
                              vertexPositions( nextfaceVertices[0], 2 ) };

  Vertex const pointNext2 = { vertexPositions( nextfaceVertices[1], 0 ),
                              vertexPositions( nextfaceVertices[1], 1 ),
                              vertexPositions( nextfaceVertices[1], 2 ) };
  Line const nextEdge02( pointNext0, pointNext2 );

  Vertex const pointNext4 = { vertexPositions( nextfaceVertices[2], 0 ),
                              vertexPositions( nextfaceVertices[2], 1 ),
                              vertexPositions( nextfaceVertices[2], 2 ) };

  Vertex const pointNext6 = { vertexPositions( nextfaceVertices[3], 0 ),
                              vertexPositions( nextfaceVertices[3], 1 ),
                              vertexPositions( nextfaceVertices[3], 2 ) };

  Line const nextEdge46( pointNext4, pointNext6 );

  // edges of previous cell
  Vertex const pointPrev1 = { vertexPositions( prevfaceVertices[0], 0 ),
                              vertexPositions( prevfaceVertices[0], 1 ),
                              vertexPositions( prevfaceVertices[0], 2 ) };

  Vertex const pointPrev3 = { vertexPositions( prevfaceVertices[1], 0 ),
                              vertexPositions( prevfaceVertices[1], 1 ),
                              vertexPositions( prevfaceVertices[1], 2 ) };

  Line const prevEdge13( pointPrev1, pointPrev3 );

  Vertex const pointPrev5 = { vertexPositions( prevfaceVertices[2], 0 ),
                              vertexPositions( prevfaceVertices[2], 1 ),
                              vertexPositions( prevfaceVertices[2], 2 ) };

  Vertex const pointPrev7 = { vertexPositions( prevfaceVertices[3], 0 ),
                              vertexPositions( prevfaceVertices[3], 1 ),
                              vertexPositions( prevfaceVertices[3], 2 ) };

  Line const prevEdge57( pointPrev5, pointPrev7 );

  // TODO: is there a better way to compare the position between two cells

  bool const isPrevCellHigher = prevEdge13.compare( nextEdge46 );
  bool const isPrevCellLower = nextEdge02.compare( prevEdge57 );

  if( isPrevCellHigher || isPrevCellLower )
  {
    return false;
  }
  else
  {
    Face faceNext( pointNext0, pointNext2, pointNext4, pointNext6 );
    Face facePrev( pointPrev1, pointPrev3, pointPrev5, pointPrev7 );

    faceNext.setIndexForFacePoints( nextfaceVertices );
    facePrev.setIndexForFacePoints( prevfaceVertices );

    adjacentFaces.emplace_back( faceNext );
    adjacentFaces.emplace_back( facePrev );
  }
  return true;
}

void updateVerticalFaceMaps( localIndex const iOwnedActiveCellPrev,
                             localIndex const iOwnedActiveCellNext,
                             ArrayOfArrays< localIndex > & ownedActiveCellToFaces,
                             ArrayOfArrays< localIndex > & faceToVertices,
                             std::vector< Vertex > & newVertices,
                             std::set< Vertex > & uniqueVerticesHelper )
{
  //  some new points might be created when two faces intersect with each other. split these points, which latter
  //  be added to coord array and update corner index

  localIndex const newFaceId = faceToVertices.size();
  faceToVertices.appendArray( newVertices.size() );
  localIndex count( 0 );
  for( auto iter = newVertices.begin(); iter != newVertices.end(); ++iter )
  {
    if( (*iter).m_localIndex < 0 )
    {
      // find  a possible newly created vertex
      if( uniqueVerticesHelper.find( (*iter)) != uniqueVerticesHelper.end() )
      {
        // already existed and get the index
        localIndex const oldLocalIdx = (*uniqueVerticesHelper.find( (*iter))).m_localIndex;
        (*iter).setIndex( oldLocalIdx );
      }
      else
      {
        // a newly created point
        localIndex const newLocalIdx = uniqueVerticesHelper.size();
        (*iter).setIndex( newLocalIdx );
        uniqueVerticesHelper.insert( (*iter) );
      }
    }

    faceToVertices( newFaceId, count ) = (*iter).m_localIndex;
    ++count;
  }

  // populate the cell-to-face map
  if( iOwnedActiveCellPrev != -1 )
  {
    ownedActiveCellToFaces.emplaceBack( iOwnedActiveCellPrev, newFaceId );
  }
  if( iOwnedActiveCellNext != -1 )
  {
    ownedActiveCellToFaces.emplaceBack( iOwnedActiveCellNext, newFaceId );
  }

}


};


void CornerPointMeshBuilder::buildCornerPointCells()
{
  // local and global dimensions
  localIndex const nX = m_dims.nX();
  localIndex const nY = m_dims.nY();

  localIndex const iMinLocal = m_dims.iMinLocal();
  localIndex const jMinLocal = m_dims.jMinLocal();
  localIndex const kMinLocal = 0;

  localIndex const iMinOverlap = m_dims.iMinOverlap();
  localIndex const jMinOverlap = m_dims.jMinOverlap();
  localIndex const iMaxOverlap = m_dims.iMaxOverlap();
  localIndex const jMaxOverlap = m_dims.jMaxOverlap();

  localIndex const nXLocal = m_dims.nXLocal();
  localIndex const nYLocal = m_dims.nYLocal();
  localIndex const nZLocal = m_dims.nZLocal() + 2; // add two layers of auxillary cells for building faces

  // Auxillary layers need to be appended to coord, zcorn, and actnum before computing vertice coordinates
  appendAuxillaryLayer();

  localIndex const nLocalCells = nXLocal*nYLocal*nZLocal;
  localIndex constexpr nCPVerticesPerCell = 8;
  localIndex const nCPVertices = nCPVerticesPerCell*nLocalCells;

  array1d< localIndex > const & actnum = m_parser.actnum();
  array1d< real64 > const & coord = m_parser.coord();
  array1d< real64 > const & zcorn = m_parser.zcorn();

  // vertex maps
  array2d< real64 > & cpVertexPositions = m_vertices.m_cpVertexPositions;
  array1d< globalIndex > & cpVertexToGlobalCPVertex = m_vertices.m_cpVertexToGlobalCPVertex;
  array1d< bool > & cpVertexIsInsidePartition = m_vertices.m_cpVertexIsInsidePartition;
  cpVertexPositions.resizeDimension< 0, 1 >( nCPVertices, 3 );
  cpVertexToGlobalCPVertex.resize( nCPVertices );
  cpVertexIsInsidePartition.resize( nCPVertices );

  // cell maps
  array1d< globalIndex > & ownedActiveCellToGlobalCell = m_cells.m_ownedActiveCellToGlobalCell;
  array1d< localIndex > & ownedActiveCellToActiveCell = m_cells.m_ownedActiveCellToActiveCell;
  array1d< localIndex > & activeCellToCell = m_cells.m_activeCellToCell;
  array1d< localIndex > & cellToCPVertices = m_cells.m_cellToCPVertices;
  ownedActiveCellToGlobalCell.reserve( nLocalCells );
  ownedActiveCellToActiveCell.reserve( nLocalCells );
  activeCellToCell.reserve( nLocalCells );
  cellToCPVertices.resize( nLocalCells );

  array1d< real64 > xPos( nCPVerticesPerCell );
  array1d< real64 > yPos( nCPVerticesPerCell );
  array1d< real64 > zPos( nCPVerticesPerCell );

  array1d< bool > cpVertexIsInside( nCPVerticesPerCell );

  // loop over all cells in the MPI domain, including overlaps
  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        bool isValid = false;
        // compute explicit local and global indices using the ijk structure
        localIndex const iLocalCell = k*nXLocal*nYLocal + j*nXLocal + i;
        globalIndex const iGlobalCell = (k+kMinLocal)*nX*nY + (j+jMinLocal)*nX + (i+iMinLocal);

        //compute the positions of the eight vertices
        isValid =
          internal::computeCellCoordinates( i, j, k,
                                            nXLocal, nYLocal,
                                            iMinOverlap, iMaxOverlap,
                                            jMinOverlap, jMaxOverlap,
                                            coord, zcorn,
                                            xPos, yPos, zPos, cpVertexIsInside );

        if( isValid )
        {
          // cell maps
          // assign local and global indices
          if( actnum( iLocalCell ) == 1 )
          {
            activeCellToCell.emplace_back( iLocalCell );
            // decide flag specifying whether the cell is in the overlap or not
            if( m_dims.columnIsInsidePartition( i, j ) )
            {
              ownedActiveCellToActiveCell.emplace_back( activeCellToCell.size()-1 );
              ownedActiveCellToGlobalCell.emplace_back( iGlobalCell );
            }
          }

          // vertex maps

          // construct the map from local cell to first CP vertex of the cell
          localIndex const iFirstVertexLocal = nCPVerticesPerCell * iLocalCell;
          globalIndex const iFirstVertexGlobal = nCPVerticesPerCell * iGlobalCell;
          cellToCPVertices( iLocalCell ) = iFirstVertexLocal;

          //      6__ __ __7
          //     /        /|
          //    /        / |
          //   /        /  |    z
          //  4__ __ __5   3    |__ x
          //  |        |  /    /
          //  |        | /    y
          //  |        |/
          //  0__ __ __1

          // save the position of the eight vertices, and assign global CP vertex indices
          localIndex const order[nCPVerticesPerCell] = { 4, 5, 6, 7, 0, 1, 2, 3 };
          for( localIndex pos = 0; pos < nCPVerticesPerCell; ++pos )
          {
            cpVertexPositions( iFirstVertexLocal + pos, 0 ) = xPos( order[pos] );
            cpVertexPositions( iFirstVertexLocal + pos, 1 ) = yPos( order[pos] );
            cpVertexPositions( iFirstVertexLocal + pos, 2 ) = -zPos( order[pos] );
            cpVertexToGlobalCPVertex( iFirstVertexLocal + pos ) = iFirstVertexGlobal + order[pos];
            cpVertexIsInsidePartition( iFirstVertexLocal + pos ) = cpVertexIsInside( order[pos] );
          }
        }
      }
    }
  }
  internal::populateReverseCellMaps( nLocalCells, activeCellToCell.size(),
                                     activeCellToCell, ownedActiveCellToActiveCell,
                                     m_cells.m_cellToActiveCell,
                                     m_cells.m_activeCellToOwnedActiveCell );
}


void CornerPointMeshBuilder::appendAuxillaryLayer()
{
  // Append top and bottom auxillary layers to zcorn and adjust actnum
  array1d< real64 > & zcorn = m_parser.zcorn();
  localIndex const nXLocal = m_dims.nXLocal() * 2;
  localIndex const nYLocal = m_dims.nYLocal() * 2;

  // step 1: Find highest and lowest points to set as up and bottom boundaries
  real64 const minZcorn = *( std::min_element( zcorn.begin(), zcorn.end() ) ) - m_offset;
  real64 const maxZcorn = *( std::max_element( zcorn.begin(), zcorn.end() ) ) + m_offset;

  // loop over all zcorns of elements at top/bottom layer in the MPI domain, including overlaps
  localIndex const k = m_dims.nZLocal() * 2; // lower layer

  for( localIndex j = 0; j < nYLocal; ++j )
  {
    for( localIndex i = 0; i < nXLocal; ++i )
    {
      // compute explicit local and global indices using the ijk structure
      // Step 2: assign values to zcorn array associated with auxillary layers
      localIndex const iLocalTopZcorn = j*nXLocal + i + 2* nYLocal * nXLocal;
      localIndex const iLocalBottomZcorn = iLocalTopZcorn + (k -1) * nYLocal * nXLocal;

      zcorn( iLocalTopZcorn -  nYLocal * nXLocal ) = zcorn( iLocalTopZcorn );
      zcorn( iLocalBottomZcorn +  nYLocal * nXLocal )= zcorn( iLocalBottomZcorn );

      zcorn( iLocalTopZcorn - 2*nYLocal * nXLocal ) = minZcorn;
      zcorn( iLocalBottomZcorn + 2*nYLocal * nXLocal ) = maxZcorn;
    }
  }

  appendCellDataForAuxillaryLayer();
}

void CornerPointMeshBuilder::appendCellDataForAuxillaryLayer()
{
  // regionID and actnum
  array1d< localIndex > & regionId = m_parser.regionId();
  array1d< localIndex > & actnum = m_parser.actnum();
  // Field data
  array1d< real64 > & poro = m_parser.poro();
  array2d< real64 > & perm = m_parser.perm();
  bool const hasPoro = ( poro.size() > 0 );
  bool const hasPerm = ( perm.size( 0 ) > 0 );

  localIndex const nXLocal = m_dims.nXLocal();
  localIndex const nYLocal = m_dims.nYLocal();
  localIndex const k = m_dims.nZLocal();

  // default values for auxillary cells

  for( localIndex j = 0; j < nYLocal; ++j )
  {
    for( localIndex i = 0; i < nXLocal; ++i )
    {
      localIndex const iLocalTopCell = j*nXLocal + i + nYLocal * nXLocal;
      localIndex const iLocalBottomCell = iLocalTopCell + (k-1) * nYLocal * nXLocal;

      localIndex const iLocalTopAuxCell = j*nXLocal + i;
      localIndex const iLocalBottomAuxCell = iLocalBottomCell + nYLocal * nXLocal;

      regionId( iLocalTopAuxCell ) = regionId( iLocalTopCell );
      actnum( iLocalTopAuxCell ) = -1;

      regionId( iLocalBottomAuxCell ) = regionId( iLocalBottomCell );
      actnum( iLocalBottomAuxCell ) = -1;

      if( hasPoro )
      {
        poro( iLocalTopAuxCell ) = 0.0;
        poro( iLocalBottomAuxCell ) = 0.0;
      }

      if( hasPerm )
      {
        for( localIndex dim = 0; dim < 3; ++dim )
        {
          perm( iLocalTopAuxCell, dim ) = 0.0;
          perm( iLocalBottomAuxCell, dim ) = 0.0;
        }
      }
    }
  }
}


void CornerPointMeshBuilder::formRegions()
{
  array1d< localIndex > const & activeCellToCell = m_cells.m_activeCellToCell;
  array1d< localIndex > const & ownedActiveCellToActiveCell = m_cells.m_ownedActiveCellToActiveCell;
  array1d< localIndex > const & regionId = m_parser.regionId();

  // Step 1: find the largest region index
  localIndex const nOwnedActiveCells = ownedActiveCellToActiveCell.size();
  localIndex localMaxRegionId = 0;
  for( localIndex iOwnedActiveCell = 0; iOwnedActiveCell < nOwnedActiveCells; ++iOwnedActiveCell )
  {
    localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
    localIndex const iCell = activeCellToCell( iActiveCell );
    localIndex const rid = regionId( iCell );
    if( rid > localMaxRegionId )
    {
      localMaxRegionId = rid;
    }
  }
  localIndex const globalMaxRegionId = MpiWrapper::max( localMaxRegionId );

  // Step 2: construct the map region to elements
  ArrayOfArrays< localIndex > & regionToOwnedActiveCell = m_regions.m_regionToOwnedActiveCell;
  localIndex const regionCapacity =
    static_cast< real64 >( nOwnedActiveCells )
    / static_cast< real64 >( globalMaxRegionId+1 );
  regionToOwnedActiveCell.resize( globalMaxRegionId+1, regionCapacity );

  for( localIndex iOwnedActiveCell = 0; iOwnedActiveCell < nOwnedActiveCells; ++iOwnedActiveCell )
  {
    localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
    localIndex const iCell = activeCellToCell( iActiveCell );

    regionToOwnedActiveCell.emplaceBack( regionId( iCell ), iOwnedActiveCell );
  }
}


void CornerPointMeshBuilder::filterVertices()
{
  array2d< real64 > const & cpVertexPositions = m_vertices.m_cpVertexPositions;
  array1d< bool > const & cpVertexIsInsidePartition = m_vertices.m_cpVertexIsInsidePartition;
  array1d< globalIndex > const & cpVertexToGlobalCPVertex = m_vertices.m_cpVertexToGlobalCPVertex;
  array1d< localIndex > & cpVertexToVertex = m_vertices.m_cpVertexToVertex;
  cpVertexToVertex.resize( cpVertexPositions.size() );

  // First step: filter the unique vertices using a set
  // the compareVertex struct has been replaced less than operator for set's sortting purpose
  std::set< Vertex > & uniqueVerticesHelper = m_vertices.m_uniqueVerticesHelper;
  std::set< Vertex >::iterator it;

  // loop over of the CP vertices (for now, including those of inactive cells)
  for( localIndex iCPVertex = 0; iCPVertex < cpVertexPositions.size( 0 ); ++iCPVertex )
  {
    // make sure cpVertices in the **exterior** pillars of the overlap are skipped
    if( cpVertexIsInsidePartition( iCPVertex ) )
    {
      Vertex v( cpVertexPositions( iCPVertex, 0 ),
                cpVertexPositions( iCPVertex, 1 ),
                cpVertexPositions( iCPVertex, 2 ) );

      // check if this vertex has already been found and inserted in the set
      it = uniqueVerticesHelper.find( v );

      // if already found, copy the local index in map and move on
      if( it != uniqueVerticesHelper.end() )
      {
        Vertex const & existingVertex = *it;
        cpVertexToVertex( iCPVertex ) = existingVertex.m_localIndex;
        if( existingVertex.m_globalIndex > cpVertexToGlobalCPVertex( iCPVertex ) )
        {
          // make sure that the global index does not depend on the order in which cpVertices are processed
          // unfortunately I cannot directly modify existingVertex.m_globalIndex (it is const)
          // TODO: find out is there is a better way to do that
          v.m_localIndex = existingVertex.m_localIndex;
          v.m_globalIndex = cpVertexToGlobalCPVertex( iCPVertex );
          uniqueVerticesHelper.erase( existingVertex );
          uniqueVerticesHelper.insert( v );
        }
      }
      // if not found yet, insert into the set
      else
      {
        v.m_localIndex = uniqueVerticesHelper.size();
        v.m_globalIndex = cpVertexToGlobalCPVertex( iCPVertex );
        cpVertexToVertex( iCPVertex ) = v.m_localIndex;
        uniqueVerticesHelper.insert( v );
      }
    }
  }

  // Second step: move the data from the set to a array2d, because we don't want the set anymore
  // TODO: check if we really need this second step
  array2d< real64 > & vertexPositions = m_vertices.m_vertexPositions;
  array1d< globalIndex > & vertexToGlobalVertex = m_vertices.m_vertexToGlobalVertex;
  vertexPositions.resizeDimension< 0, 1 >( uniqueVerticesHelper.size(), 3 );
  vertexToGlobalVertex.resize( uniqueVerticesHelper.size() );
  std::for_each( uniqueVerticesHelper.begin(), uniqueVerticesHelper.end(),
                 [&vertexPositions, &vertexToGlobalVertex]( Vertex const & v )
  {
    vertexPositions( v.m_localIndex, 0 ) = v.m_x;
    vertexPositions( v.m_localIndex, 1 ) = v.m_y;
    vertexPositions( v.m_localIndex, 2 ) = v.m_z;
    vertexToGlobalVertex( v.m_localIndex ) = v.m_globalIndex;
  } );
}


void CornerPointMeshBuilder::buildFaces()
{
  // not ready for review, in construction

  // allocate space for the maps
  localIndex const nOwnedActiveCells = m_cells.m_ownedActiveCellToActiveCell.size();
  // TODO: we might need to dynamically allocate memory when finding more intersection points, which in turn creates more faces
  localIndex const faceCapacity = 6; // This hard-coded 6 only works for conforming cartesian grid.

  ArrayOfArrays< localIndex > & ownedActiveCellToFaces = m_cells.m_ownedActiveCellToFaces;
  ownedActiveCellToFaces.resize( nOwnedActiveCells, faceCapacity );
  // TODO: reserve space for faceToNodes

  // TODO: the construction of faces involves loops in the (j,i,k) order,
  //       which is not ideal giving the layout of the arrays used here.
  //       See if it is worth permuting the arrays before this step.

  // Step 1: construct faces in the X and Y directions (between columns of cells)
  buildVerticalFaces( "X" );
  buildVerticalFaces( "Y" );

  // Step 2: construct faces in the Z direction (in a column of cells)
  buildHorizontalFaces();
}

void CornerPointMeshBuilder::buildVerticalFaces( string const & direction )
{

  //  Finds all vertical faces perpendicular to either i or j direction. Two steps are involved: first find all
  //  regular faces that contain internal, not faulted faces and boundary faces where only one cell is connected.
  //  Second, find all faulted faces and compute intersections between them.

  localIndex const nXLocal = m_dims.nXLocal();
  localIndex const nYLocal = m_dims.nYLocal();
  localIndex const nZLocal = m_dims.nZLocal();

  // cells
  array1d< localIndex > const & cellToCPVertices = m_cells.m_cellToCPVertices;
  array1d< localIndex > const & cellToActiveCell = m_cells.m_cellToActiveCell;
  array1d< localIndex > const & activeCellToOwnedActiveCell = m_cells.m_activeCellToOwnedActiveCell;
  ArrayOfArrays< localIndex > & ownedActiveCellToFaces = m_cells.m_ownedActiveCellToFaces;

  // faces
  ArrayOfArrays< localIndex > & faceToVertices = m_faces.m_faceToVertices;

  // vertices
  array1d< localIndex > const & cpVertexToVertex = m_vertices.m_cpVertexToVertex;

  // index switch. i, j corresponding to x, y
  localIndex iOffset, jOffset, iPillar, jPillar;
  localIndex iPadding = 0;
  localIndex jPadding = 0;
  bool isInternalFace = true;
  bool isXAxis = true;
  bool onFrontBoundary = false;
  bool onBackBoundary = false;

  if( direction == "X" )
  {
    iOffset = 1;
    jOffset = 0;
    iPillar = nXLocal + 1;
    jPillar = nYLocal;
  }
  else if( direction == "Y" )
  {
    iOffset = 0;
    jOffset = 1;
    iPillar = nXLocal;
    jPillar = nYLocal + 1;
    isXAxis = false;
  }
  else
  {
    GEOSX_THROW( "Wrong direction. Please check input!", InputError );
  }

  // loop over pillars in the current direction, given by the input direction.
  // loop over cells in other directions.
  // Cells at z direction has 2 auxillary cells.

  // loop over pillar or cell
  for( localIndex i = 0; i < iPillar; ++i )
  {
    // loop over pillar or cell
    for( localIndex j = 0; j < jPillar; ++j )
    {

      if( isXAxis )
      {
        // Found the boundary faces perpendicular to x axis
        onFrontBoundary = (i == 0);
        onBackBoundary = (i == nXLocal);
        isInternalFace = !onFrontBoundary && !onBackBoundary;
        iPadding = (i == nXLocal) ? 1 : 0;
      }
      else
      {
        // Found the boundary faces perpendicular to y axis
        onFrontBoundary = (j == 0);
        onBackBoundary = (j == nYLocal);
        isInternalFace = !onFrontBoundary && !onBackBoundary;
        jPadding = (j == nYLocal) ? 1 : 0;
      }

      // loop over owned cells only
      if( m_dims.columnIsInsidePartition( i - iPadding, j - jPadding ) )
      {
        // local index of the previous cell
        localIndex iCellPrev = -1;
        // local index of the next cell
        localIndex iCellNext = -1;
        // local index of the previous active cell
        localIndex iActiveCellPrev = -1;
        // local index of the next active cell
        localIndex iActiveCellNext = -1;
        // loop over cells in k direction
        for( localIndex k = 0; k < nZLocal + 2; ++k )
        {
          if( !onFrontBoundary )
          {
            iCellPrev = k*nXLocal*nYLocal + (j - jOffset)*nXLocal + i - iOffset;
            iActiveCellPrev = cellToActiveCell( iCellPrev );
          }
          if( !onBackBoundary )
          {
            iCellNext = k*nXLocal*nYLocal + j*nXLocal + i;
            iActiveCellNext = cellToActiveCell( iCellNext );
          }

          bool const prevIsActive = (iActiveCellPrev != -1);
          bool const nextIsActive = (iActiveCellNext != -1);

          // Pre-defined order for finding faces. Faces along different axis have varying orders
          // order arrays are initialized with orders along x axis

          std::vector< localIndex > orderPrev( 4 );
          std::vector< localIndex > orderNext( 4 );

          if( isXAxis )
          {
            orderPrev = { 1, 3, 7, 5 };
            orderNext = { 0, 2, 6, 4 };
          }
          else
          {
            orderPrev = { 2, 3, 7, 6 };
            orderNext = { 0, 1, 5, 4 };
          }

          if( isInternalFace )
          {
            // we need to consider the intersection between auxillary cells and internal cells
            localIndex const iFirstVertexPrev = cellToCPVertices( iCellPrev );
            localIndex const iFirstVertexNext = cellToCPVertices( iCellNext );

            // get internal faces
            // get the vertices of the previous faces
            localIndex const verticesPrevFace[4] = { cpVertexToVertex( iFirstVertexPrev + orderPrev[0] ),
                                                     cpVertexToVertex( iFirstVertexPrev + orderPrev[1] ),
                                                     cpVertexToVertex( iFirstVertexPrev + orderPrev[2] ),
                                                     cpVertexToVertex( iFirstVertexPrev + orderPrev[3] ) };
            // get the vertices of the next faces
            localIndex const verticesNextFace[4] = { cpVertexToVertex( iFirstVertexNext + orderNext[0] ),
                                                     cpVertexToVertex( iFirstVertexNext + orderNext[1] ),
                                                     cpVertexToVertex( iFirstVertexNext + orderNext[2] ),
                                                     cpVertexToVertex( iFirstVertexNext + orderNext[3] ) };

            if( std::equal( std::begin( verticesPrevFace ), std::end( verticesPrevFace ),
                            std::begin( verticesNextFace ) ))
            {
              // a conformed face is found
              if( prevIsActive && nextIsActive )
              {
                // both cells are active
                addConformingFace( verticesPrevFace,
                                   activeCellToOwnedActiveCell( iActiveCellPrev ),
                                   activeCellToOwnedActiveCell( iActiveCellNext ),
                                   ownedActiveCellToFaces,
                                   faceToVertices );
              }
              else if( prevIsActive || nextIsActive )
              {
                // either prev or next cell is inactive
                addSingleConformingFace( prevIsActive,
                                         orderPrev, orderNext,
                                         iActiveCellPrev, iActiveCellNext,
                                         iCellPrev, iCellNext );
              }
            }
            else
            {
              // a faulted face is found. Need to loop over a column of cells to find all intersection points
              addNonConformingFace( iCellPrev, k, verticesNextFace,
                                    orderPrev, iActiveCellNext );

            }
          }
          else if( prevIsActive || nextIsActive )
          {
            // this is either: an interior face between an active cell and an inactive cell
            //             or: an boundary face at the bottom
            addSingleConformingFace( prevIsActive,
                                     orderPrev, orderNext,
                                     iActiveCellPrev, iActiveCellNext,
                                     iCellPrev, iCellNext );

          }
        }
      }
    }
  }
}

void CornerPointMeshBuilder::addSingleConformingFace( bool const prevIsActive,
                                                      std::vector< localIndex > const & orderPrev,
                                                      std::vector< localIndex > const & orderNext,
                                                      localIndex const iActiveCellPrev,
                                                      localIndex const iActiveCellNext,
                                                      localIndex const iCellPrev,
                                                      localIndex const iCellNext )
{
  array1d< localIndex > const & cellToCPVertices = m_cells.m_cellToCPVertices;
  array1d< localIndex > const & cpVertexToVertex = m_vertices.m_cpVertexToVertex;
  ArrayOfArrays< localIndex > & ownedActiveCellToFaces = m_cells.m_ownedActiveCellToFaces;
  array1d< localIndex > const & activeCellToOwnedActiveCell = m_cells.m_activeCellToOwnedActiveCell;
  // faces
  ArrayOfArrays< localIndex > & faceToVertices = m_faces.m_faceToVertices;

  localIndex const iFirstVertex = prevIsActive ? cellToCPVertices( iCellPrev ) : cellToCPVertices( iCellNext );
  localIndex const iActiveCell = prevIsActive ? iActiveCellPrev : iActiveCellNext;
  std::vector< localIndex > vertexOrder( 4 );
  vertexOrder = prevIsActive ? orderPrev : orderNext;

  localIndex const vertices[4] = { cpVertexToVertex( iFirstVertex + vertexOrder[0] ),
                                   cpVertexToVertex( iFirstVertex + vertexOrder[1] ),
                                   cpVertexToVertex( iFirstVertex + vertexOrder[2] ),
                                   cpVertexToVertex( iFirstVertex + vertexOrder[3] )};

  addConformingFace( vertices,
                     activeCellToOwnedActiveCell( iActiveCell ),
                     -1,
                     ownedActiveCellToFaces,
                     faceToVertices );


}

void CornerPointMeshBuilder::addNonConformingFace( localIndex const iCellPrev,
                                                   localIndex const zIdxPrev,
                                                   localIndex const (&nextFaceVertices)[ 4 ],
                                                   std::vector< localIndex > const & orderPrev,
                                                   localIndex const iActiveCellNext )
{
  // loop over the column of cells that have the same i and j values with that of the previous cell
  // loop over cells in k direction

  // TODO: to see if there is a better way to add nonconforming faces
  localIndex const nXLocal = m_dims.nXLocal();
  localIndex const nYLocal = m_dims.nYLocal();
  localIndex const nZLocal = m_dims.nZLocal();

  // cells
  array1d< localIndex > const & cellToCPVertices = m_cells.m_cellToCPVertices;
  array1d< localIndex > const & cellToActiveCell = m_cells.m_cellToActiveCell;
  array1d< localIndex > const & activeCellToOwnedActiveCell = m_cells.m_activeCellToOwnedActiveCell;
  ArrayOfArrays< localIndex > & ownedActiveCellToFaces = m_cells.m_ownedActiveCellToFaces;

  // faces
  ArrayOfArrays< localIndex > & faceToVertices = m_faces.m_faceToVertices;

  // vertices
  array1d< localIndex > const & cpVertexToVertex = m_vertices.m_cpVertexToVertex;
  array2d< real64 > const & vertexPositions = m_vertices.m_vertexPositions;

  localIndex const firstCellAtPrevColumn = iCellPrev - zIdxPrev*nXLocal*nYLocal;
  localIndex newCellPrev = -1;

  // TODO: current strategy is to loop from the very top to the bottom, which might be quite inefficient since intersection often
  // occurs in a localized region

  for( localIndex k = 0; k < nZLocal + 2; ++k )
  {
    newCellPrev = k*nXLocal*nYLocal + firstCellAtPrevColumn;
    localIndex const iActiveCellPrev = cellToActiveCell( newCellPrev );
    localIndex const iFirstVertexPrev = cellToCPVertices( newCellPrev );
    // get internal faces
    // get the vertices of the previous faces
    localIndex const verticesPrevFace[4] = { cpVertexToVertex( iFirstVertexPrev + orderPrev[0] ),
                                             cpVertexToVertex( iFirstVertexPrev + orderPrev[1] ),
                                             cpVertexToVertex( iFirstVertexPrev + orderPrev[2] ),
                                             cpVertexToVertex( iFirstVertexPrev + orderPrev[3] ) };

    std::vector< Face > adjacentFaces;
    bool const isIntersecting = internal::checkFaceOverlap( vertexPositions,
                                                            nextFaceVertices,
                                                            verticesPrevFace,
                                                            adjacentFaces );

    if( isIntersecting && ((iActiveCellPrev != -1) || (iActiveCellNext != -1)) )
    {
      // TODO: we can check if the intersected faces satisfy certain conditions:
      //       1. area check
      //       2. cell volume check

      // TODO: we need a computational geometry function
      std::vector< Vertex > newVertices;
      bool const foundIntersection = internal::computeFaceGeometry( adjacentFaces,
                                                                    newVertices );

      localIndex const iOwnedActiveCellPrev = (iActiveCellPrev != -1) ? activeCellToOwnedActiveCell( iActiveCellPrev ) : -1;
      localIndex const iOwnedActiveCellNext = (iActiveCellNext != -1) ? activeCellToOwnedActiveCell( iActiveCellNext ) : -1;

      if( foundIntersection )
      {
        internal::updateVerticalFaceMaps( iOwnedActiveCellPrev,
                                          iOwnedActiveCellNext,
                                          ownedActiveCellToFaces,
                                          faceToVertices,
                                          newVertices,
                                          m_vertices.m_uniqueVerticesHelper );
      }
      else
      {
        GEOSX_THROW( "Cannot find intersections from two candidate nonconforming faces!", InputError );
      }
    }
    else
    {
      continue; // two faces are impossible to intersect with each other
    }
  }
}

void CornerPointMeshBuilder::buildHorizontalFaces()
{
  localIndex const nXLocal = m_dims.nXLocal();
  localIndex const nYLocal = m_dims.nYLocal();
  localIndex const nZLocal = m_dims.nZLocal() + 2;

  // cells
  array1d< localIndex > const & cellToCPVertices = m_cells.m_cellToCPVertices;
  array1d< localIndex > const & cellToActiveCell = m_cells.m_cellToActiveCell;
  array1d< localIndex > const & activeCellToOwnedActiveCell = m_cells.m_activeCellToOwnedActiveCell;
  ArrayOfArrays< localIndex > & ownedActiveCellToFaces = m_cells.m_ownedActiveCellToFaces;

  // faces
  ArrayOfArrays< localIndex > & faceToVertices = m_faces.m_faceToVertices;

  // vertices
  array1d< localIndex > const & cpVertexToVertex = m_vertices.m_cpVertexToVertex;

  // loop over owned cells only
  for( localIndex j = 0; j < nYLocal; ++j )
  {
    for( localIndex i = 0; i < nXLocal; ++i )
    {
      if( m_dims.columnIsInsidePartition( i, j ) )
      {

        // we are in column of cells, we can now loop over all the faces in this column
        localIndex const nFacesInColumn = nZLocal+1;
        for( localIndex iface = 0; iface < nFacesInColumn; ++iface )
        {
          // k-index of the cell before the face
          localIndex const kPrev = iface-1;
          // k-index of the cell after the face
          localIndex const kNext = iface;

          // local index of the previous cell
          localIndex iCellPrev = -1;
          // local index of the next cell
          localIndex iCellNext = -1;
          // local index of the previous active cell
          localIndex iActiveCellPrev = -1;
          // local index of the next active cell
          localIndex iActiveCellNext = -1;

          bool const onTopBoundary = (kPrev == -1);
          bool const onBottomBoundary = (kNext == nZLocal);
          bool const isInInterior = !onTopBoundary && !onBottomBoundary;
          if( !onTopBoundary )
          {
            iCellPrev = kPrev*nXLocal*nYLocal + j*nXLocal + i;
            iActiveCellPrev = cellToActiveCell( iCellPrev );
          }
          if( !onBottomBoundary )
          {
            iCellNext = kNext*nXLocal*nYLocal + j*nXLocal + i;
            iActiveCellNext = cellToActiveCell( iCellNext );
          }

          bool const prevIsActive = (iActiveCellPrev != -1);
          bool const nextIsActive = (iActiveCellNext != -1);
          // this is an interior cell between active cells
          if( isInInterior && prevIsActive && nextIsActive )
          {
            // get the bottom vertices of the previous face
            localIndex const iFirstVertexPrev = cellToCPVertices( iCellPrev );
            localIndex const vertices[4] = { cpVertexToVertex( iFirstVertexPrev ),
                                             cpVertexToVertex( iFirstVertexPrev + 1 ),
                                             cpVertexToVertex( iFirstVertexPrev + 3 ),
                                             cpVertexToVertex( iFirstVertexPrev + 2 ) };
            addConformingFace( vertices,
                               activeCellToOwnedActiveCell( iActiveCellPrev ),
                               activeCellToOwnedActiveCell( iActiveCellNext ),
                               ownedActiveCellToFaces,
                               faceToVertices );
          }
          // this is either: an interior face between an active cell and an inactive cell
          //             or: an boundary face at the bottom
          else if( prevIsActive || nextIsActive )
          {
            // get the bottom vertices of the previous face
            localIndex const iFirstVertex = prevIsActive ? cellToCPVertices( iCellPrev ) : cellToCPVertices( iCellNext );
            localIndex const iActiveCell = prevIsActive ? iActiveCellPrev : iActiveCellNext;
            localIndex const firstPos = prevIsActive ? 0 : 4;
            localIndex const vertices[4] = { cpVertexToVertex( iFirstVertex + firstPos ),
                                             cpVertexToVertex( iFirstVertex + firstPos + 1 ),
                                             cpVertexToVertex( iFirstVertex + firstPos + 3 ),
                                             cpVertexToVertex( iFirstVertex + firstPos + 2 ) };
            addConformingFace( vertices,
                               activeCellToOwnedActiveCell( iActiveCell ),
                               -1,
                               ownedActiveCellToFaces,
                               faceToVertices );
          }
        }
      }
    }
  }
}


void CornerPointMeshBuilder::addConformingFace( localIndex const (&faceVertices)[ 4 ],
                                                localIndex const iOwnedActiveCellPrev,
                                                localIndex const iOwnedActiveCellNext,
                                                ArrayOfArrays< localIndex > & ownedActiveCellToFaces,
                                                ArrayOfArrays< localIndex > & faceToVertices )
{
  // populate the face-to-vertices map
  localIndex const newFaceId = faceToVertices.size();
  faceToVertices.appendArray( 4 );
  faceToVertices( newFaceId, 0 ) = faceVertices[ 0 ];
  faceToVertices( newFaceId, 1 ) = faceVertices[ 1 ];
  faceToVertices( newFaceId, 2 ) = faceVertices[ 2 ];
  faceToVertices( newFaceId, 3 ) = faceVertices[ 3 ];

  // populate the cell-to-face map
  if( iOwnedActiveCellPrev != -1 )
  {
    ownedActiveCellToFaces.emplaceBack( iOwnedActiveCellPrev, newFaceId );
  }
  if( iOwnedActiveCellNext != -1 )
  {
    ownedActiveCellToFaces.emplaceBack( iOwnedActiveCellNext, newFaceId );
  }
}

REGISTER_CATALOG_ENTRY( CornerPointMeshBuilder, CornerPointMeshBuilder, string const & )

} // namespace cornerPointMesh

} // namespace geosx
