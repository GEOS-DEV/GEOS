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
 * @file NodeManager.cpp
 */

#include "NodeManager.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/BufferOps.hpp"
#include "mesh/EdgeManager.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/ToElementRelation.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"

namespace geos
{

using namespace dataRepository;

// *********************************************************************************************************************
/**
 * @return
 */
//START_SPHINX_REFPOS_REG
NodeManager::NodeManager( string const & name,
                          Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition( 0, 3 )
{
  registerWrapper( viewKeyStruct::referencePositionString(), &m_referencePosition );
  //END_SPHINX_REFPOS_REG
  this->registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );
  this->registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation );
  this->registerWrapper( viewKeyStruct::elementRegionListString(), &elementRegionList() );
  this->registerWrapper( viewKeyStruct::elementSubRegionListString(), &elementSubRegionList() );
  this->registerWrapper( viewKeyStruct::elementListString(), &elementList() );

  excludeWrappersFromPacking( { viewKeyStruct::edgeListString(),
                                viewKeyStruct::faceListString(),
                                viewKeyStruct::elementRegionListString(),
                                viewKeyStruct::elementSubRegionListString(),
                                viewKeyStruct::elementListString() } );
}


void NodeManager::resize( localIndex const newSize )
{
  m_toFacesRelation.resize( newSize, 2 * getFaceMapOverallocation() );
  m_toEdgesRelation.resize( newSize, 2 * getEdgeMapOverallocation() );
  m_toElements.m_toElementRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementSubRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementIndex.resize( newSize, 2 * getElemMapOverAllocation() );
  ObjectManagerBase::resize( newSize );
}


void NodeManager::constructGlobalToLocalMap( CellBlockManagerABC const & cellBlockManager )
{
  m_localToGlobalMap = cellBlockManager.getNodeLocalToGlobal();
  ObjectManagerBase::constructGlobalToLocalMap();
}

void NodeManager::buildSets( CellBlockManagerABC const & cellBlockManager,
                             GeometricObjectManager const & geometries )
{
  GEOS_MARK_FUNCTION;

  // Let's first copy the sets from the cell block manager.
  for( const auto & nameArray: cellBlockManager.getNodeSets() )
  {
    auto & array = m_sets.registerWrapper< SortedArray< localIndex > >( nameArray.first ).reference();
    array = nameArray.second;
  }

  // Now let's copy them from the geometric objects.
  buildGeometricSets( geometries );
}

void NodeManager::buildGeometricSets( GeometricObjectManager const & geometries )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = this->referencePosition();
  localIndex const numNodes = this->size();

  geometries.forSubGroups< SimpleGeometricObjectBase >( [&]( SimpleGeometricObjectBase const & object )
  {
    string const & name = object.getName();
    SortedArray< localIndex > & targetSet = m_sets.registerWrapper< SortedArray< localIndex > >( name ).reference();
    for( localIndex a = 0; a < numNodes; ++a )
    {
      real64 nodeCoord[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[a] );
      if( object.isCoordInObject( nodeCoord ) )
      {
        targetSet.insert( a );
      }
    }
  } );
}

void NodeManager::setDomainBoundaryObjects( FaceManager const & faceManager,
                                            EdgeManager const & edgeManager,
                                            ElementRegionManager const & elemRegionManager )
{
  arrayView1d< integer const > const isFaceOnDomainBoundary = faceManager.getDomainBoundaryIndicator();
  arrayView1d< integer const > const isEdgeOnDomainBoundary = edgeManager.getDomainBoundaryIndicator();
  arrayView1d< integer > const isNodeOnDomainBoundary = getDomainBoundaryIndicator();
  isNodeOnDomainBoundary.zero();

  // Here, we need to take into account the fact that the fractures elements
  // share the same nodes as the matrix.
  // So nodes can be considered as boundary objects while the face built on top of it may no.

  // We'll proceed as follows.
  //
  // Faces which are tagged as boundary faces will have their support nodes on the boundary as well.
  // (But keep in mind that some nodes may sit astride two faces with different statuses.)
  // Face nodes which are not on the boundary may however be tagged "boundary"
  // if they connect to a 2d element which is itself on the boundary (because those 2d elements share those nodes).
  //
  // To solve this issue, we loop over all the fracture 2d elements.
  // If they have only one neighboring matrix element,
  // this means that the fracture 2d element is on the boundary.
  // As such, we tag the support nodes of the unique neighboring face as boundary nodes.
  ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();

  forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const faceIndex )
  {
    if( isFaceOnDomainBoundary[faceIndex] == 1 )
    {
      for( localIndex const nodeIndex : faceToNodes[faceIndex] )
      {
        isNodeOnDomainBoundary[nodeIndex] = 1;
      }
    }
  } );

  auto const & edgeToNodes = edgeManager.nodeList().toViewConst();

  forAll< parallelHostPolicy >( edgeManager.size(), [=]( localIndex const edgeIndex )
  {
    if( isEdgeOnDomainBoundary[edgeIndex] == 1 )
    {
      for( localIndex const nodeIndex : edgeToNodes[edgeIndex] )
      {
        isNodeOnDomainBoundary[nodeIndex] = 1;
      }
    }
  } );

  auto const f = [&]( localIndex er,
                      SurfaceElementRegion const & region )
  {
    if( region.subRegionType() != SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      return;
    }

    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    ArrayOfArraysView< localIndex const > const elem2dToFaces = subRegion.faceList().toViewConst();
    for( int ei = 0; ei < elem2dToFaces.size(); ++ei )
    {
      if( elem2dToFaces.sizeOfArray( ei ) == 2 )
      { continue; }

      for( localIndex const & face: elem2dToFaces[ei] )
      {
        for( localIndex const & node: faceToNodes[face] )
        {
          isNodeOnDomainBoundary[node] = 1;
        }
      }
    }
  };
  elemRegionManager.forElementRegionsComplete< SurfaceElementRegion >( f );
}

void NodeManager::setGeometricalRelations( CellBlockManagerABC const & cellBlockManager,
                                           ElementRegionManager const & elemRegionManager )
{
  GEOS_MARK_FUNCTION;

  resize( cellBlockManager.numNodes() );

  m_referencePosition = cellBlockManager.getNodePositions();

  m_toEdgesRelation.base().assimilate< parallelHostPolicy >( cellBlockManager.getNodeToEdges(),
                                                             LvArray::sortedArrayManipulation::UNSORTED_NO_DUPLICATES );
  m_toFacesRelation.base().assimilate< parallelHostPolicy >( cellBlockManager.getNodeToFaces(),
                                                             LvArray::sortedArrayManipulation::UNSORTED_NO_DUPLICATES );

  ToCellRelation< ArrayOfArrays< localIndex > > const toCellBlock = cellBlockManager.getNodeToElements();
  array2d< localIndex > const blockToSubRegion = elemRegionManager.getCellBlockToSubRegionMap( cellBlockManager );
  meshMapUtilities::transformCellBlockToRegionMap< parallelHostPolicy >( blockToSubRegion.toViewConst(),
                                                                         toCellBlock,
                                                                         m_toElements );

  // TODO connect nodes to 2d elements!
  auto const connect2dElements = [&]( localIndex er,
                                      SurfaceElementRegion const & region )
  {
    if( region.subRegionType() != SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      return;
    }

    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    int const esr = 0;  // Since there's only on unique subregion, the index is always 0.
    auto const & elem2dToNodes = subRegion.nodeList();
    for( int ei = 0; ei < elem2dToNodes.size(); ++ei )
    {
      for( auto const ni: elem2dToNodes[ei] )
      {
        if( ni < 0 )
        { continue; }
        m_toElements.m_toElementRegion.emplaceBack( ni, er );
        m_toElements.m_toElementSubRegion.emplaceBack( ni, esr );
        m_toElements.m_toElementIndex.emplaceBack( ni, ei );
      }
    }
  };
  // Connecting all the 3d elements (information is already in the m_toElements mappings) and all the 2d elements.
  elemRegionManager.forElementRegionsComplete< SurfaceElementRegion >( connect2dElements );
  // TODO do the same for the EdgeManager
}

void NodeManager::setupRelatedObjectsInRelations( EdgeManager const & edgeManager,
                                                  FaceManager const & faceManager,
                                                  ElementRegionManager const & elementRegionManager )
{
  m_toEdgesRelation.setRelatedObject( edgeManager );
  m_toFacesRelation.setRelatedObject( faceManager );

  m_toElements.setElementRegionManager( elementRegionManager );
}


void NodeManager::compressRelationMaps()
{
  m_toEdgesRelation.compress();
  m_toFacesRelation.compress();
  m_toElements.m_toElementRegion.compress();
  m_toElements.m_toElementSubRegion.compress();
  m_toElements.m_toElementIndex.compress();
}


localIndex NodeManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}


localIndex NodeManager::packUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}


template< bool DO_PACKING >
localIndex NodeManager::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                            arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::edgeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toEdgesRelation.toArrayOfArraysView(),
                                               m_unmappedGlobalIndicesInToEdges,
                                               packList,
                                               this->localToGlobalMap(),
                                               m_toEdgesRelation.relatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toFacesRelation.toArrayOfArraysView(),
                                               m_unmappedGlobalIndicesInToFaces,
                                               packList,
                                               this->localToGlobalMap(),
                                               m_toFacesRelation.relatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::elementListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               this->m_toElements,
                                               packList,
                                               m_toElements.getElementRegionManager() );
  return packedSize;
}


localIndex NodeManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const )
{
  GEOS_MARK_FUNCTION;

  localIndex unPackedSize = 0;

  string temp;
  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ERROR_IF( temp != viewKeyStruct::edgeListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ERROR_IF( temp != viewKeyStruct::faceListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ERROR_IF( temp != viewKeyStruct::elementListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     this->m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}


void NodeManager::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toEdgesRelation,
                                    m_toEdgesRelation.relatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation,
                                    m_toFacesRelation.relatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );

}


void NodeManager::depopulateUpMaps( std::set< localIndex > const & receivedNodes,
                                    array2d< localIndex > const & edgesToNodes,
                                    ArrayOfArraysView< localIndex const > const & facesToNodes,
                                    ElementRegionManager const & elemRegionManager )
{

  ObjectManagerBase::cleanUpMap( receivedNodes, m_toEdgesRelation.toView(), edgesToNodes );
  ObjectManagerBase::cleanUpMap( receivedNodes, m_toFacesRelation.toView(), facesToNodes );

  for( auto const & targetIndex : receivedNodes )
  {
    std::set< std::tuple< localIndex, localIndex, localIndex > > eraseList;
    for( localIndex k=0; k<m_toElements.m_toElementRegion.sizeOfArray( targetIndex ); ++k )
    {
      localIndex const elemRegionIndex    = m_toElements.m_toElementRegion[targetIndex][k];
      localIndex const elemSubRegionIndex = m_toElements.m_toElementSubRegion[targetIndex][k];
      localIndex const elemIndex          = m_toElements.m_toElementIndex[targetIndex][k];

      CellElementSubRegion const & subRegion = elemRegionManager.getRegion( elemRegionIndex ).
                                                 getSubRegion< CellElementSubRegion >( elemSubRegionIndex );
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const downmap = subRegion.nodeList();
      bool hasTargetIndex = false;

      for( localIndex a=0; a<downmap.size( 1 ); ++a )
      {
        localIndex const compositeLocalIndex = downmap( elemIndex, a );
        if( compositeLocalIndex==targetIndex )
        {
          hasTargetIndex=true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert( std::make_tuple( elemRegionIndex, elemSubRegionIndex, elemIndex ) );
      }
    }
    for( auto const & val : eraseList )
    {
      erase( m_toElements, targetIndex, std::get< 0 >( val ), std::get< 1 >( val ), std::get< 2 >( val ) );
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, NodeManager, string const &, Group * const )

}
