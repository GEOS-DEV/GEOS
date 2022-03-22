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
#include "FaceManager.hpp"
#include "EdgeManager.hpp"
#include "ToElementRelation.hpp"
#include "BufferOps.hpp"
#include "common/TimingMacros.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
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

}


NodeManager::~NodeManager()
{}


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

/**
 * @brief Populates the node to element region and node to element subregion mappings.
 * @param [in] elementRegionMgr The ElementRegionManager associated with this mesh level. Regions and subregions come from it.
 * @param [in] n2e The node to element maps.
 * @param [in,out] n2er The face to element region map.
 * @param [in,out] n2esr The face to element subregion map.
 *
 * @warning The @p n2e, @p n2er and @p n2esr need to be allocated at the correct dimensions, but need to be empty.
 */
void populateRegions( ElementRegionManager const & elementRegionMgr,
                      ArrayOfArraysView< localIndex const > const & n2e,
                      ArrayOfArraysView< localIndex > const & n2er,
                      ArrayOfArraysView< localIndex > const & n2esr )
{
  GEOSX_ERROR_IF_NE( n2e.size(), n2er.size() );
  GEOSX_ERROR_IF_NE( n2e.size(), n2esr.size() );

  // There is an implicit convention in the (n2e, n2er, n2esr) triplet.
  // The node to element mapping (n2e) binds node indices to multiple element indices (like `n -> (e0, e1,...)`).
  // The node to regions (n2r) and sub-regions (n2sr) respectively bind
  // node indices to the regions/sub-regions: `n -> (er0, er1)` and `n -> (esr0, esr1)`.
  //
  // It is assumed in the code that triplets obtained at indices 0, 1,... of all these relations,
  // (respectively `(e0, er0, esr0)`, `(e1, er1, esr1)`,...) are consistent:
  // `e0` should belong to both `er0` and `esr0`.
  //
  // But in some configuration (multi-regions, multi-subregions, parallel computing, ghost cells),
  // it may happen that a same element index appear multiple times in a same entry of the node to elements mapping.
  // For example `n0 -> (e0, e1, e2, e0)` where `e0` appears twice.
  // These duplicated elements will in fact belong to multiple regions/sub-regions.
  // This implies a specific care when filling the nodes to regions and sub-regions mappings.
  //
  // Thus, the algorithm is the following.
  //
  // For each sub-region of each region, we consider each node of each element.
  // We list all the elements connected to this node.
  // Then we take the first element that has not already been inserted in the mappings
  // (we check if the value is still -1), and we insert it using the `er` and `esr` values,
  // before moving to another node.
  //
  // The same node index with the same elements will eventually come back during the iterative process.
  // But this time with another region/sub-region combination.
  //
  // It must be noted that we can take the duplicated elements in any order,
  // as long as the insertions are consistent!

  // The algorithm is equivalent to the algorithm described
  // in the `populateRegions` in the `FaceManager.cpp` file.
  // Instead of nodes, we'll have faces.
  //
  // Since the algorithm is quite short, and because of slight variations
  // (e.g. the different allocation between faces and nodes implementation),
  // I considered acceptable to duplicate it a bit.
  // This is surely disputable.

  // This function `f` will be applied on every sub-region.
  auto f = [&n2e, &n2er, &n2esr]( localIndex er,
                                  localIndex esr,
                                  ElementRegionBase const &,
                                  CellElementSubRegion const & subRegion ) -> void
  {
    for( localIndex iElement = 0; iElement < subRegion.size(); ++iElement )
    {
      for( localIndex iNodeLoc = 0; iNodeLoc < subRegion.numNodesPerElement(); ++iNodeLoc )
      {
        // iNodeLoc is the node index in the referential of each cell (0 to 7 for a cube, e.g.).
        // While iNode is the global index of the node.
        localIndex const & iNode = subRegion.nodeList( iElement, iNodeLoc );

        // A node may be attached to `numElementsLoc` elements.
        localIndex const numElementsLoc = n2e[iNode].size();
        for( localIndex iElementLoc = 0; iElementLoc < numElementsLoc; ++iElementLoc )
        {
          // We only consider the elements that match the mapping.
          if( n2e( iNode, iElementLoc ) != iElement )
          {
            continue;
          }

          // This loop is a small hack to back insert/allocate dummy elements
          // that will eventually be replaced by the correct values.
          for( localIndex i = n2er[iNode].size(); i < iElementLoc + 1; ++i )
          {
            // By construction n2er and n2esr have the same size, so we use the same loop.
            n2er.emplaceBack( iNode, -1 );
            n2esr.emplaceBack( iNode, -1 );
          }

          // Here we fill the mapping iff it has not already been inserted.
          if( n2er( iNode, iElementLoc ) < 0 or n2esr( iNode, iElementLoc ) < 0 )
          {
            n2er( iNode, iElementLoc ) = er;
            n2esr( iNode, iElementLoc ) = esr;

            // We only want to insert one unique index that has not been inserted,
            // so we quit the loop on indices here.
            break;
          }
        }
      }
    }
  };

  elementRegionMgr.forElementSubRegionsComplete< CellElementSubRegion >( f );
}

void NodeManager::buildRegionMaps( ElementRegionManager const & elementRegionManager )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArrays< localIndex > & toElementRegionList = m_toElements.m_toElementRegion;
  ArrayOfArrays< localIndex > & toElementSubRegionList = m_toElements.m_toElementSubRegion;
  ArrayOfArraysView< localIndex const > const & toElementList = m_toElements.m_toElementIndex.toViewConst();
  localIndex const numNodes = size();

  // Resize the node to elem map.
  toElementRegionList.resize( 0 );
  toElementSubRegionList.resize( 0 );

  // TODO Since there is a strong requirement that `toElementRegionList` and `toElementSubRegionList`
  //      share the same capacities as `toElementList`, then it's should be interesting
  //      to retrieve this information from `toElementList` instead.

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  toElementRegionList.reserve( entriesToReserve );
  toElementSubRegionList.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = numNodes + numNodes * getElemMapOverAllocation() * ( 1 + 2 * overAllocationFactor );
  toElementRegionList.reserveValues( valuesToReserve );
  toElementSubRegionList.reserveValues( valuesToReserve );

  // Append an array for each node with capacity to hold the appropriate number of elements plus some wiggle room.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    toElementRegionList.appendArray( 0 );
    toElementSubRegionList.appendArray( 0 );

    const localIndex numElementsPerNode = toElementList[nodeID].size() + getElemMapOverAllocation();
    toElementRegionList.setCapacityOfArray( nodeID, numElementsPerNode );
    toElementSubRegionList.setCapacityOfArray( nodeID, numElementsPerNode );
  }

  // Delegate the computation to a free function.
  populateRegions( elementRegionManager,
                   toElementList,
                   toElementRegionList.toView(),
                   toElementSubRegionList.toView() );
}

void NodeManager::buildSets( CellBlockManagerABC const & cellBlockManager,
                             GeometricObjectManager const & geometries )
{
  // Let's first copy the sets from the cell block manager.
  for( const auto & nameArray: cellBlockManager.getNodeSets() )
  {
    auto & array = m_sets.registerWrapper< SortedArray< localIndex > >( nameArray.first ).reference();
    array = nameArray.second;
  }

  // Now let's copy them from the geometric objects.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = this->referencePosition();
  localIndex const numNodes = this->size();

  geometries.forSubGroups< SimpleGeometricObjectBase >(
    [&]( SimpleGeometricObjectBase const & object ) -> void
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

void NodeManager::setDomainBoundaryObjects( FaceManager const & faceManager )
{
  arrayView1d< integer const > const & isFaceOnDomainBoundary = faceManager.getDomainBoundaryIndicator();
  arrayView1d< integer > const & isNodeOnDomainBoundary = getDomainBoundaryIndicator();
  isNodeOnDomainBoundary.zero();

  ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();

  forAll< parallelHostPolicy >( faceManager.size(), [&]( localIndex const k )
  {
    if( isFaceOnDomainBoundary[k] == 1 )
    {
      localIndex const numNodes = faceToNodes.sizeOfArray( k );
      for( localIndex a = 0; a < numNodes; ++a )
      {
        isNodeOnDomainBoundary[faceToNodes( k, a )] = 1;
      }
    }
  } );
}

void NodeManager::setGeometricalRelations( CellBlockManagerABC const & cellBlockManager )
{
  resize( cellBlockManager.numNodes() );

  m_referencePosition = cellBlockManager.getNodesPositions();

  m_toEdgesRelation.base() = cellBlockManager.getNodeToEdges();
  m_toFacesRelation.base() = cellBlockManager.getNodeToFaces();

  m_toElements.m_toElementIndex = cellBlockManager.getNodeToElements();
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


std::set< string > NodeManager::getPackingExclusionList() const
{
  std::set< string > result = ObjectManagerBase::getPackingExclusionList();
  result.insert( { viewKeyStruct::edgeListString(),
                   viewKeyStruct::faceListString(),
                   viewKeyStruct::elementRegionListString(),
                   viewKeyStruct::elementSubRegionListString(),
                   viewKeyStruct::elementListString() } );

  if( this->hasWrapper( "usedFaces" ) )
  {
    result.insert( "usedFaces" );
  }
  return result;
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
  localIndex unPackedSize = 0;

  string temp;
  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::edgeListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::faceListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::elementListString(), "" );
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
