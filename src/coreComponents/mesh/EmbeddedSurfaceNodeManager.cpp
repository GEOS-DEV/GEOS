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
 * @file EmbeddedSurfaceNodeManager.cpp
 */

#include "EmbeddedSurfaceNodeManager.hpp"
#include "FaceManager.hpp"
#include "EdgeManager.hpp"
#include "ToElementRelation.hpp"
#include "BufferOps.hpp"
#include "mesh/ExtrinsicMeshData.hpp"
#include "common/TimingMacros.hpp"
#include "common/MpiWrapper.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{

using namespace dataRepository;

EmbeddedSurfaceNodeManager::EmbeddedSurfaceNodeManager( string const & name,
                                                        Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition( 0, 3 )
{
  registerWrapper( viewKeyStruct::referencePositionString(), &m_referencePosition );
  this->registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );
  this->registerWrapper( viewKeyStruct::elementRegionListString(), &elementRegionList() );
  this->registerWrapper( viewKeyStruct::elementSubRegionListString(), &elementSubRegionList() );
  this->registerWrapper( viewKeyStruct::elementListString(), &elementList() );
  this->registerWrapper( viewKeyStruct::parentEdgeGlobalIndexString(), &m_parentEdgeGlobalIndex );
}


EmbeddedSurfaceNodeManager::~EmbeddedSurfaceNodeManager()
{}


void EmbeddedSurfaceNodeManager::resize( localIndex const newSize )
{
  m_toEdgesRelation.resize( newSize, 2 * getEdgeMapOverallocation() );
  m_toElements.m_toElementRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementSubRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementIndex.resize( newSize, 2 * getElemMapOverAllocation() );
  ObjectManagerBase::resize( newSize );
}


void EmbeddedSurfaceNodeManager::setEdgeMaps( EdgeManager const & embSurfEdgeManager )
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< localIndex const > const edgeToNodeMap = embSurfEdgeManager.nodeList();
  localIndex const numEdges = edgeToNodeMap.size( 0 );
  localIndex const numNodes = size();

  ArrayOfArrays< localIndex > toEdgesTemp( numNodes, EdgeManager::maxEdgesPerNode() );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeEdges = 0;

  forAll< parallelHostPolicy >( numEdges, [&]( localIndex const edgeID )
  {
    toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( edgeToNodeMap( edgeID, 0 ), edgeID );
    toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( edgeToNodeMap( edgeID, 1 ), edgeID );
    totalNodeEdges += 2;
  } );

  // Resize the node to edge map.
  m_toEdgesRelation.resize( 0 );

  // Reserve space for the number of current nodes plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  m_toEdgesRelation.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeEdges.get() + numNodes * getEdgeMapOverallocation() * ( 1 + 2 * overAllocationFactor );
  m_toEdgesRelation.reserveValues( valuesToReserve );

  // Append the individual sets.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    m_toEdgesRelation.appendSet( toEdgesTemp.sizeOfArray( nodeID ) + getEdgeMapOverallocation() );
  }

  ArrayOfSetsView< localIndex > const & toEdgesView = m_toEdgesRelation.toView();
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex * const edges = toEdgesTemp[ nodeID ];
    localIndex const numNodeEdges = toEdgesTemp.sizeOfArray( nodeID );
    localIndex const numUniqueEdges = LvArray::sortedArrayManipulation::makeSortedUnique( edges, edges + numNodeEdges );
    toEdgesView.insertIntoSet( nodeID, edges, edges + numUniqueEdges );
  } );

  m_toEdgesRelation.setRelatedObject( embSurfEdgeManager );
}

void EmbeddedSurfaceNodeManager::setElementMaps( ElementRegionManager const & elementRegionManager )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArrays< localIndex > & toElementRegionList = m_toElements.m_toElementRegion;
  ArrayOfArrays< localIndex > & toElementSubRegionList = m_toElements.m_toElementSubRegion;
  ArrayOfArrays< localIndex > & toElementList = m_toElements.m_toElementIndex;
  localIndex const numNodes = size();

  // The number of elements attached to the each node.
  array1d< localIndex > elemsPerNode( numNodes );

  // The total number of elements, the sum of elemsPerNode.
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeElems = 0;

  elementRegionManager.
    forElementSubRegions< EmbeddedSurfaceSubRegion >( [&elemsPerNode, &totalNodeElems]( EmbeddedSurfaceSubRegion const & subRegion )
  {
    EmbeddedSurfaceSubRegion::NodeMapType const & elemToNodeMap = subRegion.nodeList();
    forAll< parallelHostPolicy >( subRegion.size(), [&elemsPerNode, totalNodeElems, &elemToNodeMap ] ( localIndex const k )
    {
      localIndex const numElementNodes = elemToNodeMap.sizeOfArray( k );
      totalNodeElems += numElementNodes;
      for( localIndex a = 0; a < numElementNodes; ++a )
      {
        localIndex const nodeIndex = elemToNodeMap( k, a );
        RAJA::atomicInc< parallelHostAtomic >( &elemsPerNode[ nodeIndex ] );
      }
    } );
  } );

  // Resize the node to elem map.
  toElementRegionList.resize( 0 );
  toElementSubRegionList.resize( 0 );
  toElementList.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  toElementRegionList.reserve( entriesToReserve );
  toElementSubRegionList.reserve( entriesToReserve );
  toElementList.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeElems.get() + numNodes * getElemMapOverAllocation() * ( 1 + 2 * overAllocationFactor );
  toElementRegionList.reserveValues( valuesToReserve );
  toElementSubRegionList.reserveValues( valuesToReserve );
  toElementList.reserveValues( valuesToReserve );

  // Append an array for each node with capacity to hold the appropriate number of elements plus some wiggle room.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    toElementRegionList.appendArray( 0 );
    toElementSubRegionList.appendArray( 0 );
    toElementList.appendArray( 0 );

    toElementRegionList.setCapacityOfArray( nodeID, elemsPerNode[ nodeID ] + getElemMapOverAllocation() );
    toElementSubRegionList.setCapacityOfArray( nodeID, elemsPerNode[ nodeID ] + getElemMapOverAllocation() );
    toElementList.setCapacityOfArray( nodeID, elemsPerNode[ nodeID ] + getElemMapOverAllocation() );
  }

  // Populate the element maps.
  // Note that this can't be done in parallel because the three element lists must be in the same order.
  // If this becomes a bottleneck create a temporary ArrayOfArrays of tuples and insert into that first then copy over.
  elementRegionManager.
    forElementSubRegionsComplete< EmbeddedSurfaceSubRegion >( [&toElementRegionList, &toElementSubRegionList, &toElementList]
                                                                ( localIndex const er, localIndex const esr, ElementRegionBase const &,
                                                                EmbeddedSurfaceSubRegion const & subRegion )
  {
    EmbeddedSurfaceSubRegion::NodeMapType const & elemToNodeMap = subRegion.nodeList();
    for( localIndex k = 0; k < subRegion.size(); ++k )
    {
      for( localIndex a=0; a<elemToNodeMap.sizeOfArray( k ); ++a )
      {
        localIndex const nodeIndex = elemToNodeMap( k, a );
        toElementRegionList.emplaceBack( nodeIndex, er );
        toElementSubRegionList.emplaceBack( nodeIndex, esr );
        toElementList.emplaceBack( nodeIndex, k );
      }
    }
  } );

  this->m_toElements.setElementRegionManager( elementRegionManager );
}


void EmbeddedSurfaceNodeManager::compressRelationMaps()
{
  m_toEdgesRelation.compress();
  m_toElements.m_toElementRegion.compress();
  m_toElements.m_toElementSubRegion.compress();
  m_toElements.m_toElementIndex.compress();
}


void EmbeddedSurfaceNodeManager::appendNode( arraySlice1d< real64 const > const & pointCoord,
                                             integer const & pointGhostRank )
{
  localIndex nodeIndex =  this->size();
  this->resize( nodeIndex + 1 );
  LvArray::tensorOps::copy< 3 >( m_referencePosition[nodeIndex], pointCoord );
  m_ghostRank[ nodeIndex ] = pointGhostRank;
}

void EmbeddedSurfaceNodeManager::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::viewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementRegionListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementSubRegionListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementListString() ));
}


localIndex EmbeddedSurfaceNodeManager::packNewNodesGlobalMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packNewNodesGlobalMapsPrivate< false >( junk, packList );
}

localIndex EmbeddedSurfaceNodeManager::packNewNodesGlobalMaps( buffer_unit_type * & buffer,
                                                               arrayView1d< localIndex const > const & packList ) const
{
  return packNewNodesGlobalMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex EmbeddedSurfaceNodeManager::packNewNodesGlobalMapsPrivate( buffer_unit_type * & buffer,
                                                                      arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = bufferOps::Pack< DOPACK >( buffer, this->getName() );

  // this doesn't link without the string()...no idea why.
  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::localToGlobalMapString() ) );

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  packedSize += bufferOps::Pack< DOPACK >( buffer, rank );

  localIndex const numPackedIndices = packList.size();
  packedSize += bufferOps::Pack< DOPACK >( buffer, numPackedIndices );

  if( numPackedIndices > 0 )
  {
    // We pack 3 things:
    // 1. the global indices
    globalIndex_array globalIndices;
    globalIndices.resize( numPackedIndices );
    // 2. the ghostRank
    array1d< integer > ghostRanks;
    ghostRanks.resize( numPackedIndices );
    for( localIndex a=0; a<numPackedIndices; ++a )
    {
      globalIndices[a] = this->m_localToGlobalMap[packList[a]];
      ghostRanks[a]= this->m_ghostRank[ packList[a] ];
    }
    packedSize += bufferOps::Pack< DOPACK >( buffer, globalIndices );
    packedSize += bufferOps::Pack< DOPACK >( buffer, ghostRanks );

    // 3. the referencePosition
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & referencePosition = this->referencePosition();
    packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::referencePositionString() ) );
    packedSize += bufferOps::PackByIndex< DOPACK >( buffer, referencePosition, packList );
  }

  return packedSize;
}


localIndex EmbeddedSurfaceNodeManager::unpackNewNodesGlobalMaps( buffer_unit_type const * & buffer,
                                                                 localIndex_array & packList )
{
  GEOSX_MARK_FUNCTION;

  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOSX_ERROR_IF( groupName != this->getName(), "EmbeddedSurfaceNodeManager::unpackGlobalMaps(): group names do not match" );

  string localToGlobalString;
  unpackedSize += bufferOps::Unpack( buffer, localToGlobalString );
  GEOSX_ERROR_IF( localToGlobalString != viewKeyStruct::localToGlobalMapString(), "ObjectManagerBase::Unpack(): label incorrect" );

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int sendingRank;
  unpackedSize += bufferOps::Unpack( buffer, sendingRank );

  localIndex numUnpackedIndices;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedIndices );

  if( numUnpackedIndices > 0 )
  {
    localIndex_array unpackedLocalIndices;
    unpackedLocalIndices.resize( numUnpackedIndices );

    globalIndex_array globalIndices;
    unpackedSize += bufferOps::Unpack( buffer, globalIndices );

    array1d< integer > ghostRankOnSendingRank;
    unpackedSize += bufferOps::Unpack( buffer, ghostRankOnSendingRank );

    // Unpack referencePosition
    localIndex_array indicesOnBuffer( numUnpackedIndices );
    for( localIndex i = 0; i < numUnpackedIndices; i++ )
    {
      indicesOnBuffer[i] = i;
    }
    array2d< real64, nodes::REFERENCE_POSITION_PERM > referencePositionData( 0, 3 );
    referencePositionData.resize( numUnpackedIndices );
    arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & referencePosition = referencePositionData.toView();
    string referencePositionString;
    unpackedSize += bufferOps::Unpack( buffer, referencePositionString );
    GEOSX_ERROR_IF( referencePositionString != viewKeyStruct::referencePositionString(), "EmbeddedSurfaceNodeManager::unpackGlobalMaps(): label incorrect" );
    unpackedSize += bufferOps::UnpackByIndex( buffer, referencePosition, indicesOnBuffer );

    localIndex numNewIndices = 0;
    globalIndex_array newGlobalIndices;
    newGlobalIndices.reserve( numUnpackedIndices );
    array1d< integer > ghostRank;
    ghostRank.reserve( numUnpackedIndices );
    localIndex const oldSize = this->size();
    for( localIndex a = 0; a < numUnpackedIndices; ++a )
    {
      // check to see if the node already exists by comparing the coordinates by checking for the global
      localIndex nodeIndexOnThisRank = -1;
      real64 nodeCoord[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( referencePosition[a] );
      nodeExistsOnThisRank( nodeCoord, nodeIndexOnThisRank );
      if( nodeIndexOnThisRank > -1 )
      {
        if( m_ghostRank[nodeIndexOnThisRank] == sendingRank )
        {
          // object already exists on this rank and it's a ghost
          unpackedLocalIndices( a ) = nodeIndexOnThisRank;

          globalIndex const globalIndexOnthisRank = m_localToGlobalMap[ nodeIndexOnThisRank ];

          m_localToGlobalMap[ nodeIndexOnThisRank ] = globalIndices( a );

          //modify global to local map
          m_globalToLocalMap.erase( globalIndexOnthisRank );
          m_globalToLocalMap.insert( {globalIndices( a ), nodeIndexOnThisRank} );

        }
      }
      else
      {
        // object does not exist on this domain
        const localIndex newLocalIndex = oldSize + numNewIndices;

        // add the global index of the new object to the globalToLocal map
        m_globalToLocalMap[ globalIndices[a] ] = newLocalIndex;

        unpackedLocalIndices( a ) = newLocalIndex;

        newGlobalIndices.emplace_back( globalIndices[a] );

        if( ghostRankOnSendingRank[a] == rank )
        {
          ghostRank.emplace_back( -1 );
        }
        else
        {
          ghostRank.emplace_back( rank );
        }

        ++numNewIndices;

        GEOSX_ERROR_IF( packList.size() != 0,
                        "EmbeddedSurfaceNodeManager::unpackGlobalMaps(): packList specified, "
                        "but a new globalIndex is unpacked" );
      }
    }

    ///

    // figure out new size of object container, and resize it
    const localIndex newSize = oldSize + numNewIndices;
    this->resize( newSize );

    // add the new indices to the maps.
    for( int a=0; a<numNewIndices; ++a )
    {
      localIndex const b = oldSize + a;
      m_localToGlobalMap[b] = newGlobalIndices( a );
      m_ghostRank[b] = ghostRank( a );
    }


    packList = unpackedLocalIndices;
  }

  return unpackedSize;
}


localIndex EmbeddedSurfaceNodeManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}


localIndex EmbeddedSurfaceNodeManager::packUpDownMaps( buffer_unit_type * & buffer,
                                                       arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}


template< bool DOPACK >
localIndex EmbeddedSurfaceNodeManager::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                              arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::elementListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           this->m_toElements,
                                           packList,
                                           m_toElements.getElementRegionManager() );
  return packedSize;
}


localIndex EmbeddedSurfaceNodeManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                         localIndex_array & packList,
                                                         bool const overwriteUpMaps,
                                                         bool const )
{
  localIndex unPackedSize = 0;

  string temp;

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::elementListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     this->m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void EmbeddedSurfaceNodeManager::nodeExistsOnThisRank( real64 const (&nodeCoord)[3],
                                                       localIndex & nodeIndexOnThisRank ) const
{
  for( localIndex a=0; a<size(); a++ )
  {
    real64 const tolerance = 100.0 * std::numeric_limits< real64 >::epsilon();
    real64 distance[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_referencePosition[a] );
    LvArray::tensorOps::subtract< 3 >( distance, nodeCoord );
    if( LvArray::tensorOps::l2Norm< 3 >( distance ) < tolerance )
    {
      nodeIndexOnThisRank = a;
    }
  }
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, EmbeddedSurfaceNodeManager, string const &, Group * const )

}
