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
 * @file EdgeManager.cpp
 */

#include "EdgeManager.hpp"

#include "BufferOps.hpp"
#include "NodeManager.hpp"
#include "FaceManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;

EdgeManager::EdgeManager( string const & name,
                          Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_edgesToFractureConnectorsEdges(),
  m_fractureConnectorsEdgesToEdges(),
  m_fractureConnectorEdgesToFaceElements()
{
  this->registerWrapper( viewKeyStruct::nodeListString(), &this->m_toNodesRelation );
  this->registerWrapper( viewKeyStruct::faceListString(), &this->m_toFacesRelation );

  m_toNodesRelation.resize( 0, 2 );

  registerWrapper( viewKeyStruct::edgesTofractureConnectorsEdgesString(), &m_edgesToFractureConnectorsEdges ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of edge local indices to the fracture connector local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorEdgesToEdgesString(), &m_fractureConnectorsEdgesToEdges ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices to edge local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorsEdgesToFaceElementsIndexString(),
                   &m_fractureConnectorEdgesToFaceElements ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices face element local indices" ).
    setSizedFromParent( 0 );

}

EdgeManager::~EdgeManager()
{}

void EdgeManager::resize( localIndex const newSize )
{
  m_toFacesRelation.resize( newSize, 2 * faceMapExtraSpacePerEdge() );
  ObjectManagerBase::resize( newSize );
}

/**
 * @class EdgeBuilder
 * @brief This class stores the data necessary to construct the various edge maps.
 */
struct EdgeBuilder
{

  /**
   * @brief Constructor.
   * @param [in] n1_ the greater of the two node indices that comprise the edge.
   * @param [in] faceID_ the ID of the face this edge came from.
   * @param [in] faceLocalEdgeIndex_ the face local index of this edge.
   */
  EdgeBuilder( localIndex const n1_,
               localIndex const faceID_,
               localIndex const faceLocalEdgeIndex_ ):
    n1( int32_t( n1_ ) ),
    faceID( int32_t( faceID_ ) ),
    faceLocalEdgeIndex( int32_t( faceLocalEdgeIndex_ ) )
  {}

  /**
   * @brief Imposes an ordering on EdgeBuilders. First compares n1 and then the faceID.
   * @param [in] rhs the EdgeBuilder to compare against.
   */
  bool operator<( EdgeBuilder const & rhs ) const
  {
    if( n1 < rhs.n1 ) return true;
    if( n1 > rhs.n1 ) return false;
    return faceID < rhs.faceID;
  }

  /**
   * @brief Return true if the two EdgeBuilders share the same greatest node index.
   * @param [in] rhs the EdgeBuilder to compare against.
   */
  bool operator==( EdgeBuilder const & rhs ) const
  { return n1 == rhs.n1; }

  /**
   * @brief Return true if the two EdgeBuilders don't share the same greatest node index.
   * @param [in] rhs the EdgeBuilder to compare against.
   */
  bool operator!=( EdgeBuilder const & rhs ) const
  { return n1 != rhs.n1; }

  int32_t n1;                  // The larger of the two node indices that comprise the edge.
  int32_t faceID;              // The face the edge came from.
  int32_t faceLocalEdgeIndex;  // The face local index of the edge.
};

/**
 * @brief Populate the edgesByLowestNode map.
 * @param [in] faceToNodeMap a map that associates an ordered list of nodes with each face.
 * @param [in/out] edgesByLowestNode of size numNodes, where each sub array has been preallocated to hold
 *        *enough* space.
 * For each edge of each face, this function gets the lowest node in the edge n0, creates an EdgeBuilder
 * associated with the edge and then appends the EdgeBuilder to edgesByLowestNode[ n0 ]. Finally it sorts
 * the contents of each sub-array of edgesByLowestNode from least to greatest.
 */
void createEdgesByLowestNode( ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                              ArrayOfArraysView< EdgeBuilder > const & edgesByLowestNode )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numFaces = faceToNodeMap.size();

  // loop over all the faces.
  forAll< parallelHostPolicy >( numFaces, [&]( localIndex const faceID )
  {
    localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( faceID );

    // loop over all the nodes in the face. there will be an edge for each node.
    for( localIndex a=0; a< numNodesInFace; ++a )
    {
      // sort the nodes in order of index value.
      localIndex node0 = faceToNodeMap( faceID, a );
      localIndex node1 = faceToNodeMap( faceID, ( a + 1 ) % numNodesInFace );
      if( node0 > node1 )
        std::swap( node0, node1 );

      // And append the edge to edgesByLowestNode.
      edgesByLowestNode.emplaceBackAtomic< parallelHostAtomic >( node0, node1, faceID, a );
    }
  } );

  // Loop over all the nodes and sort the associated edges.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    EdgeBuilder * const edges = edgesByLowestNode[ nodeID ];
    std::sort( edges, edges + edgesByLowestNode.sizeOfArray( nodeID ) );
  } );
}

/**
 * @brief Return the total number of unique edges and fill in the uniqueEdgeOffsets array.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [out] uniqueEdgeOffsets an array of size numNodes + 1. After this function returns node i contains
 * edges with IDs ranging from uniqueEdgeOffsets[ i ] to uniqueEdgeOffsets[ i + 1 ] - 1.
 */
localIndex calculateTotalNumberOfEdges( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                                        arrayView1d< localIndex > const & uniqueEdgeOffsets )
{
  localIndex const numNodes = edgesByLowestNode.size();
  GEOSX_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );

  uniqueEdgeOffsets[0] = 0;

  // Loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // If there are no edges associated with this node we can skip it.
    if( numEdges == 0 )
      return;

    localIndex & numUniqueEdges = uniqueEdgeOffsets[ nodeID + 1 ];
    numUniqueEdges = 0;

    // Otherwise since edgesByLowestNode[ nodeID ] is sorted we can compare subsequent entries
    // count up the unique entries.
    localIndex j = 0;
    for(; j < numEdges - 1; ++j )
    {
      numUniqueEdges += edgesByLowestNode( nodeID, j ) != edgesByLowestNode( nodeID, j + 1 );
    }

    numUniqueEdges += j == numEdges - 1;
  } );

  // At this point uniqueEdgeOffsets[ i ] holds the number of unique edges associated with node i - 1.
  // Perform an inplace prefix-sum to get the unique edge offset.
  RAJA::inclusive_scan_inplace< parallelHostPolicy >( uniqueEdgeOffsets.begin(), uniqueEdgeOffsets.end() );

  return uniqueEdgeOffsets.back();
}

/**
 * @brief Resize the edge to face map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an containing the unique edge IDs for each node in edgesByLowestNode.
 * param [out] edgeToFaceMap the map from edges to faces. This function resizes the array appropriately.
 */
void resizeEdgeToFaceMap( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                          arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                          ArrayOfSets< localIndex > & edgeToFaceMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  array1d< localIndex > numFacesPerEdge( numUniqueEdges );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalEdgeFaces( 0.0 );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeID ];
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // loop over all the EdgeBuilders associated with the node
    localIndex j = 0;
    while( j < numEdges - 1 )
    {
      // Find the number of EdgeBuilders that describe the same edge
      localIndex numMatches = 1;
      while( edgesByLowestNode( nodeID, j ) == edgesByLowestNode( nodeID, j + numMatches ) )
      {
        ++numMatches;
        if( j + numMatches == numEdges )
          break;
      }

      // The number of matches is the number of faces associated with this edge.
      numFacesPerEdge( curEdgeID ) = numMatches;
      totalEdgeFaces += numFacesPerEdge( curEdgeID );
      ++curEdgeID;
      j += numMatches;
    }

    if( j == numEdges - 1 )
    {
      numFacesPerEdge( curEdgeID ) = 1;
      totalEdgeFaces += numFacesPerEdge( curEdgeID );
    }
  } );

  // Resize the edge to face map
  edgeToFaceMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numUniqueEdges;
  edgeToFaceMap.reserve( entriesToReserve );

  // Reserve space for the total number of edge faces + extra space for existing edges + even more space for new edges.
  localIndex const valuesToReserve = totalEdgeFaces.get() + numUniqueEdges * EdgeManager::faceMapExtraSpacePerEdge() * ( 1 + 2 * overAllocationFactor );
  edgeToFaceMap.reserveValues( valuesToReserve );

  // Append the individual sets.
  for( localIndex faceID = 0; faceID < numUniqueEdges; ++faceID )
  {
    edgeToFaceMap.appendSet( numFacesPerEdge[ faceID ] + EdgeManager::faceMapExtraSpacePerEdge() );
  }
}


/**
 * @brief Add an edge to the face to edge map, edge to face map, and edge to node map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in/out] faceToEdgeMap the map from face IDs to edge IDs.
 * @param [in/out] edgeToFacemap the map from edgeIDs to faceIDs.
 * @param [in/out] edgeToNodeMap the map from edgeIDs to nodeIDs.
 * @param [in] edgeID the ID of the edge to add.
 * @param [in] firstNodeID the ID of the first node of the edge.
 * @param [in] firstMatch the index of the first EdgeBuilder that describes this edge in edgesByLowestNode[ firstNodeID
 *].
 * @param [in] numMatches the number of EdgeBuilders that describe this edge in edgesByLowestNode[ firstNodeID ].
 */
void addEdge( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
              ArrayOfArraysView< localIndex > const & faceToEdgeMap,
              ArrayOfSetsView< localIndex > const & edgeToFaceMap,
              arrayView2d< localIndex > const & edgeToNodeMap,
              localIndex const edgeID,
              localIndex const firstNodeID,
              localIndex const firstMatch,
              localIndex const numMatches )
{
  GEOSX_ASSERT_GE( edgeToFaceMap.capacityOfSet( edgeID ), numMatches );

  // Populate the edge to node map.
  edgeToNodeMap( edgeID, 0 ) = firstNodeID;
  edgeToNodeMap( edgeID, 1 ) = edgesByLowestNode( firstNodeID, firstMatch ).n1;

  // Loop through all the matches and fill in the face to edge and edge to face maps.
  for( localIndex i = 0; i < numMatches; ++i )
  {
    localIndex const faceID = edgesByLowestNode( firstNodeID, firstMatch + i ).faceID;
    localIndex const faceLocalEdgeIndex = edgesByLowestNode( firstNodeID, firstMatch + i ).faceLocalEdgeIndex;

    faceToEdgeMap( faceID, faceLocalEdgeIndex ) = edgeID;
    edgeToFaceMap.insertIntoSet( edgeID, faceID );
  }
}

/**
 * @brief Populate the face to edge map, edge to face map, and edge to node map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an array containing the unique ID of the first edge associated with each node.
 * @param [in] faceToNodeMap the map from faces to nodes.
 * @param [in/out] faceToEdgeMap the map from face IDs to edge IDs.
 * @param [in/out] edgeToFacemap the map from edgeIDs to faceIDs.
 * @param [in/out] edgeToNodeMap the map from edgeIDs to nodeIDs.
 */
void populateMaps( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                   arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                   ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                   ArrayOfArrays< localIndex > & faceToEdgeMap,
                   ArrayOfSets< localIndex > & edgeToFaceMap,
                   arrayView2d< localIndex > const & edgeToNodeMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numFaces = faceToNodeMap.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  GEOSX_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );
  GEOSX_ERROR_IF_NE( numFaces, faceToEdgeMap.size() );
  GEOSX_ERROR_IF_NE( numUniqueEdges, edgeToFaceMap.size() );
  GEOSX_ERROR_IF_NE( numUniqueEdges, edgeToNodeMap.size( 0 ) );

  // The face to edge map has the same shape as the face to node map, so we can resize appropriately.
  localIndex totalSize = 0;
  for( localIndex faceID = 0; faceID < numFaces; ++faceID )
  {
    totalSize += faceToNodeMap.sizeOfArray( faceID );
  }

  // Resize the face to edge map
  faceToEdgeMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numFaces;
  faceToEdgeMap.reserve( entriesToReserve );

  // Reserve space for the total number of face edges + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalSize + numFaces * FaceManager::edgeMapExtraSpacePerFace() * ( 1 + 2 * overAllocationFactor );
  faceToEdgeMap.reserveValues( valuesToReserve );
  for( localIndex faceID = 0; faceID < numFaces; ++faceID )
  {
    faceToEdgeMap.appendArray( faceToNodeMap.sizeOfArray( faceID ) );
    faceToEdgeMap.setCapacityOfArray( faceToEdgeMap.size() - 1,
                                      faceToNodeMap.sizeOfArray( faceID ) + FaceManager::edgeMapExtraSpacePerFace() );
  }

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeID ];
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // loop over all the EdgeBuilders associated with the node
    localIndex j = 0;
    while( j < numEdges - 1 )
    {
      // Find the number of EdgeBuilders that describe the same edge
      localIndex numMatches = 1;
      while( edgesByLowestNode( nodeID, j ) == edgesByLowestNode( nodeID, j + numMatches ) )
      {
        ++numMatches;
        if( j + numMatches == numEdges )
          break;
      }
      // Then add the edge.
      addEdge( edgesByLowestNode, faceToEdgeMap.toView(), edgeToFaceMap.toView(), edgeToNodeMap, curEdgeID, nodeID, j, numMatches );
      ++curEdgeID;
      j += numMatches;
    }

    if( j == numEdges - 1 )
    {
      addEdge( edgesByLowestNode, faceToEdgeMap.toView(), edgeToFaceMap.toView(), edgeToNodeMap, curEdgeID, nodeID, j, 1 );
    }
  } );
}

void EdgeManager::buildEdges( NodeManager & nodeManager, FaceManager & faceManager )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = nodeManager.size();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  faceManager.edgeList().setRelatedObject( *this );
  ArrayOfArrays< localIndex > & faceToEdgeMap = faceManager.edgeList();

  m_toNodesRelation.setRelatedObject( nodeManager );
  m_toFacesRelation.setRelatedObject( faceManager );

  ArrayOfArrays< EdgeBuilder > edgesByLowestNode( numNodes, 2 * maxEdgesPerNode() );
  createEdgesByLowestNode( faceToNodeMap, edgesByLowestNode.toView() );

  array1d< localIndex > uniqueEdgeOffsets( numNodes + 1 );
  localIndex const numEdges = calculateTotalNumberOfEdges( edgesByLowestNode.toViewConst(), uniqueEdgeOffsets );

  resizeEdgeToFaceMap( edgesByLowestNode.toViewConst(),
                       uniqueEdgeOffsets,
                       m_toFacesRelation );

  resize( numEdges );

  populateMaps( edgesByLowestNode.toViewConst(),
                uniqueEdgeOffsets,
                faceToNodeMap,
                faceToEdgeMap,
                m_toFacesRelation,
                m_toNodesRelation );

  // make sets from nodesets
  auto const & nodeSets = nodeManager.sets().wrappers();
  for( int i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    createSet( setName );
  }

  // Then loop over them in parallel.
  forAll< parallelHostPolicy >( nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    SortedArrayView< localIndex const > const targetSet = nodeManager.sets().getReference< SortedArray< localIndex > >( setName ).toViewConst();
    constructSetFromSetAndMap( targetSet, m_toNodesRelation, setName );
  } );

  setDomainBoundaryObjects( faceManager );
}

void EdgeManager::buildEdges( localIndex const numNodes,
                              ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                              ArrayOfArrays< localIndex > & faceToEdgeMap )
{
  ArrayOfArrays< EdgeBuilder > edgesByLowestNode( numNodes, 2 * maxEdgesPerNode() );
  createEdgesByLowestNode( faceToNodeMap, edgesByLowestNode.toView() );

  array1d< localIndex > uniqueEdgeOffsets( numNodes + 1 );
  localIndex const numEdges = calculateTotalNumberOfEdges( edgesByLowestNode.toViewConst(), uniqueEdgeOffsets );

  resizeEdgeToFaceMap( edgesByLowestNode.toViewConst(),
                       uniqueEdgeOffsets,
                       m_toFacesRelation );

  resize( numEdges );

  populateMaps( edgesByLowestNode.toViewConst(),
                uniqueEdgeOffsets,
                faceToNodeMap,
                faceToEdgeMap,
                m_toFacesRelation,
                m_toNodesRelation );
}


void EdgeManager::setDomainBoundaryObjects( FaceManager const & faceManager )
{
  // get the "isDomainBoundary" field from the faceManager. This should have
  // been set already!
  arrayView1d< integer const > const & isFaceOnDomainBoundary = faceManager.getDomainBoundaryIndicator();

  // get the "isDomainBoundary" field from for *this, and set it to zero
  arrayView1d< integer > const & isEdgeOnDomainBoundary = this->getDomainBoundaryIndicator();
  isEdgeOnDomainBoundary.zero();

  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

  // loop through all faces
  for( localIndex kf=0; kf<faceManager.size(); ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isFaceOnDomainBoundary[kf] == 1 )
    {
      localIndex const numFaceEdges = faceToEdgeMap.sizeOfArray( kf );

      // loop over all nodes connected to face, and set isNodeDomainBoundary
      for( localIndex a = 0; a < numFaceEdges; ++a )
      {
        isEdgeOnDomainBoundary( faceToEdgeMap( kf, a ) ) = 1;
      }
    }
  }
}

bool EdgeManager::hasNode( const localIndex edgeID, const localIndex nodeID ) const
{
  return m_toNodesRelation( edgeID, 0 ) == nodeID || m_toNodesRelation( edgeID, 1 ) == nodeID;
}

//localIndex EdgeManager::FindEdgeFromNodeIDs(const localIndex nodeA, const localIndex nodeB, const NodeManager *
// nodeManager)
//{
//  localIndex val = std::numeric_limits<localIndex>::max();
//
//  if (nodeA == nodeB)
//    return (val);
//
//  for( SortedArray<localIndex>::const_iterator iedge=nodeManager->m_nodeToEdgeMap[nodeA].begin() ;
//       iedge!=nodeManager->m_nodeToEdgeMap[nodeA].end() ; ++iedge )
//  {
//    if (hasNode(*iedge, nodeB))
//      val = *iedge;
//  }
//  return(val);
//}

void EdgeManager::setIsExternal( FaceManager const & faceManager )
{
  // get the "isExternal" field from the faceManager...This should have been
  // set already!
  arrayView1d< integer const > const & isExternalFace = faceManager.isExternal();

  ArrayOfArraysView< localIndex const > const & faceToEdges = faceManager.edgeList().toViewConst();

  // get the "isExternal" field from for *this, and set it to zero
  m_isExternal.zero();

  // loop through all faces
  for( localIndex kf=0; kf<faceManager.size(); ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isExternalFace[kf] == 1 )
    {
      // loop over all nodes connected to face, and set isNodeDomainBoundary
      localIndex const numEdges = faceToEdges.sizeOfArray( kf );
      for( localIndex a = 0; a < numEdges; ++a )
      {
        m_isExternal[ faceToEdges( kf, a ) ] = 1;
      }
    }
  }
}


void EdgeManager::extractMapFromObjectForAssignGlobalIndexNumbers( NodeManager const & nodeManager,
                                                                   std::vector< std::vector< globalIndex > > & globalEdgeNodes )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numEdges = size();

  arrayView2d< localIndex const > const edgeNodes = this->nodeList();
  arrayView1d< integer const > const isDomainBoundary = this->getDomainBoundaryIndicator();

  globalEdgeNodes.resize( numEdges );

  forAll< parallelHostPolicy >( numEdges, [&]( localIndex const edgeID )
  {
    std::vector< globalIndex > & curEdgeGlobalNodes = globalEdgeNodes[ edgeID ];

    if( isDomainBoundary( edgeID ) )
    {
      curEdgeGlobalNodes.resize( 2 );

      for( localIndex a = 0; a < 2; ++a )
      {
        curEdgeGlobalNodes[ a ]= nodeManager.localToGlobalMap()( edgeNodes[ edgeID ][ a ] );
      }

      std::sort( curEdgeGlobalNodes.begin(), curEdgeGlobalNodes.end() );
    }
  } );
}


void EdgeManager::connectivityFromGlobalToLocal( const SortedArray< localIndex > & indices,
                                                 const map< globalIndex, localIndex > & nodeGlobalToLocal,
                                                 const map< globalIndex, localIndex > & GEOSX_UNUSED_PARAM( faceGlobalToLocal ) )
{


  for( localIndex const ke : indices )
  {
    for( localIndex a=0; a<m_toNodesRelation.size( 1 ); ++a )
    {
      const globalIndex gnode = m_toNodesRelation( ke, a );
      const localIndex lnode = stlMapLookup( nodeGlobalToLocal, gnode );
      m_toNodesRelation( ke, a ) = lnode;
    }
  }

//  array1d<SortedArray<localIndex>>* const edgesToFlowFaces =
// GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
//  if( edgesToFlowFaces != NULL )
//  {
//    for( SortedArray<localIndex>::const_iterator ke=indices.begin() ; ke!=indices.end() ; ++ke )
//    {
//      SortedArray<localIndex>& edgeToFlowFaces = (*edgesToFlowFaces)[*ke];
//      SortedArray<localIndex> newSet;
//      for( SortedArray<localIndex>::iterator faceIndex=edgeToFlowFaces.begin() ;
// faceIndex!=edgeToFlowFaces.end() ; ++faceIndex )
//      {
//
//        std::map<globalIndex,localIndex>::const_iterator MapIter =
// faceGlobalToLocal.find( static_cast<globalIndex>(*faceIndex) );
//        if( MapIter!=faceGlobalToLocal.end()  )
//        {
////          const localIndex faceLocalIndex = stlMapLookup( faceGlobalToLocal,
// static_cast<globalIndex>(*faceIndex) );
//          newSet.insert( MapIter->second );
//
//        }
//      }
//      edgeToFlowFaces = newSet;
//    }
//
//  }

}

localIndex EdgeManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}

localIndex EdgeManager::packUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex EdgeManager::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = localToGlobalMap();
  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();

  localIndex packedSize = bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::nodeListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_toNodesRelation.base().toViewConst(),
                                           m_unmappedGlobalIndicesInToNodes,
                                           packList,
                                           localToGlobal,
                                           nodeLocalToGlobal );


  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_toFacesRelation.base().toArrayOfArraysView(),
                                           m_unmappedGlobalIndicesInToFaces,
                                           packList,
                                           localToGlobal,
                                           faceLocalToGlobal );

  return packedSize;
}



localIndex EdgeManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  // GEOSX_MARK_FUNCTION;

  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOSX_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_toNodesRelation.relatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOSX_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void EdgeManager::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation.base(),
                                    m_toFacesRelation.relatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );
}

void EdgeManager::compressRelationMaps()
{
  m_toFacesRelation.compress();
}

void EdgeManager::depopulateUpMaps( std::set< localIndex > const & receivedEdges,
                                    ArrayOfArraysView< localIndex const > const & facesToEdges )
{
  ObjectManagerBase::cleanUpMap( receivedEdges, m_toFacesRelation.toView(), facesToEdges );
}


} /// namespace geosx
