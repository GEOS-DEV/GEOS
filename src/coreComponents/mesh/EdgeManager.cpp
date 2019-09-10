/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;

EdgeManager::EdgeManager( std::string const & name,
                          Group * const parent ):
  ObjectManagerBase(name,parent),
  m_edgesToFractureConnectorsEdges(),
  m_fractureConnectorsEdgesToEdges(),
  m_fractureConnectorEdgesToFaceElements()
{
  this->registerWrapper(viewKeyStruct::nodeListString, &this->m_toNodesRelation, 0 );
  this->registerWrapper(viewKeyStruct::faceListString, &this->m_toFacesRelation, 0 );

  m_toNodesRelation.resize( 0, 2 );
  // TODO Auto-generated constructor stub


  registerWrapper( viewKeyStruct::edgesTofractureConnectorsEdgesString, &m_edgesToFractureConnectorsEdges, 0 )->
    setPlotLevel(PlotLevel::NOPLOT)->
    setDescription( "A map of edge local indices to the fracture connector local indices.")->
    setSizedFromParent(0);

  registerWrapper( viewKeyStruct::fractureConnectorEdgesToEdgesString, &m_fractureConnectorsEdgesToEdges, 0 )->
    setPlotLevel(PlotLevel::NOPLOT)->
    setDescription( "A map of fracture connector local indices to edge local indices.")->
    setSizedFromParent(0);

  registerWrapper( viewKeyStruct::fractureConnectorsEdgesToFaceElementsIndexString,
                       &m_fractureConnectorEdgesToFaceElements, 0 )->
    setPlotLevel(PlotLevel::NOPLOT)->
    setDescription( "A map of fracture connector local indices face element local indices")->
    setSizedFromParent(0);

}

EdgeManager::~EdgeManager()
{
  // TODO Auto-generated destructor stub
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
               localIndex const faceLocalEdgeIndex_ ) :
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
    if ( n1 < rhs.n1 ) return true;
    if ( n1 > rhs.n1 ) return false;
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
void createEdgesByLowestNode( arrayView1d< arrayView1d< localIndex const > const > const & faceToNodeMap,
                              ArrayOfArraysView< EdgeBuilder > const & edgesByLowestNode )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numFaces = faceToNodeMap.size();
  
  // loop over all the faces.
  forall_in_range< parallelHostPolicy >( 0, numFaces, [&]( localIndex const faceID ) 
  {
    localIndex const numNodesInFace = faceToNodeMap[ faceID ].size();

    // loop over all the nodes in the face. there will be an edge for each node.
    for( localIndex a=0 ; a< numNodesInFace ; ++a )
    {
      // sort the nodes in order of index value.
      localIndex node0 = faceToNodeMap[ faceID ][ a ];
      localIndex node1 = faceToNodeMap[ faceID ][ ( a + 1 ) % numNodesInFace ];
      if( node0 > node1 ) std::swap( node0, node1 );

      // And append the edge to edgesByLowestNode.
      edgesByLowestNode.atomicAppendToArray( RAJA::atomic::auto_atomic{}, node0, EdgeBuilder( node1, faceID, a ) );
    }
  } );

  // Loop over all the nodes and sort the associated edges.
  forall_in_range< parallelHostPolicy >( 0, numNodes, [&]( localIndex const nodeID )
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
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  GEOS_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );

  uniqueEdgeOffsets[0] = 0;

  // Loop over all the nodes.
  forall_in_range< parallelHostPolicy >( 0, numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // If there are no edges associated with this node we can skip it.
    if ( numEdges == 0 ) return;

    localIndex & numUniqueEdges = uniqueEdgeOffsets[ nodeID + 1 ];
    numUniqueEdges = 0;

    // Otherwise since edgesByLowestNode[ nodeID ] is sorted we can compare subsequent entries
    // count up the unique entries.
    localIndex j = 0;
    for ( ; j < numEdges - 1; ++j )
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
 * @brief Add an edge to the face to edge map, edge to face map, and edge to node map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in/out] faceToEdgeMap the map from face IDs to edge IDs.
 * @param [in/out] edgeToFacemap the map from edgeIDs to faceIDs.
 * @param [in/out] edgeToNodeMap the map from edgeIDs to nodeIDs.
 * @param [in] edgeID the ID of the edge to add.
 * @param [in] firstNodeID the ID of the first node of the edge.
 * @param [in] firstMatch the index of the first EdgeBuilder that describes this edge in edgesByLowestNode[ firstNodeID ].
 * @param [in] numMatches the number of EdgeBuilders that describe this edge in edgesByLowestNode[ firstNodeID ].
 */
void addEdge( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
              arrayView1d< array1d< localIndex > > const & faceToEdgeMap,
              arrayView1d< SortedArray< localIndex > > const & edgeToFaceMap,
              arrayView2d< localIndex > const & edgeToNodeMap,
              localIndex const edgeID,
              localIndex const firstNodeID,
              localIndex const firstMatch,
              localIndex const numMatches )
{
  // Populate the edge to node map.
  edgeToNodeMap( edgeID, 0 ) = firstNodeID;
  edgeToNodeMap( edgeID, 1 ) = edgesByLowestNode( firstNodeID, firstMatch ).n1;

  // Loop through all the matches and fill in the face to edge and edge to face maps.
  for ( localIndex i = 0; i < numMatches; ++i )
  {
    localIndex const faceID = edgesByLowestNode( firstNodeID, firstMatch + i ).faceID;
    localIndex const faceLocalEdgeIndex = edgesByLowestNode( firstNodeID, firstMatch + i ).faceLocalEdgeIndex;

    faceToEdgeMap[ faceID ][ faceLocalEdgeIndex ] = edgeID;
    edgeToFaceMap[ edgeID ].insert( faceID );
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
                   arrayView1d< arrayView1d< localIndex const > const > const & faceToNodeMap, 
                   arrayView1d< array1d< localIndex > > const & faceToEdgeMap,
                   arrayView1d< SortedArray< localIndex > > const & edgeToFaceMap,
                   arrayView2d< localIndex > const & edgeToNodeMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numFaces = faceToNodeMap.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  GEOS_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );
  GEOS_ERROR_IF_NE( numFaces, faceToEdgeMap.size() );
  GEOS_ERROR_IF_NE( numUniqueEdges, edgeToFaceMap.size() );
  GEOS_ERROR_IF_NE( numUniqueEdges, edgeToNodeMap.size( 0 ) );

  // The face to edge map has the same shape as the face to node map, so we can resize appropriately.
  GEOSX_MARK_BEGIN("Reserving space in faceToEdgeMap");
  for ( localIndex faceID = 0; faceID < numFaces; ++faceID )
  {
    faceToEdgeMap[ faceID ].resize( faceToNodeMap[ faceID ].size() );
  }
  GEOSX_MARK_END("Reserving space in faceToEdgeMap");

  // Need to be smarter about this. Should precalculate the number of faces associated with each edge
  // and store it in an array which would be used here.
  GEOSX_MARK_BEGIN("Reserving space in edgeToFaceMap");
  for ( localIndex edgeID = 0; edgeID < numUniqueEdges; ++edgeID )
  {
    edgeToFaceMap[ edgeID ].reserve( 10 );
  }
  GEOSX_MARK_END("Reserving space in edgeToFaceMap");

  // loop over all the nodes.
  forall_in_range< parallelHostPolicy >( 0, numNodes, [&]( localIndex const nodeID )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeID ];
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );
    
    // loop over all the EdgeBuilders associated with the node
    localIndex j = 0;
    while (j < numEdges - 1)
    {
      // Find the number of EdgeBuilders that describe the same edge
      localIndex numMatches = 1;
      while ( edgesByLowestNode( nodeID, j ) == edgesByLowestNode( nodeID, j + numMatches ) )
      {
        ++numMatches;
        if ( j + numMatches == numEdges ) break;
      }

      // Then add the edge.
      addEdge( edgesByLowestNode, faceToEdgeMap, edgeToFaceMap, edgeToNodeMap, curEdgeID, nodeID, j, numMatches );
      ++curEdgeID;
      j += numMatches;
    }

    if ( j == numFaces - 1 )
    {
      addEdge( edgesByLowestNode, faceToEdgeMap, edgeToFaceMap, edgeToNodeMap, curEdgeID, nodeID, j, 1 );
    }
  } );
}

void EdgeManager::BuildEdges( FaceManager * const faceManager, NodeManager * const nodeManager )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = nodeManager->size();
  arrayView1d< arrayView1d< localIndex const > const > const & faceToNodeMap = faceManager->nodeList().toViewConst();

  ArrayOfArrays<EdgeBuilder> edgesByLowestNode( numNodes, 2 * maxEdgesPerNode() );
  createEdgesByLowestNode( faceToNodeMap, edgesByLowestNode );
  
  array1d< localIndex > uniqueEdgeOffsets( numNodes + 1 );
  localIndex const numEdges = calculateTotalNumberOfEdges( edgesByLowestNode, uniqueEdgeOffsets );

  OrderedVariableOneToManyRelation& faceToEdgeMap = faceManager->edgeList();
  faceToEdgeMap.SetRelatedObject( this );

  m_toNodesRelation.SetRelatedObject( nodeManager );
  m_toFacesRelation.SetRelatedObject( faceManager );

  GEOSX_MARK_BEGIN("EdgeManager::resize");
  resize(numEdges);
  GEOSX_MARK_END("EdgeManager::resize");

  populateMaps( edgesByLowestNode,
                uniqueEdgeOffsets,
                faceToNodeMap, 
                faceToEdgeMap,
                m_toFacesRelation,
                m_toNodesRelation );

  // make sets from nodesets
  auto const & nodeSets = nodeManager->sets()->wrappers();
  for ( int i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    CreateSet( setName );
  }

  // Then loop over them in parallel.
  GEOSX_MARK_BEGIN("Set construction");
  forall_in_range<parallelHostPolicy>( 0, nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    const set<localIndex>& targetSet = nodeManager->sets()->getReference<set<localIndex>>( setName );
    ConstructSetFromSetAndMap( targetSet, m_toNodesRelation, setName );
  } );
  GEOSX_MARK_END("Set construction");

  SetDomainBoundaryObjects( faceManager );
}


/// Calculates the midpoint of the edge
//void EdgeManager::EdgeCenter(const NodeManager * nodeManager, localIndex edge,
// R1Tensor& center) const
//{
//
//  const array1d< R1Tensor >& refPosition =
// nodeManager->GetFieldData<FieldInfo::referencePosition>();
//  const array1d< R1Tensor >& displacement =
// nodeManager->GetFieldData<FieldInfo::displacement>();
//  FixedOneToManyRelation const & m_toNodesRelation =
// this->getWrapper<FixedOneToManyRelation>(string("edgesToNodes")).reference();
////  UnorderedVariableOneToManyRelation const & toFacesRelation =
// this->getWrapper<UnorderedVariableOneToManyRelation>(string("edgesToFaces")).reference();
//
//  if (m_toNodesRelation.Dimension(1) >= 2)
//  {
//    const localIndex& node0 = m_toNodesRelation(edge,0);
//    center =  refPosition[node0];
//    center += displacement[node0];
//    const localIndex& node1 = m_toNodesRelation(edge,1);
//    center += refPosition[node1];
//    center += displacement[node1];
//    center *= 0.5;
//  }
//  else
//  {
//    const localIndex& node0 = m_toNodesRelation(edge,0);
//    center =  refPosition[node0];
//  }
//}

//
///// Calculates the vector from node 0 to node 1
//void EdgeManager::EdgeVector(const NodeManager * nodeManager, localIndex edge,
// R1Tensor& v) const{
//  const array1d< R1Tensor >& refPosition =
// nodeManager->GetFieldData<FieldInfo::referencePosition>();
//  const array1d< R1Tensor >& displacement =
// nodeManager->GetFieldData<FieldInfo::displacement>();
//  FixedOneToManyRelation const & m_toNodesRelation =
// this->getWrapper<FixedOneToManyRelation>(string("edgesToNodes")).reference();
//
//  const localIndex& node1 = m_toNodesRelation(edge,1);
//  v =  refPosition[node1];
//  v += displacement[node1];
//  const localIndex& node0 = m_toNodesRelation(edge,0);
//  v -= refPosition[node0];
//  v -= displacement[node0];
//}
//
///// Returns the length of the edge
//realT EdgeManager::EdgeLength(const NodeManager * nodeManager, localIndex
// edge) const{
//  const array1d< R1Tensor >& refPosition =
// nodeManager->GetFieldData<FieldInfo::referencePosition>();
//  const array1d< R1Tensor >& displacement =
// nodeManager->GetFieldData<FieldInfo::displacement>();
//  FixedOneToManyRelation const & m_toNodesRelation =
// this->getWrapper<FixedOneToManyRelation>(string("edgesToNodes")).reference();
//
//  if (m_toNodesRelation.Dimension(1) >= 2)
//  {
//  const localIndex& node0 = m_toNodesRelation(edge,0);
//  const localIndex& node1 = m_toNodesRelation(edge,1);
//  R1Tensor v =  refPosition[node0];
//  v += displacement[node0];
//  v -= refPosition[node1];
//  v -= displacement[node1];
//  return v.L2_Norm();
//  }
//  else
//  {
//    return 1.0;
//  }
//}


void EdgeManager::SetDomainBoundaryObjects( const ObjectDataStructureBaseT * const referenceObject )
{
  GEOSX_MARK_FUNCTION;

  // make sure that the reference object is a faceManger object
//  referenceObject->CheckObjectType( ObjectDataStructureBaseT::FaceManager );

  // cast the referenceObject into a faceManager
  const FaceManager * faceManager = static_cast<const FaceManager *>( referenceObject);

  // get the "isDomainBoundary" field from the faceManager->..This should have
  // been set already!
  const array1d<integer>& isFaceOnDomainBoundary = faceManager->getReference< array1d<integer> >( ObjectManagerBase::viewKeyStruct::domainBoundaryIndicatorString );

  // get the "isDomainBoudnary" field from for *this, and set it to zero
  array1d<integer>& isNodeOnDomainBoundary = this->getReference< array1d<integer> >(ObjectManagerBase::viewKeyStruct::domainBoundaryIndicatorString);
  isNodeOnDomainBoundary = 0;

  // loop through all faces
  for( localIndex kf=0 ; kf<faceManager->size() ; ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isFaceOnDomainBoundary[kf] == 1 )
    {
      localIndex_array const & faceToEdges = faceManager->edgeList()(kf);

      // loop over all nodes connected to face, and set isNodeDomainBoundary
      for( localIndex_array::const_iterator a=faceToEdges.begin() ; a!=faceToEdges.end() ; ++a )
      {
        isNodeOnDomainBoundary(*a) = 1;
      }
    }
  }
}

bool EdgeManager::hasNode( const localIndex edgeID, const localIndex nodeID ) const
{

  if( m_toNodesRelation(edgeID,0) == nodeID || m_toNodesRelation(edgeID,1) == nodeID )
  {
    return true;
  }
  else
    return false;
}

//localIndex EdgeManager::FindEdgeFromNodeIDs(const localIndex nodeA, const localIndex nodeB, const NodeManager * nodeManager)
//{
//  localIndex val = std::numeric_limits<localIndex>::max();
//
//  if (nodeA == nodeB)
//    return (val);
//
//  for( set<localIndex>::const_iterator iedge=nodeManager->m_nodeToEdgeMap[nodeA].begin() ;
//       iedge!=nodeManager->m_nodeToEdgeMap[nodeA].end() ; ++iedge )
//  {
//    if (hasNode(*iedge, nodeB))
//      val = *iedge;
//  }
//  return(val);
//}

void EdgeManager::SetIsExternal( FaceManager const * const faceManager )
{
  // get the "isExternal" field from the faceManager->..This should have been
  // set already!
  array1d<integer> const & isExternalFace = faceManager->isExternal();
  array1d<integer>& isExternal = this->isExternal();

  array1d<localIndex_array> const & facesToEdges = faceManager->edgeList();

  // get the "isExternal" field from for *this, and set it to zero
  isExternal = 0;


  // loop through all faces
  for( localIndex kf=0 ; kf<faceManager->size() ; ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isExternalFace[kf] == 1 )
    {
      localIndex_array const & faceToEdges = facesToEdges[kf];

      // loop over all nodes connected to face, and set isNodeDomainBoundary

      for( auto a : faceToEdges )
      {
        isExternal[a] = 1;
      }
    }
  }
}


void EdgeManager::ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const nodeManager,
                                                                   std::vector< std::vector< globalIndex > > & globalEdgeNodes )
{
  GEOSX_MARK_FUNCTION;
  nodeManager->CheckTypeID( typeid( NodeManager ) );

  localIndex const numEdges = size();

  arrayView2d< localIndex const > const & edgeNodes = this->nodeList();
  arrayView1d< integer const > const & isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);

  globalEdgeNodes.resize( numEdges );

  forall_in_range< parallelHostPolicy >( 0, numEdges, [&]( localIndex const edgeID ) 
  {
    std::vector< globalIndex > & curEdgeGlobalNodes = globalEdgeNodes[ edgeID ];

    if( isDomainBoundary( edgeID ) )
    {
      curEdgeGlobalNodes.resize( 2 );

      for ( localIndex a = 0; a < 2 ; ++a )
      {
        curEdgeGlobalNodes[ a ]= nodeManager->m_localToGlobalMap( edgeNodes[ edgeID ][ a ] );
      }

      std::sort( curEdgeGlobalNodes.begin(), curEdgeGlobalNodes.end() );
    }
  } );
}

#if 0
// It seems like this function is not used anywhere
void EdgeManager::SplitEdge( const localIndex indexToSplit,
                             const localIndex parentNodeIndex,
                             const localIndex childNodeIndex[2],
                             array1d<set<localIndex>>& nodesToEdges )
{

  localIndex newEdgeIndex[2];

  bool didSplit = SplitObject( indexToSplit, newEdgeIndex );

  // copy the parent edges edgeToNodes relation, replacing the parentNode with
  // one of the new nodes.

  // loop over each new edge
  for( int ke=0 ; ke<2 ; ++ke )
  {


    // loop over each node on the edge
    for( int a=0 ; a<2 ; ++a )
    {
      // nodeIndex is the node on the parent edge
      const localIndex& nodeIndex = m_m_toNodesRelation(indexToSplit,a);

      // modify the edgesToNodes for new edges. They should point at a
      // combination of one of the parents nodes,
      // one or two of the new nodes that were split.

      // if the nodeIndex==parentNodeIndex then the parent node is the target
      // node
      if( nodeIndex == parentNodeIndex )
      {
        // adding the child node to the edgesToNodes map of the child edge
        m_m_toNodesRelation(newEdgeIndex[ke],a) = childNodeIndex[ke];
        std::cout<<"    m_edgesToNodes("<<newEdgeIndex[ke]<<","<<a<<") = "<<childNodeIndex[ke]<<std::endl;

        nodesToEdges(childNodeIndex[ke]).insert(newEdgeIndex[ke]);
//        std::cout<<"
//
// nodesToEdges("<<childNodeIndex[ke]<<").insert("<<newEdgeIndex[ke]<<")"<<std::endl;
      }
      else
      {
        if( didSplit )
        {
          m_m_toNodesRelation(newEdgeIndex[ke],a) = nodeIndex;
          std::cout<<"    m_edgesToNodes("<<newEdgeIndex[ke]<<","<<a<<") = "<<nodeIndex<<std::endl;

          nodesToEdges(nodeIndex).insert(newEdgeIndex[ke]);
//          std::cout<<"
//    nodesToEdges("<<nodeIndex<<").insert("<<newEdgeIndex[ke]<<")"<<std::endl;
        }
      }
    }
  }



}


template< typename T_indices >
unsigned int EdgeManager::PackEdges( const T_indices& sendedges,
                                     const NodeManager& nodeManager,
                                     const FaceManager& faceManager,
                                     bufvector& buffer,
                                     const bool packConnectivityToGlobal,
                                     const bool packFields,
                                     const bool packMaps,
                                     const bool packSets  ) const
{

  unsigned int sizeOfPacked = 0;

  const std::string label = "EdgeData";
  sizeOfPacked += buffer.Pack(label);

  // pack data in the base
  sizeOfPacked += ObjectDataStructureBaseT::PackBaseObjectData( buffer, sendedges, packFields, packMaps, packSets, packConnectivityToGlobal );

  // pack the edge specific data
  for( typename T_indices::const_iterator edgeIndex=sendedges.begin() ; edgeIndex!=sendedges.end() ; ++edgeIndex )
  {
    const localIndex* const nodelist = m_m_toNodesRelation[*edgeIndex];

    int numnodes =m_m_toNodesRelation.Dimension(1);
    sizeOfPacked += buffer.Pack(numnodes);

    for( int a=0 ; a<numnodes ; ++a )
    {
      globalIndex gnode = GLOBALINDEX_MAX;
      if( packConnectivityToGlobal )
      {
        gnode = nodeManager->m_localToGlobalMap(nodelist[a]);
      }
      else
      {
        gnode = nodelist[a];
      }
      sizeOfPacked += buffer.Pack(gnode);
    }
  }


  const array1d<set<localIndex>>* const edgesToFlowFaces = GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
  if( edgesToFlowFaces != NULL )
  {
    for( typename T_indices::const_iterator edgeIndex=sendedges.begin() ; edgeIndex!=sendedges.end() ; ++edgeIndex )
    {
      const set<localIndex>& edgeToFlowFaces = (*edgesToFlowFaces)[*edgeIndex];
//      std::cout<<"pack edgeIndex, size =
// "<<this->m_localToGlobalMap(*edgeIndex)<<", "<<edgeToFlowFaces.size()<<" ";


      if( packConnectivityToGlobal )
      {
        sizeOfPacked += buffer.PackGlobal( edgeToFlowFaces, faceManager->m_localToGlobalMap );
//        std::cout<<"global"<<std::endl;
      }
      else
      {
        sizeOfPacked += buffer.Pack( edgeToFlowFaces );
//        std::cout<<"local"<<std::endl;
      }
    }
  }


  return sizeOfPacked;
}
template unsigned int EdgeManager::PackEdges( const set<localIndex>&, const NodeManager&, const FaceManager&, bufvector&, const bool, const bool, const bool,
                                              const bool ) const;
template unsigned int EdgeManager::PackEdges( const localIndex_array&, const NodeManager&, const FaceManager&, bufvector&, const bool, const bool, const bool,
                                              const bool ) const;



unsigned int EdgeManager::UnpackEdges( const char*& buffer,
                                       const NodeManager& nodeManager,
                                       const FaceManager& faceManager,
                                       localIndex_array& edgeReceiveLocalIndices,
                                       const bool unpackConnectivityToLocal,
                                       const bool unpackFields,
                                       const bool unpackMaps,
                                       const bool unpackSets  )
{
  unsigned int sizeOfUnpacked = 0;

  const std::string label = "EdgeData";
  std::string temp;

  sizeOfUnpacked += bufvector::Unpack( buffer, temp );
  if( label.compare(temp)!=0 )
  {
    throw GPException("EdgeManager::UnpackEdges: buffer location incorrect\n");
  }

  localIndex_array newEdgeLocalIndices;
  // unpack data from base object
  sizeOfUnpacked += ObjectDataStructureBaseT::UnpackBaseObjectData( buffer, edgeReceiveLocalIndices, newEdgeLocalIndices, unpackFields, unpackMaps, unpackSets,
                                                                    unpackConnectivityToLocal );

  const localIndex_array::size_type numReceivedEdges = edgeReceiveLocalIndices.size();
  // unpack face specific data
  for( localIndex_array::size_type ke=0 ; ke<numReceivedEdges ; ++ke )
  {
    int numnodes;
    sizeOfUnpacked += bufvector::Unpack( buffer, numnodes );
    for( int a=0 ; a<numnodes ; ++a )
    {
      globalIndex gnode;
      sizeOfUnpacked += bufvector::Unpack( buffer, gnode );

      if( unpackConnectivityToLocal )
      {
        const localIndex lnode = stlMapLookup( nodeManager->m_globalToLocalMap, gnode );
        m_m_toNodesRelation(edgeReceiveLocalIndices(ke),a) = lnode;
      }
      else
      {
        m_m_toNodesRelation(edgeReceiveLocalIndices(ke),a) = gnode;
      }
    }

  }


  array1d<set<localIndex>>* const edgesToFlowFaces = GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
  if( edgesToFlowFaces != NULL )
  {
    for( localIndex_array::iterator edgeIndex=edgeReceiveLocalIndices.begin() ; edgeIndex!=edgeReceiveLocalIndices.end() ; ++edgeIndex )
    {
      set<localIndex> newEdgeToFlowFaces;
      set<localIndex>& edgeToFlowFaces = (*edgesToFlowFaces)[*edgeIndex];


      if( unpackConnectivityToLocal )
      {
//        std::cout<<"global unpack edgeIndex, size =
// "<<this->m_localToGlobalMap(*edgeIndex)<<", ";
        sizeOfUnpacked += bufvector::UnpackGlobal( buffer, faceManager->m_globalToLocalMap, newEdgeToFlowFaces );
      }
      else
      {
//        std::cout<<"local unpack edgeIndex, size =
// "<<this->m_localToGlobalMap(*edgeIndex)<<", ";
        sizeOfUnpacked += bufvector::Unpack( buffer, newEdgeToFlowFaces );

        set<globalIndex> globals;
        for( set<localIndex>::iterator i=edgeToFlowFaces.begin() ; i!=edgeToFlowFaces.end() ; ++i )
        {
          globalIndex gi = faceManager->m_localToGlobalMap[*i];
          globals.insert(gi);
        }
        edgeToFlowFaces.clear();
        edgeToFlowFaces.insert( globals.begin(), globals.end() );
      }
//      std::cout<<edgeToFlowFaces.size()<<std::endl;

      edgeToFlowFaces.insert( newEdgeToFlowFaces.begin(), newEdgeToFlowFaces.end() );
    }
  }


  return sizeOfUnpacked;

}
#endif



void EdgeManager::ConnectivityFromGlobalToLocal( const set<localIndex>& indices,
                                                 const map<globalIndex,localIndex>& nodeGlobalToLocal,
                                                 const map<globalIndex,localIndex>& faceGlobalToLocal )
{


  for( set<localIndex>::const_iterator ke=indices.begin() ; ke!=indices.end() ; ++ke )
  {
    for( localIndex a=0 ; a<m_toNodesRelation.size(1) ; ++a )
    {
      const globalIndex gnode = m_toNodesRelation(*ke,a);
      const localIndex lnode = stlMapLookup( nodeGlobalToLocal, gnode );
      m_toNodesRelation(*ke,a) = lnode;
    }
  }

//  array1d<set<localIndex>>* const edgesToFlowFaces =
// GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
//  if( edgesToFlowFaces != NULL )
//  {
//    for( set<localIndex>::const_iterator ke=indices.begin() ; ke!=indices.end() ; ++ke )
//    {
//      set<localIndex>& edgeToFlowFaces = (*edgesToFlowFaces)[*ke];
//      set<localIndex> newSet;
//      for( set<localIndex>::iterator faceIndex=edgeToFlowFaces.begin() ;
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

// Fu note on 20130416: Looks like this was temporary.  This function is now
// taken care of by element region.
// We will keep it here for a while and delete it later.
//void EdgeManager::UpdateEdgeExternalityFromSplit( const FaceManager&
// faceManager,
//                                                 const set<localIndex>& newEdgeIndices,
//                                                 const set<localIndex>&
// modifiedEdgeIndices )
//{
//  set<localIndex> allEdges;
//  allEdges.insert( newEdgeIndices.begin(), newEdgeIndices.end() );
//  allEdges.insert( modifiedEdgeIndices.begin(), modifiedEdgeIndices.end() );
//
//
//  for( set<localIndex>::const_iterator edgeIndex=allEdges.begin() ;
// edgeIndex!=allEdges.end() ; ++edgeIndex )
//  {
//
//    for( set<localIndex>::const_iterator iface=m_toFacesRelation[*edgeIndex].begin() ;
//        iface!=m_toFacesRelation[*edgeIndex].end() ; ++iface )
//    {
//      if (faceManager->m_isExternal[*iface] == 1)
//      {
//        m_isExternal[*edgeIndex] =1;
//      }
//    }
//
//  }
//}

void EdgeManager::AddToEdgeToFaceMap( const FaceManager * faceManager,
                                      const localIndex_array& newFaceIndices )
{
  OrderedVariableOneToManyRelation const & facesToEdges = faceManager->edgeList();
  // loop over all faces in list
  for( localIndex kf=0 ; kf<newFaceIndices.size() ; ++kf )
  {
    localIndex lfi = newFaceIndices(kf);
    // now iterate over the faceToNodeMap (i.e. all nodes in the faceToNodeMap)
    for( localIndex_array::const_iterator a=facesToEdges[lfi].begin() ;
         a != facesToEdges[lfi].end() ; ++a )
    {
      // enter the value of the face index into the nodeToFace map
      m_toFacesRelation[*a].insert(lfi);
    }
  }
}

//void EdgeManager::SetLayersFromDomainBoundary(const NodeManager * const nodeManager)
//{
//
//  array1d<integer>& layersEdge = this->GetFieldData<int>("LayersFromDomainBoundary");
//  const array1d<integer>& layersNode = nodeManager->GetFieldData<int>("LayersFromDomainBoundary");
//
//  for (localIndex ie = 0 ; ie!= this->size() ; ++ie)
//  {
//    for( unsigned int a=0 ; a<m_toNodesRelation.Dimension(1) ; ++a )
//    {
//      if (a==0)
//      {
//        layersEdge[ie] = layersNode[m_toNodesRelation(ie,a)];
//      }
//      else
//      {
//        layersEdge[ie] = std::max(layersNode[m_toNodesRelation(ie,a)], layersEdge[ie]);
//      }
//    }
//  }
//}




localIndex EdgeManager::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

localIndex EdgeManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template<bool DOPACK>
localIndex EdgeManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d<localIndex const> const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::nodeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_toNodesRelation.Base(),
                                         m_unmappedGlobalIndicesInToNodes,
                                         packList,
                                         m_localToGlobalMap,
                                         m_toNodesRelation.RelatedObjectLocalToGlobal() );


  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::faceListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_toFacesRelation.Base(),
                                         m_unmappedGlobalIndicesInToFaces,
                                         packList,
                                         m_localToGlobalMap,
                                         m_toFacesRelation.RelatedObjectLocalToGlobal() );

  return packedSize;
}



localIndex EdgeManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const overwriteDownMaps )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->m_globalToLocalMap,
                                     m_toNodesRelation.RelatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->m_globalToLocalMap,
                                     m_toFacesRelation.RelatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void EdgeManager::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_toFacesRelation,
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );
}


void EdgeManager::depopulateUpMaps( std::set<localIndex> const & receivedEdges,
                                    array1d< array1d< localIndex > > const & facesToEdges )
{
  ObjectManagerBase::CleanUpMap( receivedEdges, m_toFacesRelation, facesToEdges );
}


}
