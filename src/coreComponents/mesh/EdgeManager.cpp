/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @author settgast1
 * @date Jun 22, 2011
 */

#include "EdgeManager.hpp"
#include "NodeManager.hpp"
#include "FaceManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{

EdgeManager::EdgeManager( std::string const & name,
                          ManagedGroup * const parent ):
  ObjectManagerBase(name,parent)
{
  this->RegisterViewWrapper(viewKeyStruct::nodeListString, &this->m_toNodesRelation, 0 );
  this->RegisterViewWrapper(viewKeyStruct::faceListString, &this->m_toFacesRelation, 0 );

  m_toNodesRelation.resize( 0, 2 );
  // TODO Auto-generated constructor stub

}

EdgeManager::~EdgeManager()
{
  // TODO Auto-generated destructor stub
}

void EdgeManager::BuildEdges( FaceManager * const faceManager, NodeManager * const nodeManager )
{
  GEOSX_MARK_FUNCTION;
  if (faceManager->size() == 0 || nodeManager->size() == 0)
    return;

  localIndex numEdges = 0;
  OrderedVariableOneToManyRelation& faceToEdgeMap = faceManager->edgeList();
  faceToEdgeMap.SetRelatedObject( this );
  array1d<localIndex_array>& faceToNodeMap = faceManager->nodeList();

  m_toNodesRelation.SetRelatedObject( nodeManager );
  m_toFacesRelation.SetRelatedObject( faceManager );

  m_toNodesRelation.resize(8 * nodeManager->size());
  m_toFacesRelation.resize(8 * nodeManager->size());

  UnorderedVariableOneToManyRelation & nodeToEdgeMap = nodeManager->edgeList();
  // this will be used to hold a list pairs that store the 2nd node in the edge,
  // and the edge number.
  // they are keyed the edges lowest node.
  array1d<localIndex_array> edgesByLowestNode(nodeManager->size());

  // loop over all the faces
  for( localIndex kf=0 ; kf<faceManager->size() ; ++kf )
  {
    const localIndex numNodesInFace = faceToNodeMap[kf].size();
    const localIndex_array& nodeList = faceToNodeMap[kf];

    localIndex node0, node1, temp;

    // loop over all the nodes in the face. there will be an edge for each
    // node.
    for( localIndex a=0 ; a<numNodesInFace ; ++a )
    {
      // sort the nodes in order of index value
      node0 = nodeList[a];
      node1 = nodeList[ ( a + 1 ) % numNodesInFace ];

      if( node0 > node1 )
      {
        temp = node0;
        node0 = node1;
        node1 = temp;
      }

      // check to see if the edge is in the edgesByLowestNode array (i.e. it
      // is already registered)

      bool duplicate = false;
      localIndex_array& edgesWithSameFirstNode = edgesByLowestNode[node0];
      const localIndex numMatches = edgesWithSameFirstNode.size();
      for( localIndex i = 0; i < numMatches; ++i )
      {
        const localIndex existingEdgeIndex = edgesWithSameFirstNode[i];
        const localIndex existingNode1 = m_toNodesRelation(existingEdgeIndex, 1);
        if( existingNode1 == node1 )
        {
          duplicate = true;
          faceToEdgeMap(kf).push_back(existingEdgeIndex);
          m_toFacesRelation(existingEdgeIndex).insert(kf);
          break;
        }
      }

      // if the edge is not duplicate, then we will assign a new pair into the
      // edgesByLowestNode array
      if( !duplicate )
      {
        edgesByLowestNode[node0].push_back( numEdges );
        faceToEdgeMap(kf).push_back(numEdges);
        m_toFacesRelation(numEdges).insert(kf);
        
        m_toNodesRelation(numEdges, 0) = node0;
        m_toNodesRelation(numEdges, 1) = node1;

        nodeToEdgeMap(node0).insert(numEdges);
        nodeToEdgeMap(node1).insert(numEdges);
        
        numEdges += 1;
      }
    }
  }

  resize(numEdges);

  // make sets from nodesets

  auto const & nodeSets = nodeManager->GetGroup(dataRepository::keys::sets)->wrappers();
  for ( int i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    CreateSet( setName );
  }

  // Then loop over them in parallel.
  raja::forall_in_range<parallelHostPolicy>( 0, nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    const set<localIndex>& targetSet = nodeManager->GetGroup(dataRepository::keys::sets)->getReference<set<localIndex>>( setName );
    ConstructSetFromSetAndMap( targetSet, m_toNodesRelation, setName );
  } );


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
                                                                   array1d<globalIndex_array>& edgesToNodes )
{

  FixedOneToManyRelation const & edgeNodes = this->nodeList();
  integer_array const & isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);

  nodeManager->CheckTypeID( typeid( NodeManager ) );


  edgesToNodes.clear();
  edgesToNodes.resize(size());
  for( localIndex kf=0 ; kf<size() ; ++kf )
  {

    if( isDomainBoundary(kf) != 0 )
    {
      globalIndex_array temp;

      for( localIndex a=0 ; a<edgeNodes.size(1) ; ++a )
      {
        const globalIndex gnode = nodeManager->m_localToGlobalMap( edgeNodes[kf][a] );
        temp.push_back( gnode );
      }
      std::sort( temp.begin(), temp.end() );
      edgesToNodes[kf] = temp;
    }
  }
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
                                                 const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                                 const std::map<globalIndex,localIndex>& faceGlobalToLocal )
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
                                         packList,
                                         m_localToGlobalMap,
                                         m_toNodesRelation.RelatedObjectLocalToGlobal() );


  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::faceListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_toFacesRelation.Base(),
                                         packList,
                                         m_localToGlobalMap,
                                         m_toFacesRelation.RelatedObjectLocalToGlobal() );


  return packedSize;
}



localIndex EdgeManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF( nodeListString != viewKeyStruct::nodeListString, "");

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     this->m_globalToLocalMap,
                                     m_toNodesRelation.RelatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF( faceListString != viewKeyStruct::faceListString, "");

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     this->m_globalToLocalMap,
                                     m_toFacesRelation.RelatedObjectGlobalToLocal(),
                                     false );

  return unPackedSize;
}

}
