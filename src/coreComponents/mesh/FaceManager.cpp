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
 * @file FaceManager.cpp
 * @author settgast1
 */

#include "FaceManager.hpp"
#include "ElementRegionManager.hpp"
#include "NodeManager.hpp"
#include "BufferOps.hpp"

namespace geosx
{
using namespace dataRepository;

/**
 *
 * @return
 */
FaceManager::FaceManager( string const &, ManagedGroup * const parent ):
  ObjectManagerBase("FaceManager",parent)
{

  this->RegisterViewWrapper( viewKeyStruct::nodeListString, &m_nodeList, false );
  this->RegisterViewWrapper( viewKeyStruct::edgeListString, &m_edgeList, false );
//  m_nodeList.SetRelatedObject( parent->getGroup<NodeManager>(MeshLevel::groupStructKeys::nodeManagerString));

  this->RegisterViewWrapper( viewKeyStruct::elementRegionListString,
                             &(elementRegionList()),
                             false );

  this->RegisterViewWrapper( viewKeyStruct::elementSubRegionListString,
                             &(elementSubRegionList()),
                             false );

  this->RegisterViewWrapper( viewKeyStruct::elementListString,
                             &(elementList()),
                             false );

  this->RegisterViewWrapper( viewKeyStruct::faceCenterString, &m_faceCenter, false);

  m_toElements.resize(0,2);

  //0-based; note that the following field is ALSO 0
  //for faces that are not external faces, so check isExternal before using
//  this->AddKeylessDataField<localIndex>("externalFaceIndex", true, true);
//
//  this->AddKeylessDataField<R1Tensor>("FaceCenter",true,true);
}

/**
 *
 * @return
 */
FaceManager::~FaceManager()
{}



void FaceManager::FillDocumentationNode()
{
  ObjectManagerBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( "InternalMesh" );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "a mesh generator" );


//  docNode->AllocateChildNode( viewKeys.elementRegionList.Key(),
//                              viewKeys.elementRegionList.Key(),
//                              -1,
//                              "integer_array",
//                              "integer_array",
//                              "List containing the element regions of the faces",
//                              "List containing the element regions of the faces",
//                              "1",
//                              "",
//                              1,
//                              0,
//                              0 );



}


//
void FaceManager::BuildFaces( NodeManager * const nodeManager, ElementRegionManager * const elementManager )
{

  localIndex_array tempNodeList;
  array1d<localIndex_array> tempFaceToNodeMap;

  localIndex numFaces = 0;
  array1d<localIndex_array> facesByLowestNode;

  array2d<localIndex> & elemRegionList = this->elementRegionList();
  array2d<localIndex> & elemSubRegionList = this->elementSubRegionList();
  array2d<localIndex> & elemList = this->elementList();

  m_toElements.setElementRegionManager( elementManager );

  elemRegionList.resize( 4*nodeManager->size() );
  elemSubRegionList.resize( 4*nodeManager->size() );
  elemList.resize( 4*nodeManager->size() );

  elemRegionList = -1;
  elemSubRegionList = -1;
  elemList = -1;

  for( typename dataRepository::indexType kReg=0 ; kReg<elementManager->numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = elementManager->GetRegion(kReg);

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(kSubReg);

      for( localIndex ke=0 ; ke<subRegion->size() ; ++ke )
      {
        // kelf = k'th element local face index
        for( localIndex kelf=0 ; kelf<subRegion->numFacesPerElement() ; ++kelf )
        {
          // get the nodes associated with the local face
          subRegion->GetFaceNodes( ke, kelf, tempNodeList );
          FixedOneToManyRelation & elemsToFaces = subRegion->faceList();

          //Special treatment for the triangle faces of prisms.
          if (tempNodeList[tempNodeList.size()-1] == std::numeric_limits<localIndex>::max())
            tempNodeList.pop_back();

          // sort the nodes
          std::sort(tempNodeList.begin(), tempNodeList.end() );

          // get the lowest node index from the list for simplicity
          const localIndex& lowNode = tempNodeList[0];

          // now check to see if the lowest node index has an entry in the
          // facesByLowestNode vector
          if( facesByLowestNode.size() < (lowNode+1) )
          {
            // the node has not been entered, so add it.
            facesByLowestNode.resize(lowNode+1);

            // this a new face, so add it,
            AddNewFace( kReg, kSubReg, ke, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, *subRegion );
          }
          else
          {
            // does the node have an entry? If not, then this has to be a new
            // face.
            if( facesByLowestNode[lowNode].empty() )
            {
              // this a new face, so add it,
              AddNewFace( kReg, kSubReg, ke, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, *subRegion );
            }
            else
            {
              // the node does have an entry, so it is possible that the facet
              // has already be assigned a number

              // make a flag to indicate whether the face is a
              // duplicate...assume that it isn't unless this is disproved.
              bool duplicate = false;


              // there are faces in facesByLowestNode, so lets loop over them
              // and check for duplicates
              for( localIndex_array::iterator existingFaceIndex = facesByLowestNode[lowNode].begin() ;
                   existingFaceIndex != facesByLowestNode[lowNode].end() ; ++existingFaceIndex )
              {
                // this is the nodelist of the face that we are testing agains
                const localIndex_array& existingFaceNodelist = tempFaceToNodeMap[*existingFaceIndex];

                // test to see if the size of the nodelists are the same....
                if( existingFaceNodelist.size() == tempNodeList.size() )
                {
                  // since the size is the same, then we should test the
                  // nodes...they are sorted, so
                  // the std::equal() algorithm will work for this.
                  if( std::equal( existingFaceNodelist.begin(), existingFaceNodelist.end(), tempNodeList.begin() ) )
                  {
                    // they are equal!
                    duplicate = true;

                    // add the element to the faceToElement map
//                        m_toElementsRelation[*existingFaceIndex].push_back(
// std::pair<ElementRegionT*, localIndex>(
// const_cast<ElementRegionT*>(&elementRegion), ke) );


                    if( elemRegionList[*existingFaceIndex][0] == -1 )
                    {
                      elemRegionList[*existingFaceIndex][0]    = kReg;
                      elemSubRegionList[*existingFaceIndex][0] = kSubReg;
                      elemList[*existingFaceIndex][0]          = ke;
                    }
                    else
                    {
                      elemRegionList[*existingFaceIndex][1]    = kReg;
                      elemSubRegionList[*existingFaceIndex][1] = kSubReg;
                      elemList[*existingFaceIndex][1]          = ke;
                    }


                    // add the face to the elementToFaceMap for the element
                    // region.
                    elemsToFaces(ke,kelf) = *existingFaceIndex;

                    // now remove the entry from the face that we were checking
                    // against from the facesByLowestNode list...
                    // because it is no longer possible that it will have
                    // another element that has this face.
                    facesByLowestNode[lowNode].erase( existingFaceIndex );

                    // break the loop
                    break;
                  }
                }
              }
              if( !duplicate )
              {
                // the face is not a duplicate of any in the facesByLowestNode
                // list, so we need to add a new face.
                AddNewFace( kReg, kSubReg, ke, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, *subRegion );
              }
            }
          }
        }
      }
    }

  }

  // resize the data vectors according to the number of faces indicated by the
  // size of tempFaceToNodeMap
  this->resize(tempFaceToNodeMap.size());

  // set m_FaceToNodeMap
  array1d< localIndex_array > & faceToNodes = m_nodeList;
  faceToNodes = tempFaceToNodeMap;

  auto const & nodeSets = nodeManager->GetGroup(string("Sets"))->wrappers();

  // make sets from nodesets
  for( auto const & setWrapper : nodeSets )
  {
    std::string const & setName = setWrapper.second->getName();
    const set<localIndex>& targetSet = nodeManager->GetGroup(keys::sets)->getReference<set<localIndex>>( setName );
    this->ConstructSetFromSetAndMap( targetSet, faceToNodes, setName );
  }


  // sort the face node lists
  SortAllFaceNodes(*nodeManager,*elementManager);

  SetDomainBoundaryObjects( nodeManager );
//  array1d<R1Tensor>& faceCenter = this->GetFieldData<R1Tensor>( "FaceCenter" );
//  for( localIndex k=0 ; k<DataLengths() ; ++k )
//  {
//    FaceCenter( nodeManager, k , faceCenter[k] );
//
//  }
//
//
//  // Figure out if we need to write arbitrary polygons to silo.  Have to do
// this here because cannot do allreduce in the write silo file.
//  int maxNodePerFace(-100), minNodePerFace(1000),
// writeArbitraryPolygonLocal(0);
//  for (localIndex kf = 0; kf < DataLengths(); ++kf)
//  {
//    maxNodePerFace = std::max(maxNodePerFace,
// int(m_toNodesRelation[kf].size()));
//    minNodePerFace = std::min(minNodePerFace,
// int(m_toNodesRelation[kf].size()));
//  }
//  if (maxNodePerFace != minNodePerFace || maxNodePerFace > 4)
// writeArbitraryPolygonLocal = 1;
//  MPI_Allreduce(&writeArbitraryPolygonLocal, &m_writeArbitraryPolygon, 1,
// MPI_INT, MPI_MAX, MPI_COMM_GEOSX);

}



void FaceManager::AddNewFace( localIndex const & kReg,
                              localIndex const & kSubReg,
                              localIndex const & ke,
                              localIndex const & kelf,
                              localIndex & numFaces,
                              array1d<localIndex_array>& facesByLowestNode,
                              localIndex_array& tempNodeList,
                              array1d<localIndex_array>& tempFaceToNodeMap,
                              CellBlockSubRegion & elementRegion )
{
  // and add the face to facesByLowestNode[]
  facesByLowestNode[tempNodeList[0]].push_back(numFaces);

  // add the face to the elementToFaceMap
  elementRegion.faceList()(ke,kelf) = numFaces;

  // add the nodes to the faceToNodeMap
  tempFaceToNodeMap.push_back(tempNodeList);


  // add to the element information to the faceToElementMap
//  tempFaceToElemEntry.push_back( std::pair<ElementRegionT*, localIndex>(
// const_cast<ElementRegionT*>(&elementRegion), ke) );
//  m_toElementsRelation.push_back(tempFaceToElemEntry);
  auto & toElementRegionList = elementRegionList();
  auto & toElementSubRegionList = elementSubRegionList();
  auto & toElementList = elementList();

  if( toElementRegionList[numFaces][0] == -1 )
  {
    toElementRegionList[numFaces][0]    = kReg;
    toElementSubRegionList[numFaces][0] = kSubReg;
    toElementList[numFaces][0]          = ke;
  }
  else
  {
    toElementRegionList[numFaces][1]    = kReg;
    toElementSubRegionList[numFaces][1] = kSubReg;
    toElementList[numFaces][1]          = ke;
  }

  // now increment numFaces to reflect the number of faces rather than the index
  // of the new face
  ++numFaces;

}


void FaceManager::SetDomainBoundaryObjects( NodeManager * const nodeManager )
{
  // Set value of domainBounaryIndicator to one if it is found to have only one elements that it
  // is connected to.
  integer_array & faceDomainBoundaryIndicator = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);
  faceDomainBoundaryIndicator = 0;

  array2d<localIndex> const & elemRegionList = this->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = this->elementSubRegionList();
  array2d<localIndex> const & elemList = this->elementList();

  for( localIndex kf=0 ; kf<size() ; ++kf )
  {
    if( elemRegionList[kf][1] == -1 )
    {
      faceDomainBoundaryIndicator(kf) = 1;
    }
  }

  integer_array & nodeDomainBoundaryIndicator = nodeManager->getReference<integer_array>(nodeManager->viewKeys.domainBoundaryIndicator);
  nodeDomainBoundaryIndicator = 0;

  OrderedVariableOneToManyRelation const & faceToNodesMap = this->nodeList();

  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( faceDomainBoundaryIndicator[k] == 1 )
    {
      arrayView1d<localIndex const> nodelist = faceToNodesMap[k];
      for( localIndex a=0 ; a<nodelist.size() ; ++a )
      {
        nodeDomainBoundaryIndicator[nodelist[a]] = 1;
      }
    }
  }

}

//
//void FaceManager::SetIsExternal()
//{
//  integer_array const & isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();
//  m_isExternal = 0;
//  for( localIndex k=0 ; k<m_numFaces ; ++k )
//  {
//    if( isDomainBoundary[k]==1 && ( m_matchedBoundaryFaces.find(k)==m_matchedBoundaryFaces.end() ) )
//    {
//      m_isExternal[k] = 1;
//      m_externalFaces.insert(k);
//    }
//  }
//}

//void
//FaceManager::
//SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject )
//{
//  array1d< localIndex_array > const & faceToNodes = this->getReference< array1d< localIndex_array > >( viewKeys.nodeList );
//  globalIndex_array const & nodalGlobalIndex = compositionalObject->m_localToGlobalMap;
//  integer_array const & isDomainBoundary = this->getReference<integer_array>(viewKeys.isDomainBoundary);
//
//  mpiBuffer buffer;
//
//  localIndex numFaces;
//  for( localIndex k=0 ; k<size() ; ++k )
//  {
//    if( isDomainBoundary[k] == 1 )
//    {
//    }
//  }
//  for( localIndex k=0 ; k<size() ; ++k )
//  {
//    if( isDomainBoundary[k] == 1 )
//    {
//      CommBufferOps::Pack( buffer, )
//    }
//  }
//}


void FaceManager::SortAllFaceNodes( NodeManager const & nodeManager,
                                    ElementRegionManager const & elemManager )
{
  array2d<localIndex> & elemRegionList = this->elementRegionList();
  array2d<localIndex> & elemSubRegionList = this->elementSubRegionList();
  array2d<localIndex> & elemList = this->elementList();

  for(localIndex kf =0 ; kf < size() ; ++kf )
  {
    ElementRegion const * const elemRegion     = elemManager.GetRegion( elemRegionList[kf][0] );
    CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion( elemSubRegionList[kf][0] );
    R1Tensor elementCenter = subRegion->GetElementCenter( elemList[kf][0],
                                                          nodeManager,
                                                          true );
    SortFaceNodes( nodeManager, elementCenter, kf );
  }
}

void FaceManager::SortFaceNodes( NodeManager const & nodeManager,
                                 R1Tensor const & elementCenter,
                                 const localIndex faceIndex )
{

  arrayView1d<localIndex> faceNodes = nodeList()[faceIndex];
  const localIndex firstNodeIndex = faceNodes[0];
  const localIndex numFaceNodes = faceNodes.size();

  r1_array const & X = nodeManager.referencePosition();

  // get face center (average vertex location) and store node coordinates
  array1d<R1Tensor> faceCoords(numFaceNodes);
  R1Tensor fc;
  for( localIndex n =0 ; n < numFaceNodes ; ++n)
  {
    localIndex nd = faceNodes[n];
    faceCoords[n] = X[nd];
    fc += faceCoords[n];
  }
  fc /= realT(numFaceNodes);

  // find center of element 0
//  R1Tensor ec;
//  elemManager.GetRegion( elemRegionList[faceIndex][0] );
//  ElementIdPair eid = m_toElementsRelation[faceIndex][0];
//  ec = eid.first->GetElementCenter(eid.second,nodeManager);


  R1Tensor ex, ey, ez;
  // Approximate face normal direction (unscaled)
  if (numFaceNodes == 2)  //2D only.
  {
    ex = X[faceNodes[1]];
    ex -= X[faceNodes[0]];
    ey = elementCenter;
    ey -= fc;

  }
//  else if (eid.first->m_numFacesPerElement == 1)
//  {
//    //  The original/default algorithm does not work for shell elements where
// the face is the element itself
//    //  In the new algorithm, we construct ez based on two vectors ex and ey
// in the plane.
//    //Because ex and ey are generally not perpendicular to each other, we have
// to replace ey after we get ez.
//
//      ex = faceCoords[0];
//      ex -= fc;
//      ex /= ex.L2_Norm();
//
//      ey = faceCoords[1];
//      ey -= fc;
//      ey /= ey.L2_Norm();
//      ez.Cross(ex, ey);
//      if (ez.L2_Norm() < 0.01)  // Node 0, 1, and face center are roughly
// along one straight line.  We use another node to construct the vectors.
//      {
//        ey = faceCoords[2];
//        ey -= fc;
//        ey /= ey.L2_Norm();
//      }
//
//      ez.Cross(ex, ey); ez /= ez.L2_Norm();
//      ey.Cross(ez,ex); ey /= ey.L2_Norm();
//  }
  else
  {
    ez = fc;
    ez -=elementCenter;

    /// Approximate in-plane axis
    ex = faceCoords[0];
    ex -= fc;
    ex /= ex.L2_Norm();
    ey.Cross(ez,ex);
    ey /= ey.L2_Norm();
  }


  if (numFaceNodes > 2)
  {
    /// Sort nodes counterclockwise around face center
    array1d< std::pair<realT,int> > thetaOrder(numFaceNodes);
    for( localIndex n =0 ; n < numFaceNodes ; ++n)
    {
      R1Tensor v = faceCoords[n];
      v -= fc;
      thetaOrder[n] = std::pair<realT,int>(atan2(Dot(v,ey),Dot(v,ex)),faceNodes[n]);
    }

    sort(thetaOrder.begin(), thetaOrder.end());

    // Reorder nodes on face
    for( localIndex n =0 ; n < numFaceNodes ; ++n)
    {
      faceNodes[n] = thetaOrder[n].second;
    }

    localIndex_array tempFaceNodes(numFaceNodes);
    localIndex firstIndexIndex = 0;
    for( localIndex n =0 ; n < numFaceNodes ; ++n)
    {
      tempFaceNodes[n] = thetaOrder[n].second;
      if( tempFaceNodes[n] == firstNodeIndex )
      {
        firstIndexIndex = n;
      }
    }
    for( localIndex n=0 ; n < numFaceNodes ; ++n)
    {
      const localIndex index = firstIndexIndex+n < numFaceNodes ? firstIndexIndex+n : firstIndexIndex+n-numFaceNodes;
      faceNodes[n] = tempFaceNodes[index];
    }

  }
  else
  {
    ez.Cross(ex, ey);
    // The element should be on the right hand side of the vector from node 0 to
    // node 1.
    // This ensure that the normal vector of an external face points to outside
    // the element.
    if (ez[2] > 0)
    {
      localIndex itemp = faceNodes[0];
      faceNodes[0] = faceNodes[1];
      faceNodes[1] = itemp;

    }
  }

}


void FaceManager::ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const & nodeManager,
                                                                   array1d<globalIndex_array>& faceToNodes )
{

  OrderedVariableOneToManyRelation const & faceNodes = this->nodeList();
  integer_array const & isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);

  nodeManager.CheckTypeID( typeid( NodeManager ) );


  faceToNodes.clear();
  for( localIndex kf=0 ; kf<size() ; ++kf )
  {

    if( isDomainBoundary(kf) != 0 )
    {
      globalIndex_array temp;

      for( localIndex a=0 ; a<faceNodes[kf].size() ; ++a )
      {
        const globalIndex gnode = nodeManager.m_localToGlobalMap( faceNodes[kf][a] );
        temp.push_back( gnode );
      }
      std::sort( temp.begin(), temp.end() );
//      temp.insert( temp.begin(), this->m_localToGlobalMap[kf] );
      faceToNodes.push_back(temp);
    }
  }
}



void FaceManager::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementSubRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementListString));
}


localIndex FaceManager::PackUpDownMapsSize( localIndex_array const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

localIndex FaceManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                 localIndex_array const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template<bool DOPACK>
localIndex FaceManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                        localIndex_array const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::nodeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                             m_nodeList,
                                             packList,
                                             this->m_localToGlobalMap,
                                             m_nodeList.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::edgeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                             m_edgeList,
                                             packList,
                                             this->m_localToGlobalMap,
                                             m_edgeList.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::elementListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         this->m_toElements,
                                         packList,
                                         m_toElements.getElementRegionManager() );


  return packedSize;
}



localIndex FaceManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                 localIndex_array const & packList )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF( nodeListString==viewKeyStruct::nodeListString, "");

  unPackedSize += bufferOps::Unpack( buffer,
                                         m_nodeList,
                                         packList,
                                         this->m_globalToLocalMap,
                                         m_nodeList.RelatedObjectGlobalToLocal() );

  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF( edgeListString==viewKeyStruct::edgeListString, "");

  unPackedSize += bufferOps::Unpack( buffer,
                                         m_edgeList,
                                         packList,
                                         this->m_globalToLocalMap,
                                         m_edgeList.RelatedObjectGlobalToLocal() );


  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF( elementListString==viewKeyStruct::elementListString, "");

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager() );



  return unPackedSize;
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, FaceManager, std::string const &, ManagedGroup * const )

}
