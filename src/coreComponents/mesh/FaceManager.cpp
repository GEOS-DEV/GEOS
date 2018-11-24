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
#include "common/TimingMacros.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

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
  m_toElements.m_toElementRegion.setDefaultValue(-1);
  m_toElements.m_toElementSubRegion.setDefaultValue(-1);
  m_toElements.m_toElementIndex.setDefaultValue(-1);

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

  docNode->AllocateChildNode( viewKeyStruct::faceAreaString,
                              viewKeyStruct::faceAreaString,
                              -1,
                              "real64_array",
                              "real64_array",
                              "Face surface area",
                              "Face surface area",
                              "",
                              this->getName(),
                              1,
                              0,
                              1 );

  docNode->AllocateChildNode( viewKeyStruct::faceCenterString,
                              viewKeyStruct::faceCenterString,
                              -1,
                              "r1_array",
                              "r1_array",
                              "Face centroid coordinates",
                              "Face centroid coordinates",
                              "",
                              this->getName(),
                              1,
                              0,
                              1 );

  docNode->AllocateChildNode( viewKeyStruct::faceNormalString,
                              viewKeyStruct::faceNormalString,
                              -1,
                              "r1_array",
                              "r1_array",
                              "Face normal ",
                              "Face normal",
                              "",
                              this->getName(),
                              1,
                              0,
                              1 );


}

void FaceManager::BuildFaces( NodeManager * const nodeManager, ElementRegionManager * const elementManager )
{
  GEOSX_MARK_FUNCTION;
  localIndex numFaces = 0;
  localIndex_array tempNodeList;
  array1d<localIndex_array> facesByLowestNode(nodeManager->size());

  array2d<localIndex> & elemRegionList = elementRegionList();
  array2d<localIndex> & elemSubRegionList = elementSubRegionList();
  array2d<localIndex> & elemList = elementList();
  OrderedVariableOneToManyRelation & node_list = nodeList();

  m_toElements.setElementRegionManager( elementManager );

  elemRegionList.resize( 20*nodeManager->size() );
  elemSubRegionList.resize( 20*nodeManager->size() );
  elemList.resize( 20*nodeManager->size() );  // We need to reserve a lot more space for tets.  These get resized later.
  node_list.resize( 20*nodeManager->size() );

  elemRegionList = -1;
  elemSubRegionList = -1;
  elemList = -1;

  for( typename dataRepository::indexType kReg=0 ; kReg<elementManager->numRegions() ; ++kReg  )
  {
    ElementRegion * const elemRegion = elementManager->GetRegion(kReg);

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion * const subRegion = elemRegion->GetSubRegion(kSubReg);
      localIndex const numFacesPerElement = subRegion->numFacesPerElement();
      array2d<localIndex> & elemsToFaces = subRegion->faceList();

      for( localIndex ke=0 ; ke<subRegion->size() ; ++ke )
      {
        // kelf = k'th element local face index
        for( localIndex kelf=0 ; kelf< numFacesPerElement ; ++kelf )
        {
          // get the nodes associated with the local face
          subRegion->GetFaceNodes( ke, kelf, tempNodeList );

          //Special treatment for the triangle faces of prisms.
          if (tempNodeList[tempNodeList.size() - 1] == std::numeric_limits<localIndex>::max())
          {
            tempNodeList.pop_back();
          }

          // sort the nodes
          std::sort(tempNodeList.begin(), tempNodeList.end());

          // get the lowest node index from the list for simplicity
          const localIndex lowNode = tempNodeList[0];

          // Check if the face has already been added
          bool duplicate = false;

          // Iterate through the possible matches
          localIndex_array& matching_faces = facesByLowestNode[lowNode];
          const localIndex n_matches = matching_faces.size();
          for( localIndex cur_match = 0; cur_match < n_matches; ++cur_match )
          {
            const localIndex existingFaceIndex = matching_faces[cur_match];

            // this is the nodelist of the face that we are testing against
            const localIndex_array& existingFaceNodelist = node_list[existingFaceIndex];

            // since the size is the same, then we should test the
            // nodes...they are sorted, so
            // the std::equal() algorithm will work for this.
            if( existingFaceNodelist.size() == tempNodeList.size() && 
                std::equal( existingFaceNodelist.begin(), existingFaceNodelist.end(), tempNodeList.begin() ) )
            {
              // now remove the entry from the face that we were checking
              // against from the facesByLowestNode list...
              // because it is no longer possible that it will have
              // another element that has this face.
              matching_faces[cur_match] = matching_faces.back();
              matching_faces.pop_back();

              // they are equal!
              duplicate = true;

              // add the element to the faceToElement map
              elemRegionList[existingFaceIndex][1]    = kReg;
              elemSubRegionList[existingFaceIndex][1] = kSubReg;
              elemList[existingFaceIndex][1]          = ke;

              // add the face to the elementToFaceMap for the element
              // region.
              elemsToFaces[ke][kelf] = existingFaceIndex;

              // break the loop
              break;
            }
          }
          
          if( !duplicate )
          {
            // the face is not a duplicate of any in the facesByLowestNode
            // list, so we need to add a new face.

            const localIndex curFace = numFaces++;

            // add the nodes to the faceToNodeMap
            node_list[curFace] = tempNodeList;

//            GEOS_LOG_RANK( "face #"<<curFace<<" ("<<tempNodeList<<")");

            // and add the face to facesByLowestNode[]
            matching_faces.push_back(curFace);

            elemRegionList[curFace][0] = kReg;
            elemSubRegionList[curFace][0] = kSubReg;
            elemList[curFace][0] = ke;

            // add the face to the elementToFaceMap
            elemsToFaces[ke][kelf] = curFace;
          }
        }
      }
    }
  }

  // resize the data vectors according to the number of faces
  resize(numFaces);

  // make sets from nodesets
  // First create the sets
  auto const & nodeSets = nodeManager->GetGroup(string("Sets"))->wrappers();
  for ( int i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    CreateSet( setName );
  }

  // Then loop over them in parallel.
  forall_in_range<parallelHostPolicy>( 0, nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    const set<localIndex>& targetSet = nodeManager->GetGroup(keys::sets)->getReference<set<localIndex>>( setName );
    ConstructSetFromSetAndMap( targetSet, m_nodeList, setName );
  } );

  // sort the face node lists
  SortAllFaceNodes( nodeManager, elementManager);

  SetDomainBoundaryObjects( nodeManager );






  real64_array & faceArea  = getReference<real64_array>( viewKeyStruct::
                                                         faceAreaString);

  r1_array & faceNormal = getReference<r1_array>( viewKeyStruct::
                                                     faceNormalString);

  r1_array & faceCenter = getReference<r1_array>( viewKeyStruct::
                                                      faceCenterString);

  r1_array const & X = nodeManager->referencePosition();


  // loop over faces and calculate faceArea, faceNormal and faceCenter
  for (localIndex kf = 0; kf < this->size(); ++kf)
  {
    faceArea[kf] = computationalGeometry::Centroid_3DPolygon(m_nodeList[kf],
                                                             X,
                                                             faceCenter[kf],
                                                             faceNormal[kf]);
  }

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

  forall_in_range<parallelHostPolicy>( 0, size(), [&]( localIndex const kf ) -> void
  {
    if( elemRegionList[kf][1] == -1 )
    {
      faceDomainBoundaryIndicator(kf) = 1;
    }
  } );

  integer_array & nodeDomainBoundaryIndicator = nodeManager->getReference<integer_array>(nodeManager->viewKeys.domainBoundaryIndicator);
  nodeDomainBoundaryIndicator = 0;

  OrderedVariableOneToManyRelation const & faceToNodesMap = this->nodeList();

  forall_in_range<parallelHostPolicy>( 0, size(), [&]( localIndex const k ) mutable -> void
  {
    if( faceDomainBoundaryIndicator[k] == 1 )
    {
      arrayView1d<localIndex> const& nodelist = faceToNodesMap[k];
      for( localIndex a=0 ; a< nodelist.size() ; ++a )
      {
        nodeDomainBoundaryIndicator[nodelist[a]] = 1;
      }
    }
  } );
}


void FaceManager::SetIsExternal()
{
  integer_array const &
  isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);

  m_isExternal = 0;
  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( isDomainBoundary[k]==1 )
    {
      m_isExternal[k] = 1;
    }
  }
}

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

localIndex FaceManager::getMaxFaceNodes() const
{
  localIndex max_size = 0;
  OrderedVariableOneToManyRelation const & faces_to_nodes = nodeList();
  for(localIndex kf =0 ; kf < size() ; ++kf )
  {
    const localIndex size = faces_to_nodes[kf].size();
    if ( size > max_size )
    {
      max_size = size;
    }
  }

  return max_size;
}

void FaceManager::SortAllFaceNodes( NodeManager const * const nodeManager,
                                    ElementRegionManager const *const elemManager )
{
  array2d<localIndex> const & elemRegionList = elementRegionList();
  array2d<localIndex> const & elemSubRegionList = elementSubRegionList();
  array2d<localIndex> const & elemList = elementList();
  r1_array const & X = nodeManager->referencePosition();

  const indexType max_face_nodes = getMaxFaceNodes();
  GEOS_ERROR_IF( max_face_nodes >= MAX_FACE_NODES, "More nodes on a face than expected!" );

  forall_in_range<parallelHostPolicy>( 0, size(), [&]( localIndex const kf ) -> void
  {
    ElementRegion const * const elemRegion     = elemManager->GetRegion( elemRegionList[kf][0] );
    CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion( elemSubRegionList[kf][0] );
    R1Tensor const elementCenter = subRegion->GetElementCenter( elemList[kf][0], *nodeManager );
    
    const localIndex numFaceNodes = nodeList()[kf].size();
    arrayView1d<localIndex> & faceNodes = nodeList()[kf];
    SortFaceNodes( X, elementCenter, faceNodes, numFaceNodes );
  } );
}

void FaceManager::SortFaceNodes( arrayView1d<R1Tensor> const & X,
                                 R1Tensor const & elementCenter,
                                 arrayView1d<localIndex> & faceNodes,
                                 localIndex const numFaceNodes )
{
  localIndex const firstNodeIndex = faceNodes[0];

  // get face center (average vertex location) and store node coordinates
  R1Tensor const * face_coords[MAX_FACE_NODES];
  R1Tensor fc(0);
  for( localIndex n =0 ; n < numFaceNodes ; ++n)
  {
    localIndex nd = faceNodes[n];
    face_coords[n] = &(X[nd]);
    fc += *(face_coords[n]);
  }
  fc /= realT(numFaceNodes);

  R1Tensor ex, ey, ez;
  // Approximate face normal direction (unscaled)

  if (numFaceNodes == 2)  //2D only.
  {
    ex = X[faceNodes[1]];
    ex -= X[faceNodes[0]];
    ey = elementCenter;
    ey -= fc;

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
  else
  {
    ez = fc;
    ez -= elementCenter;

    /// Approximate in-plane axis
    ex = *(face_coords[0]);
    ex -= fc;
    ex /= ex.L2_Norm();
    ey.Cross(ez, ex);
    ey /= ey.L2_Norm();

    std::pair<realT, localIndex> thetaOrder[MAX_FACE_NODES];

    /// Sort nodes counterclockwise around face center
    for( localIndex n =0 ; n < numFaceNodes ; ++n)
    {
      R1Tensor v = *(face_coords[n]);
      v -= fc;
      thetaOrder[n] = std::pair<realT, localIndex>(atan2(Dot(v,ey),Dot(v,ex)),faceNodes[n]);
    }

    sort(thetaOrder, thetaOrder + numFaceNodes);

    // Reorder nodes on face
    for( localIndex n =0 ; n < numFaceNodes ; ++n)
    {
      faceNodes[n] = thetaOrder[n].second;
    }

    localIndex tempFaceNodes[MAX_FACE_NODES];

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
}


void FaceManager::ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const nodeManager,
                                                                   array1d<globalIndex_array>& faceToNodes )
{

  OrderedVariableOneToManyRelation const & faceNodes = this->nodeList();
  integer_array const & isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);

  nodeManager->CheckTypeID( typeid( NodeManager ) );


  faceToNodes.clear();
  faceToNodes.resize(size());
  for( localIndex kf=0 ; kf<size() ; ++kf )
  {

    if( isDomainBoundary(kf) != 0 )
    {
      globalIndex_array temp;

      for( localIndex a=0 ; a<faceNodes[kf].size() ; ++a )
      {
        const globalIndex gnode = nodeManager->m_localToGlobalMap( faceNodes[kf][a] );
        temp.push_back( gnode );
      }
      std::sort( temp.begin(), temp.end() );
      faceToNodes[kf] = temp;
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


localIndex FaceManager::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

localIndex FaceManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template<bool DOPACK>
localIndex FaceManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d<localIndex const> const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::nodeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_nodeList.Base(),
                                         packList,
                                         this->m_localToGlobalMap,
                                         m_nodeList.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::edgeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_edgeList.Base(),
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
                                          localIndex_array & packList )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF( nodeListString != viewKeyStruct::nodeListString, "");

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_nodeList,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->m_globalToLocalMap,
                                     m_nodeList.RelatedObjectGlobalToLocal() );

  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF( edgeListString != viewKeyStruct::edgeListString, "");

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_edgeList,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->m_globalToLocalMap,
                                     m_edgeList.RelatedObjectGlobalToLocal() );


  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF( elementListString != viewKeyStruct::elementListString, "");

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     false );



  return unPackedSize;
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, FaceManager, std::string const &, ManagedGroup * const )

}
