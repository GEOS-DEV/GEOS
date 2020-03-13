/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MeshUtilities.cpp
 *
 */

#include "MeshUtilities.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "SimpleGeometricObjects/SimpleGeometricObjectBase.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/NodeManager.hpp"

namespace geosx
{
using namespace dataRepository;

MeshUtilities::MeshUtilities()
{
  // TODO Auto-generated constructor stub

}

MeshUtilities::~MeshUtilities()
{
  // TODO Auto-generated destructor stub
}



void MeshUtilities::GenerateNodesets( dataRepository::Group const * geometries,
                                      NodeManager * const nodeManager )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager->referencePosition();
  localIndex const numNodes = nodeManager->size();
  Group * sets = nodeManager->sets();

  for( int i = 0; i < geometries->GetSubGroups().size(); ++i )
  {
    SimpleGeometricObjectBase const * const object = geometries->GetGroup< SimpleGeometricObjectBase >( i );
    if( object!=nullptr )
    {
      string name = object->getName();
      SortedArray< localIndex > & targetSet = sets->registerWrapper< SortedArray< localIndex > >( name )->reference();
      for( localIndex a=0; a<numNodes; ++a )
      {
        if( object->IsCoordInObject( X[a] ))
        {
          targetSet.insert( a );
        }
      }
    }

  }
}
//
//void MeshUtilities::GenerateFasesetsAndAssociatedNodesets(
// TICPP::HierarchicalDataNode& hdn,
//                                                           FaceManagerT *
// faceManager,
//                                                           NodeManager *
// nodeManager )
//{
//
////  std::map< std::string, SortedArray<localIndex> >& nodeSets = nodeManager->m_Sets;
//  std::map< std::string, SortedArray<localIndex> >& faceSets = faceManager->m_Sets;
////  array1d<R1Tensor>& X = *(nodeManager->m_refposition);
//
//  //We calculate face centers here. This is cheaper than calculating it when
// we loop through faces.
//
//  faceManager->AddKeylessDataField( FieldInfo::R1TensorField, "FaceCenter",
// true, true );
//  faceManager->AddKeylessDataField( FieldInfo::R1TensorField, "FaceNormal",
// true, true );
//  array1d<R1Tensor>& faceCenter = faceManager->GetFieldData<R1Tensor>(
// "FaceCenter" );
//  array1d<R1Tensor>& faceNormal = faceManager->GetFieldData<R1Tensor>(
// "FaceNormal" );
//
//  for( localIndex kf=0 ; kf<faceManager->DataLengths() ; ++kf )
//  {
//    faceManager->FaceCenter( nodeManager, kf, faceCenter[kf] );
//    faceManager->FaceNormal( nodeManager, kf, faceNormal[kf] );
//  }
//
//  for(TICPP::HierarchicalDataNode* hdnNode = hdn.Next(true); hdnNode; hdnNode
// = hdn.Next() )
//  {
//    std::string header = "Faceset";
//    if( !( header.compare(hdnNode->Heading() ) ) )
//    {
//      /*
//      SimpleGeometricObjectBase::Types type =
// SimpleGeometricObjectBase::IntToType(hdnNode->GetAttributeValue<int>("type")
// );
//      */
//
//      std::string name = hdnNode->GetAttributeString("name");
////      SortedArray<localIndex>& currentNodeset = nodeSets[name];
//      SortedArray<localIndex>& currentFaceset = faceSets[name];
//
//      SimpleGeometricObjectBase* object;
//
//      // new allocation method
//      TICPP::HierarchicalDataNode* geometryNode =
// hdnNode->GetChild("Geometry");
//      if(geometryNode){
//      // treating geometry node as a union boolean geometry allows multiple
// objects to be defined within it:
//        /**
//          <Geometry>
//              <Box  ... />
//              <Sphere .. />
//              <Intersection>
//                <Cylinder .. />
//                <Not>
//                  <Cylinder .. />
//                </Not>
//              </Intersection>
//          <Geometry>
//          **/
//        object =
// SimpleGeometricObjectBase::Allocate(SimpleGeometricObjectBase::unionGeometry);
//        object->ReadXML( *geometryNode );
//      } else {
//        // old allocation method for backwards compatability
//        std::string geometricObjectTypeStr =
// hdnNode->GetAttributeStringOrDefault("type", "Box");
//        SimpleGeometricObjectBase::Types type =
// fromString<SimpleGeometricObjectBase::Types>(geometricObjectTypeStr);
//        object = SimpleGeometricObjectBase::Allocate( type );
//
//        object->ReadXML( *hdnNode );
//      }
//
//      R1Tensor zeroVector;
//      zeroVector *= 0.0;
//      R1Tensor normalVector =
// hdnNode->GetAttributeOrDefault<R1Tensor>("normalVector", zeroVector);
//      if (normalVector.Normalize() == 0.0)
//        throw GPException("A normal vector must be given to define face
// sets");
//
//      realT withinAngle = hdnNode->GetAttributeOrDefault<realT>("withinAngle",
// 1.0);
//      realT tol = std::cos(withinAngle / 180.0 * 3.14159265);
//
//      for( localIndex kf=0 ; kf<faceManager->DataLengths() ; ++kf )
//      {
//        if( object->IsCoordInObject( faceCenter[kf] ) &&
// std::fabs(Dot(faceNormal[kf], normalVector)) > tol)
//        {
//          currentFaceset.insert(kf);
////          for( localIndex a=0 ; a<faceManager->m_toNodesRelation[kf].size()
// ; ++a )
////            currentNodeset.insert(faceManager->m_toNodesRelation[kf][a]);
//        }
//      }
//      SimpleGeometricObjectBase::Deallocate( object );
//    }
//  }
//}
//
//
//void MeshUtilities::GenerateElementsets( TICPP::HierarchicalDataNode& hdn,
//                                         const NodeManager& nodeManager,
//                                      ElementManagerT& elementManager )
//{
//  for(TICPP::HierarchicalDataNode* hdnNode = hdn.Next(true); hdnNode; hdnNode
// = hdn.Next() )
//  {
//    std::string header = "Elementset";
//    if( !( header.compare(hdnNode->Heading() ) ) )
//    {
//      SimpleGeometricObjectBase* object;
//
//      // new allocation method
//      TICPP::HierarchicalDataNode* geometryNode =
// hdnNode->GetChild("Geometry");
//      if(geometryNode){
//      // treating geometry node as a union boolean geometry allows multiple
// objects to be defined within it:
//        /**
//          <Geometry>
//              <Box  ... />
//              <Sphere .. />
//              <Intersection>
//                <Cylinder .. />
//                <Not>
//                  <Cylinder .. />
//                </Not>
//              </Intersection>
//          <Geometry>
//          **/
//        object =
// SimpleGeometricObjectBase::Allocate(SimpleGeometricObjectBase::unionGeometry);
//        object->ReadXML( *geometryNode );
//      }
//      else
//      {
//        std::string geometricObjectTypeStr =
// hdnNode->GetAttributeStringOrDefault("type", "Box");
//        SimpleGeometricObjectBase::Types type =
// fromString<SimpleGeometricObjectBase::Types>(geometricObjectTypeStr);
//        object = SimpleGeometricObjectBase::Allocate( type );
//
//        object->ReadXML( *hdnNode );
//      }
//
//      for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::iterator
// elementRegionIter =
//          elementManager.m_ElementRegions.begin(); elementRegionIter !=
// elementManager.m_ElementRegions.end(); ++elementRegionIter) //Loop over
// regions
//      {
//        ElementRegionT& elemRegion = elementRegionIter->second;
//        localIndex numEle = elemRegion.DataLengths();
//
//        std::map< std::string, SortedArray<localIndex> >& sets = elemRegion.m_Sets;
//
//        std::string name = hdnNode->GetAttributeString("name");
//        SortedArray<localIndex>& set = sets[name];
//
//        for (localIndex iElm=0 ; iElm<numEle ; ++iElm)
//        {
//          auto elmCenter = elemRegion.GetElementCenter(iElm, nodeManager);
//
//          if( object->IsCoordInObject( elmCenter ) )
//          {
//            set.insert(iElm);
//          }
//        }
//      }
//      SimpleGeometricObjectBase::Deallocate( object );
//    }
//  }
//}
}
