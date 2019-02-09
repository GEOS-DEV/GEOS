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

#ifndef __GEOS_RAJA_WRAPPER__HPP
#define __GEOS_RAJA_WRAPPER__HPP

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"

namespace geosx
{
   
//-------------
//Here we want to unpack data and 
//use one of the templated loops above
//------------
template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_nodes( MeshLevel const * const mesh, LAMBDA && body)
{  
  NodeManager const * const nodeManager = mesh->getNodeManager();  
  forall_in_range<POLICY> (0, nodeManager->size(), body);
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_faces( MeshLevel const * const mesh, LAMBDA && body)
{
  FaceManager const * const faceManager = mesh->getFaceManager();
  forall_in_range<POLICY> (0, faceManager->size(), body);
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_elems( MeshLevel const * const mesh, LAMBDA && body)
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  dataRepository::ManagedGroup const * const elementRegions = elemManager->GetGroup(dataRepository::keys::elementRegions);

  for( auto & region : elementRegions->GetSubGroups() )
  {
    dataRepository::ManagedGroup const * const elementSubRegions = region.second->GetGroup(ElementRegion::viewKeyStruct::elementSubRegions);
    for( auto const & iterSubRegions : elementSubRegions->GetSubGroups() )
    {
      CellElementSubRegion const * const elementSubRegion = elementSubRegions->GetGroup<CellElementSubRegion>(iterSubRegions.first);
      
      forall_in_range<POLICY>(0,elementSubRegion->size(), body);

    }
  }
}
template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_elems( MeshLevel const * const mesh, const localIndex *setList, localIndex listLen, LAMBDA && body)
{

  ElementRegionManager const * const elemManager = mesh->getElemManager();
  dataRepository::ManagedGroup const * const elementRegions = elemManager->GetGroup(dataRepository::keys::elementRegions);
  
  for( auto const & region : elementRegions->GetSubGroups() )
    {
    dataRepository::ManagedGroup const * const elementSubRegions = region.second->GetGroup(ElementRegion::viewKeyStruct::elementSubRegions);
    for( auto & iterCellBlocks : elementSubRegions->GetSubGroups() )
    {
      CellElementSubRegion const * const elementSubRegion = elementSubRegions->GetGroup<CellElementSubRegion>(iterCellBlocks.first);

      forall_in_set<POLICY>(setList, listLen, body);
    }
  }
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_elems_in_subRegion( ElementSubRegionBase const * const subRegion, LAMBDA && body)
{
  forall_in_range<POLICY>(0,subRegion->size(), body);
}

  
template<class POLICY=elemPolicy,typename LAMBDA=void>
void forAllElemsInMesh( MeshLevel const * const mesh, LAMBDA && lambdaBody)
{

  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      ElementSubRegionBase const * const elementSubRegion = elemRegion->GetSubRegion(esr);

      forall_in_range<POLICY>(0, elementSubRegion->size(),
                              [=](localIndex index) mutable -> void
                              {
                                lambdaBody(er,esr,index);
                              });
    }
  }
}


template<typename NUMBER=real64,class EXEC_POLICY=elemPolicy,class REDUCE_POLICY=reducePolicy,typename LAMBDA=void>
NUMBER sum_in_range(localIndex const begin, const localIndex end, LAMBDA && body)
{
  ReduceSum<REDUCE_POLICY, NUMBER> sum(NUMBER(0));
  
  forall_in_range(begin, end, GEOSX_LAMBDA (localIndex index) mutable -> void
  {
      sum += body(index);
  });
  
  return sum.get();
}


template<typename NUMBER=real64,class EXEC_POLICY=elemPolicy,class REDUCE_POLICY=reducePolicy,typename LAMBDA=void>
NUMBER sum_in_set(localIndex const * const indexList, const localIndex len, LAMBDA && body)
{
  ReduceSum<REDUCE_POLICY, NUMBER> sum(NUMBER(0));
  forall_in_set(indexList, GEOSX_LAMBDA (localIndex index) mutable -> void
   {
     sum += body(index);
   });
  
  return sum.get();
}

template<class EXEC_POLICY=elemPolicy, class REDUCE_POLICY=reducePolicy, typename LAMBDA=void>
real64 sumOverElemsInMesh( MeshLevel const * const mesh, LAMBDA && lambdaBody)
{
  real64 sum = 0.0;

  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    elemRegion->forElementSubRegionsIndex( [&]( localIndex const esr, auto const * const elementSubRegion )
    {
      auto ebody = [=](localIndex index) mutable -> real64
      {
        return lambdaBody(er,esr,index);
      };

      sum += sum_in_range<real64,EXEC_POLICY,REDUCE_POLICY>(0, elementSubRegion->size(), ebody);
    });
  }

  return sum;
}


}


#define FOR_ELEMS_IN_SUBREGION( SUBREGION, INDEX )  \
    for_elems_in_subRegion( SUBREGION, GEOSX_LAMBDA ( localIndex const INDEX ) mutable -> void
#endif



//
// The following code is commented out as it may serve future purpose
//
/*
template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_elems_by_constitutive( MeshLevel const * const mesh,
                               constitutive::ConstitutiveManager * const constitutiveManager,
                               FiniteElementSpaceManager const * const feDiscretizationManager,
                               LAMBDA && body )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  dataRepository::ManagedGroup const * const elementRegions = elemManager->GetGroup(dataRepository::keys::elementRegions);


  for( auto const & regionPair : elementRegions->GetSubGroups() )
  {
    dataRepository::ManagedGroup const * const elementRegion = regionPair.second;
    auto const & numMethodName = elementRegion->getReference<string>(dataRepository::keys::numericalMethod);
    FiniteElementSpace const * const feDiscretization = feDiscretizationManager->GetGroup<FiniteElementSpace>(numMethodName);

    dataRepository::ManagedGroup const * const elementSubRegions = elementRegion->GetGroup(dataRepository::keys::elementSubRegions);
    for( auto & iterCellBlocks : elementSubRegions->GetSubGroups() )
    {
      CellBlockSubRegion const * cellBlock = elementSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);

      //auto const & dNdX = cellBlock->getData< multidimensionalArray::ManagedArray< R1Tensor, 3 > >(keys::dNdX);
      arrayView3d<R1Tensor> const & dNdX = cellBlock->getReference< array3d<R1Tensor> >(dataRepository::keys::dNdX);
      
      arrayView2d<real64> const & detJ = cellBlock->getReference< array2d<real64> >(dataRepository::keys::detJ);

      auto const & constitutiveMap = cellBlock->getReference< std::pair< array2d<localIndex>,array2d<localIndex> > >(CellBlockSubRegion::viewKeyStruct::constitutiveMapString);
//      RAJA::View< localIndex const, RAJA::Layout<2> > constitutiveMapView( reinterpret_cast<localIndex const*>(constitutiveMap.second.data()),
//                                                                           constitutiveMap.second.size(0),
//                                                                           constitutiveMap.second.size(1) );
      arrayView2d<localIndex> const & constitutiveMapView = constitutiveMap.second;

      auto const & constitutiveGrouping = cellBlock->getReference< map< string, localIndex_array > >(CellBlockSubRegion::viewKeyStruct::constitutiveGroupingString);
      arrayView2d<localIndex> const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();// getData<array2d<localIndex>>(keys::nodeList);

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      for( auto const & constitutiveGroup : constitutiveGrouping )
      {
        string const constitutiveName = constitutiveGroup.first;
//      localIndex_array const & elementList = constitutiveGroup.second;
        arrayView1d<localIndex> const & elementList = constitutiveGroup.second;

        constitutive::ConstitutiveBase * constitutiveModel = constitutiveManager->GetGroup<constitutive::ConstitutiveBase>( constitutiveName );

        void * constitutiveModelData;
        constitutiveModel->SetParamStatePointers(constitutiveModelData);

        constitutive::ConstitutiveBase::UpdateFunctionPointer
        constitutiveUpdate = constitutiveModel->GetStateUpdateFunctionPointer();

        ///Local copies --------
        arrayView1d<real64> const & meanStress = constitutiveModel->GetGroup(std::string("StateData"))->getReference<real64_array>(std::string("MeanStress"));
        arrayView1d<R2SymTensor> const & devStress = constitutiveModel->GetGroup(std::string("StateData"))->getReference<r2Sym_array>(std::string("DeviatorStress"));
        //------------------------

        //Element loop is packed with parameters...
        auto ebody = [=](localIndex index) mutable -> void
          {body(index,
                numNodesPerElement,
                elemsToNodes,
                feDiscretization->m_finiteElement->n_quadrature_points(),
                dNdX,
                constitutiveMapView,
                detJ,
                devStress,
                meanStress,
                constitutiveUpdate,
                constitutiveModelData
                ); };

        
        forall_in_set<POLICY>(elementList.data(), elementList.size(), ebody);

      }
    }
  }
}

#define FOR_ELEMS_FOR_CONSTITUTIVE( mesh, constitutiveManager, feDiscretizationManager)\
    for_elems_by_constitutive( mesh,\
    constitutiveManager,\
    feDiscretizationManager,\
    GEOSX_LAMBDA( localIndex const k,\
    localIndex const numNodesPerElement,\
    arrayView2d<localIndex> const elemsToNodes,\
    localIndex const numQuadraturePoints,\
    arrayView3d<R1Tensor> const & dNdX,\
    arrayView2d<localIndex> const constitutiveMapView,\
    arrayView2d<real64> const detJ,\
    arrayView1d<R2SymTensor> devStress,\
    arrayView1d<real64> meanStress,\
    constitutive::ConstitutiveBase::UpdateFunctionPointer constitutiveUpdate,\
    void * constitutiveModelData\
    ) mutable -> void
*/
