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

#ifndef __GEOS_RAJA_WRAPPER__HPP
#define __GEOS_RAJA_WRAPPER__HPP

#include "../rajaInterface/GEOS_RAJA_Policies.hpp"

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/FiniteElementSpace.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"


#define GEOSX_LAMBDA [=]

namespace geosx{

/*
template<class POLICY=elemPolicy,typename LAMBDA=void>
void forall_in_range( localIndex const begin, const localIndex end, LAMBDA && body)
{
  RAJA::RangeSegment seg(begin, end);
  RAJA::forall<POLICY>( seg , [=] (localIndex index) mutable -> void
  {
    body(index);
  } );
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void forall_in_set(localIndex const * const indexList, const localIndex len, LAMBDA && body)
{
  RAJA::TypedListSegment<localIndex> listSeg(indexList, len, RAJA::Unowned);
  RAJA::forall<POLICY>( listSeg , [=] (localIndex index) mutable -> void
  {
    body(index);
  } );
}

*/

//================
//Hevy Computation
//================
#define FORALL( INDEX, begin, end )  \
  forall_in_range<RAJA::omp_parallel_for_exec>(begin, end,   \
    GEOSX_LAMBDA ( localIndex const INDEX ) mutable -> void

#define FORALL_IN_SET(INDEX, arr, arrLen)                 \
  forall_in_set<RAJA::omp_parallel_for_exec> (arr, arrLen, \
  GEOSX_LAMBDA ( localIndex const INDEX ) mutable -> void

//======================
//Streaming Computations
//======================
#define FORALL_NODES(INDEX, begin, end) \
  forall_in_range<RAJA::loop_exec>(begin, end,   \
  GEOSX_LAMBDA ( localIndex const INDEX ) mutable -> void

#define FORALL_NODES_IN_SET(INDEX, arr, arrLen) \
  forall_in_set<RAJA::loop_exec> (arr, arrLen,  \
  GEOSX_LAMBDA ( localIndex const INDEX ) mutable -> void

  
//Randy's code
//template<class POLICY=elemPolicy,typename LAMBDA=void>
//void forall_in_range( dataRepository::ManagedGroup const * const group, LAMBDA && body)
//{
//FORALL(a, 0, group->size()){
//body(a);
//});
//}

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
void for_nodes( MeshLevel const * const mesh, const localIndex *setList, localIndex listLen, LAMBDA && body){
  forall_in_set<POLICY> (setList, listLen, body);
}


template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_faces( MeshLevel const * const mesh, LAMBDA && body)
{
  FaceManager const * const faceManager = mesh->getFaceManager();
  forall_in_range<POLICY> (0, faceManager->size(), body);
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_faces( MeshLevel const * const mesh, const localIndex *setList, localIndex listLen, LAMBDA && body){
  forall_in_set<POLICY> (setList, listLen, body);
}

//How to unpack these?
//template<class POLICY=elemPolicy,typename LAMBDA=void>
//void for_edges( MeshLevel const * const mesh, LAMBDA && body);
 

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_elems( MeshLevel const * const mesh, LAMBDA && body)
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  dataRepository::ManagedGroup const * const elementRegions = elemManager->GetGroup(dataRepository::keys::elementRegions);

  for( auto & region : elementRegions->GetSubGroups() )
  {
    dataRepository::ManagedGroup const * const cellBlockSubRegions = region.second->GetGroup(dataRepository::keys::cellBlockSubRegions);
    for( auto const & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
    {
      CellBlockSubRegion const * const cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);
      
      forall_in_range<POLICY>(0,cellBlock->size(), body);

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
    dataRepository::ManagedGroup const * const cellBlockSubRegions = region.second->GetGroup(dataRepository::keys::cellBlockSubRegions);
    for( auto & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
    {
      CellBlockSubRegion const * const cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);

      forall_in_set<POLICY>(setList, listLen, body);
    }
  }
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_elems_in_subRegion( CellBlockSubRegion const * const subRegion, LAMBDA && body)
{
  forall_in_range<POLICY>(0,subRegion->size(), body);
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_elems_by_constitutive( MeshLevel const * const mesh,
                               constitutive::ConstitutiveManager * const constitutiveManager,
                               FiniteElementSpaceManager const * const feSpaceManager,
                               LAMBDA && body )
{
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  dataRepository::ManagedGroup const * const elementRegions = elemManager->GetGroup(dataRepository::keys::elementRegions);


  for( auto const & regionPair : elementRegions->GetSubGroups() )
  {
    dataRepository::ManagedGroup const * const elementRegion = regionPair.second;
    auto const & numMethodName = elementRegion->getData<string>(dataRepository::keys::numericalMethod);
    FiniteElementSpace const * const feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);

    dataRepository::ManagedGroup const * const cellBlockSubRegions = elementRegion->GetGroup(dataRepository::keys::cellBlockSubRegions);
    for( auto & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
    {
      CellBlockSubRegion const * cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);

      //auto const & dNdX = cellBlock->getData< multidimensionalArray::ManagedArray< R1Tensor, 3 > >(keys::dNdX);
      multidimensionalArray::ManagedArray<R1Tensor, 3> const & dNdX = cellBlock->getReference< multidimensionalArray::ManagedArray<R1Tensor, 3> >(dataRepository::keys::dNdX);
      
      array_view<real64,2> const & detJ            = cellBlock->getReference< array2d<real64> >(dataRepository::keys::detJ).View();

      auto const & constitutiveMap = cellBlock->getReference< std::pair< array2d<localIndex>,array2d<localIndex> > >(CellBlockSubRegion::viewKeyStruct::constitutiveMapString);
//      RAJA::View< localIndex const, RAJA::Layout<2> > constitutiveMapView( reinterpret_cast<localIndex const*>(constitutiveMap.second.data()),
//                                                                           constitutiveMap.second.size(0),
//                                                                           constitutiveMap.second.size(1) );
      array_view<localIndex,2> constitutiveMapView = constitutiveMap.second.View();

      auto const & constitutiveGrouping = cellBlock->getReference< map< string, localIndex_array > >(CellBlockSubRegion::viewKeyStruct::constitutiveGroupingString);
      array_view<localIndex,2> const elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference().View();// getData<array2d<localIndex>>(keys::nodeList);

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      for( auto const & constitutiveGroup : constitutiveGrouping )
      {
        string const constitutiveName = constitutiveGroup.first;
//      localIndex_array const & elementList = constitutiveGroup.second;
        array_view<localIndex,1> const elementList = constitutiveGroup.second.View();

        constitutive::ConstitutiveBase * constitutiveModel = constitutiveManager->GetGroup<constitutive::ConstitutiveBase>( constitutiveName );

        void * constitutiveModelData;
        constitutiveModel->SetParamStatePointers(constitutiveModelData);

        constitutive::ConstitutiveBase::UpdateFunctionPointer
        constitutiveUpdate = constitutiveModel->GetStateUpdateFunctionPointer();

        ///Local copies --------
        array_view<real64,1> meanStress    = constitutiveModel->GetGroup(std::string("StateData"))->getReference<real64_array>(std::string("MeanStress")).View();
        array_view<R2SymTensor,1> devStress    = constitutiveModel->GetGroup(std::string("StateData"))->getReference<r2Sym_array>(std::string("DeviatorStress")).View();
        //------------------------

        //Element loop is packed with parameters...
        auto ebody = [=](localIndex index) mutable -> void
          {body(index,
                numNodesPerElement,
                elemsToNodes,
                feSpace->m_finiteElement->n_quadrature_points(),
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

#define FOR_ELEMS_FOR_CONSTITUTIVE( mesh, constitutiveManager, feSpaceManager)\
    for_elems_by_constitutive( mesh,\
    constitutiveManager,\
    feSpaceManager,\
    GEOSX_LAMBDA( localIndex const k,\
    localIndex const numNodesPerElement,\
    array_view<localIndex,2> const elemsToNodes,\
    localIndex const numQuadraturePoints,\
    multidimensionalArray::ManagedArray<R1Tensor, 3> const & dNdX,\
    array_view<localIndex,2> const constitutiveMapView,\
    array_view<real64,2> const detJ,\
    array_view<R2SymTensor,1> devStress,\
    array_view<real64,1> meanStress,\
    constitutive::ConstitutiveBase::UpdateFunctionPointer constitutiveUpdate,\
    void * constitutiveModelData\
    ) mutable -> void




template<class POLICY=elemPolicy,typename LAMBDA=void>
void forAllElemsInMesh( MeshLevel const * const mesh, LAMBDA && lambdaBody)
{

  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);

      auto ebody = [=](localIndex index) mutable -> void
      {
        lambdaBody(er,esr,index);
      };

      forall_in_range<POLICY>(0, cellBlockSubRegion->size(), ebody);
    }
  }
}

template<typename NUMBER=real64,class EXEC_POLICY=elemPolicy,class REDUCE_POLICY=reducePolicy,typename LAMBDA=void>
NUMBER sum_in_range(localIndex const begin, const localIndex end, LAMBDA && body)
{
  RAJA::ReduceSum<REDUCE_POLICY, NUMBER> sum(NUMBER(0));
  RAJA::RangeSegment seg(begin, end);
  RAJA::forall<EXEC_POLICY>( seg , [=] (localIndex index) mutable -> void
  {
    sum += body(index);
  } );
  return sum.get();
}

template<typename NUMBER=real64,class EXEC_POLICY=elemPolicy,class REDUCE_POLICY=reducePolicy,typename LAMBDA=void>
NUMBER sum_in_set(localIndex const * const indexList, const localIndex len, LAMBDA && body)
{
  RAJA::ReduceSum<REDUCE_POLICY, NUMBER> sum(NUMBER(0));
  RAJA::TypedListSegment<localIndex> listSeg(indexList, len, RAJA::Unowned);
  RAJA::forall<EXEC_POLICY>( listSeg , [=] (localIndex index) mutable -> void
  {
    sum += body(index);
  } );
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
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);

      auto ebody = [=](localIndex index) mutable -> real64
      {
        return lambdaBody(er,esr,index);
      };

      sum += sum_in_range<real64,EXEC_POLICY,REDUCE_POLICY>(0, cellBlockSubRegion->size(), ebody);
    }
  }

  return sum;
}


}


#define FOR_ELEMS_IN_SUBREGION( SUBREGION, INDEX )  \
    for_elems_in_subRegion( SUBREGION, GEOSX_LAMBDA ( localIndex const INDEX ) mutable -> void
#endif

#define END_FOR );



