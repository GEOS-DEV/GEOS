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

#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/FiniteElementSpace.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"

namespace geosx
{

//using namespace raja;
  
  
//-------------
//Here we want to unpack data and 
//use one of the templated loops above
//------------
template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_nodes( MeshLevel const * const mesh, LAMBDA && body)
{  
  NodeManager const * const nodeManager = mesh->getNodeManager();  
  ::geosx::forall_in_range<POLICY> (0, nodeManager->size(), body);
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_faces( MeshLevel const * const mesh, LAMBDA && body)
{
  FaceManager const * const faceManager = mesh->getFaceManager();
  ::geosx::forall_in_range<POLICY> (0, faceManager->size(), body);
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_faces( MeshLevel const * const mesh, const localIndex *setList, localIndex listLen, LAMBDA && body){
  ::geosx::forall_in_set<POLICY> (setList, listLen, body);
}


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
      
      ::geosx::forall_in_range<POLICY>(0,cellBlock->size(), body);

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

      ::geosx::forall_in_set<POLICY>(setList, listLen, body);
    }
  }
}

template<class POLICY=elemPolicy,typename LAMBDA=void>
void for_elems_in_subRegion( CellBlockSubRegion const * const subRegion, LAMBDA && body)
{
  ::geosx::forall_in_range<POLICY>(0,subRegion->size(), body);
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
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);

      forall_in_range<POLICY>(0, cellBlockSubRegion->size(),
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
  
  ::geosx::forall_in_range(begin, end, GEOSX_LAMBDA (localIndex index) mutable -> void
  {      
      sum += body(index);
  });

  
  return sum.get();
}


template<typename NUMBER=real64,class EXEC_POLICY=elemPolicy,class REDUCE_POLICY=reducePolicy,typename LAMBDA=void>
NUMBER sum_in_set(localIndex const * const indexList, const localIndex len, LAMBDA && body)
{
  ReduceSum<REDUCE_POLICY, NUMBER> sum(NUMBER(0));
  ::geosx::forall_in_set(indexList, GEOSX_LAMBDA (localIndex index) mutable -> void
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



