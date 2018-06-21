// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * BoundaryConditionManager.cpp
 *
 *  Created on: May 26, 2017
 *      Author: rrsettgast
 */

#include "BoundaryConditionManager.hpp"
#include "BoundaryConditionBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"

#include "mesh/MeshBody.hpp"

#include "finiteElement/FiniteElementManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "codingUtilities/StringUtilities.hpp"

#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
BoundaryConditionManager::BoundaryConditionManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  // TODO Auto-generated constructor stub

}


BoundaryConditionManager * BoundaryConditionManager::get()
{
  static BoundaryConditionManager bcman("bcMan",nullptr);
  return &bcman;
}


BoundaryConditionManager::~BoundaryConditionManager()
{
  // TODO Auto-generated destructor stub
}

void BoundaryConditionManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<BoundaryConditionBase> bc = BoundaryConditionBase::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup(childName, std::move(bc) );
}

void BoundaryConditionManager::ApplyBoundaryCondition( dataRepository::ManagedGroup * object,
                                                       std::string const & fieldName,
                                                       real64 const time )
{
  dataRepository::ManagedGroup const * sets = object->GetGroup(dataRepository::keys::sets);

  // iterate over all boundary conditions.
  forSubGroups<BoundaryConditionBase>( [&]( BoundaryConditionBase * bc ) -> void
    {
      if( time >= bc->GetStartTime() && time < bc->GetEndTime() && ( bc->GetFieldName()==fieldName) )
      {
        string_array setNames = bc->GetSetNames();
        for( auto & setName : setNames )
        {
          dataRepository::ViewWrapper<lSet> const * const setWrapper = sets->getWrapper<lSet>(setName);
          if( setWrapper != nullptr )
          {
            lSet const & set = setWrapper->reference();
            bc->ApplyBounaryConditionDefaultMethod<rtTypes::equateValue>(set,time, object, fieldName);
          }
        }
      }
    });
}

void BoundaryConditionManager::ApplyInitialConditions( ManagedGroup * domain ) const
{

//  forSubGroups<BoundaryConditionBase>( [&] ( BoundaryConditionBase const * bc
// )-> void
//  {
  for( auto & subGroup : this->GetSubGroups() )
  {
    BoundaryConditionBase const * bc = subGroup.second->group_cast<BoundaryConditionBase const *>();
    if( bc->initialCondition() )
    {

      if( bc->GetElementRegion().empty() )
      {
        string_array objectPath = stringutilities::Tokenize( bc->GetFieldName(), "/");
        localIndex const pathLength = objectPath.size();
        ManagedGroup * currentGroup = domain;
        for( integer a=0 ; a<(pathLength-1) ; ++a )
        {
          currentGroup = currentGroup->GetGroup(objectPath[a]);
        }
        string const fieldName = objectPath[pathLength-1];

        dataRepository::ManagedGroup const * setGroup = currentGroup->GetGroup(dataRepository::keys::sets);
        string_array setNames = bc->GetSetNames();
        for( auto & setName : setNames )
        {
          dataRepository::ViewWrapper<lSet> const * const setWrapper = setGroup->getWrapper<lSet>(setName);
          if( setWrapper != nullptr )
          {
            lSet const & set = setWrapper->reference();
            bc->ApplyBounaryConditionDefaultMethod<rtTypes::equateValue>( set, 0.0, currentGroup, fieldName );
          }
        }
      }
      else
      {
        ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
//        ConstitutiveManager::constitutiveMaps const & constitutiveMaps =
// constitutiveManager->GetMaps(0);
        typename ManagedGroup::subGroupMap::LookupMapType const & constitutiveIndexLookup = constitutiveManager->GetSubGroups().keys();

        lSet targetSet;



        string_array targetPath = stringutilities::Tokenize( bc->GetFieldName(), "/");
        localIndex const targetPathLength = targetPath.size();
        ManagedGroup * targetGroup = domain;
        for( integer a=0 ; a<(targetPathLength-1) ; ++a )
        {
          targetGroup = targetGroup->GetGroup(targetPath[a]);
        }
        string const materialName = targetPath[targetPathLength-3];
        string const fieldName = targetPath[targetPathLength-1];


        FiniteElementManager const * numericalMethodManager = domain->getParent()->GetGroup<FiniteElementManager>(keys::finiteElementManager);
        FiniteElementSpaceManager const * feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);

        // Get element Region
        string const elementRegionName = bc->GetElementRegion();
        ManagedGroup * ElementRegionManager = ManagedGroup::group_cast<DomainPartition*>(domain)->getMeshBody(0)->getMeshLevel(0)->getElemManager();
        ManagedGroup * ElementRegions = ElementRegionManager->GetGroup(keys::elementRegions);
        ElementRegion * elementRegion = ElementRegions->GetGroup<ElementRegion>(elementRegionName);
        // ManagedGroup * elementSubRegions =
        // elementRegion->GetGroup(dataRepository::keys::cellBlockSubRegions);



        auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);
        FiniteElementSpace const * feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);

        string_array setNames = bc->GetSetNames();

        for( auto & subRegionIter : elementRegion->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups() )
        {
          CellBlockSubRegion * subRegion = subRegionIter.second->group_cast<CellBlockSubRegion *>();
//        elementRegion->forCellBlocks( [&] ( CellBlockSubRegion * subRegion )
// -> void
//        {
          auto const & constitutiveMap = subRegion->getReference< std::pair< Array2dT<localIndex>,Array2dT<localIndex> > >(CellBlockSubRegion::viewKeyStruct::constitutiveMapString);
          ManagedGroup const * sets = subRegion->GetGroup(keys::sets);

          for( auto & setName : setNames )
          {
            dataRepository::ViewWrapper<lSet> const * const setWrapper = sets->getWrapper<lSet>(setName);
            if( setWrapper != nullptr )
            {
              lSet const & set = setWrapper->reference();
              localIndex materialIndex = constitutiveIndexLookup.at(materialName);
              for( auto const & k : set )
              {
                if( constitutiveMap.first(k,0) == materialIndex )
                {
                  for( auto q=0 ; q < feSpace->m_finiteElement->n_quadrature_points() ; ++q )
                  {
                    targetSet.insert(constitutiveMap.second(k,q));
                  }
                }
              }
            }
          }

        }

        // ManagedGroup const * elementSets =
        // elementRegion->GetGroup(dataRepository::keys::sets);


        bc->ApplyBounaryConditionDefaultMethod<rtTypes::equateValue>( targetSet, 0.0, targetGroup, fieldName );
      }
    }
  }
}

} /* namespace geosx */
