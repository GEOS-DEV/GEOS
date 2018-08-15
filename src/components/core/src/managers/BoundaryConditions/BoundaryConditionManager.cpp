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

#include "managers/NumericalMethodsManager.hpp"
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
          dataRepository::ViewWrapper<set<localIndex>> const * const setWrapper = sets->getWrapper<set<localIndex>>(setName);
          if( setWrapper != nullptr )
          {
            set<localIndex> const & set = setWrapper->reference();
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
      string_array const objectPath = stringutilities::Tokenize( bc->GetObjectPath(), "/");
      localIndex const targetPathLength = objectPath.size();
      ManagedGroup * targetGroup = domain;

      if( objectPath[0]=="ElementRegion" )
      {
        ElementRegionManager * const
        elementRegionManager = domain->group_cast<DomainPartition*>()->
                               getMeshBody(0)->getMeshLevel(0)->getElemManager();
        ElementRegion * const elementRegion = elementRegionManager->GetRegion(objectPath[1]);

        CellBlockSubRegion * specifiedSubRegion = nullptr;
        if( objectPath.size() > 2 )
        {
          specifiedSubRegion = elementRegion->GetSubRegion(objectPath[2]);
        }

        if( bc->GetConstitutivePath().empty() )
        {

          string fieldName = bc->GetFieldName();

          if( objectPath.size()>3 )
          {
            GEOS_ASSERT( !( !bc->GetFieldName().empty() && !objectPath[3].empty() ) ,
                         "fieldName specified in both fieldName entry ("<<bc->GetFieldName()
                         <<") and objectPath ("<<objectPath[3]<<")");

            GEOS_ASSERT( !( bc->GetFieldName().empty() && objectPath[3].empty() ),
                         "fieldName not specified in either fieldName entry or objectPath");

            if( !objectPath[3].empty() )
            {
              fieldName = objectPath[3];
            }
          }
          else
          {
            GEOS_ASSERT( !bc->GetFieldName().empty(),
                         "fieldName not specified in either fieldName entry or objectPath" );
          }

          for( auto & subRegionIter : elementRegion->GetSubRegions() )
          {
            CellBlockSubRegion * subRegion = subRegionIter.second->group_cast<CellBlockSubRegion *>();

            if( specifiedSubRegion==nullptr || specifiedSubRegion == subRegion)
            {
              targetGroup = subRegion;

              dataRepository::ManagedGroup const * setGroup = targetGroup->GetGroup(dataRepository::keys::sets);
              string_array setNames = bc->GetSetNames();
              for( auto & setName : setNames )
              {
                dataRepository::ViewWrapper<set<localIndex>> const * const setWrapper = setGroup->getWrapper<set<localIndex>>(setName);
                if( setWrapper != nullptr )
                {
                  set<localIndex> const & set = setWrapper->reference();
                  bc->ApplyBounaryConditionDefaultMethod<rtTypes::equateValue>( set, 0.0, targetGroup, fieldName );
                }
              }
            }
          }
        } // if( bc->GetConstitutivePath().empty() )
        else
        {
          string_array const constitutivePath = stringutilities::Tokenize( bc->GetConstitutivePath(), "/");
          ConstitutiveManager * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

          typename ManagedGroup::subGroupMap::LookupMapType const & constitutiveIndexLookup = constitutiveManager->GetSubGroups().keys();

          string const materialName = constitutivePath[0];
          string const paramOrState = constitutivePath[1];
          targetGroup = constitutiveManager->GetGroup<ConstitutiveBase>(materialName)
                                           ->GetGroup(paramOrState);

          string fieldName = bc->GetFieldName();

          if( constitutivePath.size()>2 )
          {
            GEOS_ASSERT( !( !bc->GetFieldName().empty() && !constitutivePath[2].empty() ) ,
                         "fieldName specified in both fieldName entry ("<<bc->GetFieldName()
                         <<") and constitutivePath ("<<constitutivePath[2]<<")");

            GEOS_ASSERT( !( bc->GetFieldName().empty() && constitutivePath[2].empty() ),
                         "fieldName not specified in either fieldName entry or constitutivePath");

            if( !constitutivePath[2].empty() )
            {
              fieldName = constitutivePath[2];
            }
          }
          else
          {
            GEOS_ASSERT( !bc->GetFieldName().empty(),
                         "fieldName not specified in either fieldName entry or constitutivePath" );
          }

          set<localIndex> targetSet;


          NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
          FiniteElementSpaceManager const * feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);


          auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);
          FiniteElementSpace const * feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);

          string_array setNames = bc->GetSetNames();

          for( auto & subRegionIter : elementRegion->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups() )
          {
            CellBlockSubRegion * subRegion = subRegionIter.second->group_cast<CellBlockSubRegion *>();

            if( specifiedSubRegion==nullptr || specifiedSubRegion == subRegion)
            {
              auto const & constitutiveMap = subRegion->getReference< std::pair< array2d<localIndex>,array2d<localIndex> > >(CellBlockSubRegion::viewKeyStruct::constitutiveMapString);
              ManagedGroup const * sets = subRegion->GetGroup(keys::sets);

              for( auto & setName : setNames )
              {
                dataRepository::ViewWrapper<set<localIndex>> const * const setWrapper = sets->getWrapper<set<localIndex>>(setName);
                if( setWrapper != nullptr )
                {
                  set<localIndex> const & set = setWrapper->reference();
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
          }
          bc->ApplyBounaryConditionDefaultMethod<rtTypes::equateValue>( targetSet, 0.0, targetGroup, fieldName );
        }
      } // if( objectPath[0]=="ElementRegion" )
      else
      {

        targetGroup = domain->group_cast<DomainPartition*>()->
                      getMeshBody(0)->getMeshLevel(0)->GetGroup(objectPath[0]);

        string fieldName = bc->GetFieldName();

        if( objectPath.size()>1 )
        {
          GEOS_ASSERT( !( !bc->GetFieldName().empty() && !objectPath[1].empty() ) ,
                       "fieldName specified in both fieldName entry ("<<bc->GetFieldName()
                       <<") and objectPath ("<<objectPath[1]<<")");

          GEOS_ASSERT( !( bc->GetFieldName().empty() && objectPath[1].empty() ),
                       "fieldName not specified in either fieldName entry or objectPath");

          if( !objectPath[1].empty() )
          {
            fieldName = objectPath[1];
          }
        }
        else
        {
          GEOS_ASSERT( !bc->GetFieldName().empty(),
                       "fieldName not specified in either fieldName entry or objectPath" );
        }


        dataRepository::ManagedGroup const * setGroup = targetGroup->GetGroup(dataRepository::keys::sets);
        string_array setNames = bc->GetSetNames();
        for( auto & setName : setNames )
        {
          dataRepository::ViewWrapper<set<localIndex>> const * const setWrapper = setGroup->getWrapper<set<localIndex>>(setName);
          if( setWrapper != nullptr )
          {
            set<localIndex> const & set = setWrapper->reference();
            bc->ApplyBounaryConditionDefaultMethod<rtTypes::equateValue>( set, 0.0, targetGroup, fieldName );
          }
        }
      }

    }
  }
}

} /* namespace geosx */
