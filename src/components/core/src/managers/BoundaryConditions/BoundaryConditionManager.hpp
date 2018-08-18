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
 * BoundaryConditionManager.hpp
 *
 *  Created on: May 26, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#include "common/DataTypes.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "managers/DomainPartition.hpp"
#include "BoundaryConditionBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const boundaryConditionMananger( "BoundaryConditionMananger" );
}
}


class BoundaryConditionManager : public dataRepository::ManagedGroup
{
public:
  BoundaryConditionManager( string const & name, dataRepository::ManagedGroup * const parent );
  virtual ~BoundaryConditionManager() override;

  static BoundaryConditionManager * get();

  virtual void CreateChild( string const & childKey, string const & childName ) override;


  void ApplyBoundaryConditionToField( real64 const time,
                                      dataRepository::ManagedGroup * domain,
                                      string const & fieldPath,
                                      string const & fieldName ) const
  {
    ApplyBoundaryConditionToField( time, domain, fieldPath, fieldName,
                                   [&]( BoundaryConditionBase const * const,
                                        set<localIndex> const & ){} );
  }

  template< typename LAMBDA >
  void ApplyBoundaryConditionToField( real64 const time,
                                      dataRepository::ManagedGroup * domain,
                                      string const & fieldPath,
                                      string const & fieldName,
                                      LAMBDA && lambda ) const;

  void ApplyInitialConditions( dataRepository::ManagedGroup * domain ) const;



  template< typename LAMBDA >
  void ApplyBoundaryCondition( real64 const time,
                               dataRepository::ManagedGroup * domain,
                               string const & fieldPath,
                               string const & fieldName,
                               LAMBDA && lambda ) const
  {
    for( auto & subGroup : this->GetSubGroups() )
    {
      BoundaryConditionBase const * bc = subGroup.second->group_cast<BoundaryConditionBase const *>();
      int const isInitialCondition = bc->initialCondition();

      if( ( isInitialCondition && fieldPath=="") || ( bc->GetObjectPath() == fieldPath ) )
      {
        string_array const targetPath = stringutilities::Tokenize( bc->GetObjectPath(), "/" );
        std::cout<<"objectPath = "<<bc->GetObjectPath()<<std::endl;
        localIndex const targetPathLength = targetPath.size();
        string const targetName = bc->GetFieldName();
        std::cout<<"targetName = "<<targetName<<std::endl;

        if( ( isInitialCondition && fieldName=="" ) ||
            ( time >= bc->GetStartTime() && time < bc->GetEndTime() && targetName==fieldName ) )
        {

          MeshLevel * const meshLevel = domain->group_cast<DomainPartition*>()->
                                        getMeshBody( 0 )->getMeshLevel( 0 );

          dataRepository::ManagedGroup * targetGroup = meshLevel;

          string processedPath;
          for( localIndex pathLevel=0 ; pathLevel<targetPathLength ; ++pathLevel )
          {
            std::cout<<targetPath[pathLevel]<<std::endl;

            targetGroup = targetGroup->GetGroup( targetPath[pathLevel] );
            processedPath += "/" + targetPath[pathLevel];
            std::cout<<"processedPath="<<processedPath<<std::endl;

            GEOS_ASSERT( targetGroup != nullptr,
                         "ApplyBoundaryCondition(): Last entry in objectPath ("<<processedPath<<") is not found" )
          }

          dataRepository::ManagedGroup const * setGroup = targetGroup->GetGroup( dataRepository::keys::sets );
          string_array setNames = bc->GetSetNames();
          for( auto & setName : setNames )
          {
            dataRepository::ViewWrapper<set<localIndex> > const * const setWrapper = setGroup->getWrapper<set<localIndex> >( setName );
            if( setWrapper != nullptr )
            {
              set<localIndex> const & targetSet = setWrapper->reference();
              lambda( bc, setName, targetSet, targetGroup, targetName );
            }
          }
        }
      }
    }
  }
};

template< typename LAMBDA >
void
BoundaryConditionManager::
ApplyBoundaryConditionToField( real64 const time,
                               dataRepository::ManagedGroup * domain,
                               string const & fieldPath,
                               string const & fieldName,
                               LAMBDA && lambda ) const
{
  BoundaryConditionBase const * bcBase = nullptr;
  set<localIndex> const * targetSetCopy = nullptr;
  ApplyBoundaryCondition( time, domain, fieldPath, fieldName,
                          [&]( BoundaryConditionBase const * const bc,
                               string const &,
                               set<localIndex> const & targetSet,
                               ManagedGroup * const targetGroup,
                               string const & targetField )
    {
      bc->ApplyBoundaryConditionToField<BcEqual>( targetSet, time, targetGroup, targetField );
      lambda( bc, targetSet );
    } );


}


} /* namespace geosx */

#endif /*
          SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
        */
