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
 * @file BoundaryConditionManager.hpp
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

/**
 * @class BoundaryConditionManager
 * This class contains the boundary condition objects and provides an interface for administering
 * the boundary conditions. The class is a singleton.
 */
class BoundaryConditionManager : public dataRepository::ManagedGroup
{
public:


  /**
   * @brief Singleton getter returns a pointer to the Singleton instance of
   *        BoundaryConditionManager.
   * @return a pointer to the singleton BoundaryConditionManager
   */
  static BoundaryConditionManager * get();

  /**
   * @brief create a new BoundaryConditionBase object as a child of this group.
   * @param childKey the catalog key of the new BoundaryConditionBase derived type to create
   * @param childName the name of the new BoundaryConditionBase object in the repository
   */
  virtual void CreateChild( string const & childKey, string const & childName ) override;


  /**
   * @brief Function to apply a boundary condition directly to a field variable.
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the boundary condition will be evaluated. For instance if the
   *             boundary condition is a time dependent function, this is the evaluation time.
   * @param domain The DomainParition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the BC string that specifies the
   *                  path from the actual BC specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified in the BC, and that will be
   *                  sufficient to proceed with the application of the BC.
   * @param fieldName The name of the field/variable that the boundary condition will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the BC specification.
   *
   * Interface function that applies a boundary condition directly to a field variable by looping
   * through all BoundaryConditionBase objects present in the BoundaryConditionManager. Searches
   * for the string specified in fieldPath as a substring in the objectPath specified in the input
   * file, checks if fieldName is equal to the fieldName specified in the input file, and check if
   * the time parameter falls within the beginTime and endTime of the BoundaryConditionBase object,
   * and calls BoundaryConditionBase::ApplyBoundaryConditionToField().
   *
   */
  void ApplyBoundaryConditionToField( real64 const time,
                                      dataRepository::ManagedGroup * domain,
                                      string const & fieldPath,
                                      string const & fieldName ) const
  {
    ApplyBoundaryConditionToField( time, domain, fieldPath, fieldName,
                                   [&]( BoundaryConditionBase const * const,
                                        set<localIndex> const & ){} );
  }

  /**
   * @brief Function to apply a boundary condition directly to a field variable and applies a lambda
   *        for any post BC operations that are needed.
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the boundary condition will be evaluated. For instance if the
   *             boundary condition is a time dependent function, this is the evaluation time.
   * @param domain The DomainParition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the BC string that specifies the
   *                  path from the actual BC specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified in the BC, and that will be
   *                  sufficient to proceed with the application of the BC.
   * @param fieldName The name of the field/variable that the boundary condition will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the BC specification.
   * @param lambda A lambda function that defines any operations that should be performed
   *               after application of the boundary condition.
   *
   * Interface function that applies a boundary condition directly to a field variable by looping
   * through all BoundaryConditionBase objects present in the BoundaryConditionManager. Searches
   * for the string specified in fieldPath as a substring in the objectPath specified in the input
   * file, checks if fieldName is equal to the fieldName specified in the input file, and check if
   * the time parameter falls within the beginTime and endTime of the BoundaryConditionBase object,
   * and calls BoundaryConditionBase::ApplyBoundaryConditionToField(), and calls the lambda function
   * to apply any operations required for completing the application of the BC in addition to
   * setting the target field.
   */
  template< typename LAMBDA >
  void ApplyBoundaryConditionToField( real64 const time,
                                      dataRepository::ManagedGroup * domain,
                                      string const & fieldPath,
                                      string const & fieldName,
                                      LAMBDA && lambda ) const;

  /**
   * @brief function to apply initial conditions
   * @param domain the DomainParition object
   */
  void ApplyInitialConditions( dataRepository::ManagedGroup * domain ) const;


  /**
   * @brief This function is the main driver for the application of boundary conditions.
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the boundary condition will be evaluated. For instance if the
   *             boundary condition is a time dependent function, this is the evaluation time.
   * @param domain The DomainParition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the BC string that specifies the
   *                  path from the actual BC specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified in the BC, and that will be
   *                  sufficient to proceed with the application of the BC.
   * @param fieldName The name of the field/variable that the boundary condition will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the BC specification.
   * @param lambda A lambda function that defines the application of the boundary condition.
   *
   * This function loops through all available boundary conditions, checks to see if the BC/IC
   * should be applied, and applies the condition. More specifically, this function simply checks
   * values of fieldPath,fieldName, against each BoundaryConditionBase object contained in the
   * BoundaryConditionManager and decides on whether or not to call the user defined lambda.
   */
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

      if( ( isInitialCondition && fieldPath=="" ) ||
          ( !isInitialCondition && bc->GetObjectPath().find(fieldPath) != string::npos ) )
      {
        string_array const targetPath = stringutilities::Tokenize( bc->GetObjectPath(), "/" );
//        std::cout<<"objectPath = "<<bc->GetObjectPath()<<std::endl;
        localIndex const targetPathLength = targetPath.size();
        string const targetName = bc->GetFieldName();
//        std::cout<<"targetName = "<<targetName<<std::endl;

        if( ( isInitialCondition && fieldName=="" ) ||
            ( !isInitialCondition && time >= bc->GetStartTime() && time < bc->GetEndTime() && targetName==fieldName ) )
        {
//          std::cout<<bc->getName()<<std::endl;

          MeshLevel * const meshLevel = domain->group_cast<DomainPartition*>()->
                                        getMeshBody( 0 )->getMeshLevel( 0 );

          dataRepository::ManagedGroup * targetGroup = meshLevel;

          string processedPath;
          for( localIndex pathLevel=0 ; pathLevel<targetPathLength ; ++pathLevel )
          {
//            std::cout<<targetPath[pathLevel]<<std::endl;

            targetGroup = targetGroup->GetGroup( targetPath[pathLevel] );
            processedPath += "/" + targetPath[pathLevel];
//            std::cout<<"processedPath="<<processedPath<<std::endl;

            GEOS_ERROR_IF( targetGroup == nullptr,
                         "ApplyBoundaryCondition(): Last entry in objectPath ("<<processedPath<<") is not found" );
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

private:
  /**
   * @brief private constructor for the singleton BoundaryConditionManager.
   * @param name The name of the BoundaryConditionManager in the data repository.
   * @param parent The parent of BoundaryConditionManager in the data repository.
   */
  BoundaryConditionManager( string const & name, dataRepository::ManagedGroup * const parent );
  virtual ~BoundaryConditionManager() override;

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
