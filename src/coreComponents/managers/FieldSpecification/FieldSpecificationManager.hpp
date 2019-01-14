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

/**
 * @file BoundaryConditionManager.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#include "managers/FieldSpecification/FieldSpecificationBase.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const boundaryConditionManager( "BoundaryConditionManager" );
}
}

/**
 * @class FieldSpecificationManager
 * This class contains the field objects and provides an interface for administering
 * to specify them The class is a singleton.
 */
class FieldSpecificationManager : public dataRepository::ManagedGroup
{
public:


  /**
   * @brief Singleton getter returns a pointer to the Singleton instance of
   *        BoundaryConditionManager.
   * @return a pointer to the singleton FieldSpecificationManager
   */
  static FieldSpecificationManager * get();

  /**
   * @brief create a new FieldSpecificationBase object as a child of this group.
   * @param childKey the catalog key of the new FieldSpecificationBase derived type to create
   * @param childName the name of the new FieldSpecificationBase object in the repository
   */
  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;


  /**
   * @brief Function to apply a value directly to a field variable.
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the value will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param domain The DomainPartition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the string that specifies the
   *                  path from the actual field specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified, and that will be
   *                  sufficient to proceed with the application of the.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the field specification.
   *
   * Interface function that applies a boundary condition directly to a field variable by looping
   * through all FieldSpecificationBase objects present in the FieldSpecificationManager. Searches
   * for the string specified in fieldPath as a substring in the objectPath specified in the input
   * file, checks if fieldName is equal to the fieldName specified in the input file, and check if
   * the time parameter falls within the beginTime and endTime of the FieldSpecificationBase object,
   * and calls FieldSpecificationBase::ApplyFieldValue().
   *
   */
  void ApplyFieldValue( real64 const time,
                                      dataRepository::ManagedGroup * domain,
                                      string const & fieldPath,
                                      string const & fieldName ) const
  {
	  ApplyFieldValue( time, domain, fieldPath, fieldName,
                                   [&]( FieldSpecificationBase const * const,
                                        set<localIndex> const & ){} );
  }

  /**
   * @brief Function to apply a value directly to a field variable and applies a lambda
   *        for any post operations that are needed.
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the field will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param domain The DomainPartition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the string that specifies the
   *                  path from the actual field specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified , and that will be
   *                  sufficient to proceed with the application of the value.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the field specification.
   * @param lambda A lambda function that defines any operations that should be performed
   *               after application of value to the field.
   *
   * Interface function that applies a value directly to a field by looping
   * through all FieldSpecificationBase objects present in the FieldSpecificationManager. Searches
   * for the string specified in fieldPath as a substring in the objectPath specified in the input
   * file, checks if fieldName is equal to the fieldName specified in the input file, and check if
   * the time parameter falls within the beginTime and endTime of the FieldSpecificationBase object,
   * and calls FieldSpecificationBase::ApplyFieldValue(), and calls the lambda function
   * to apply any operations required for completing the application of the value to the field in addition to
   * setting the target field.
   */
  template< typename LAMBDA >
  void ApplyFieldValue( real64 const time,
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
   * @brief This function is the main driver for the field applications
   * @tparam LAMBDA The type of the lambda function
   * @param time The time at which the field will be evaluated. For instance if the
   *             field is a time dependent function, this is the evaluation time.
   * @param domain The DomainPartition object.
   * @param fieldPath The path to the object that contains the variable described in fieldName. This
   *                  path need not be the complete path, but rather the check that is performed is
   *                  that the fieldPath specified is contained in the string that specifies the
   *                  path from the actual field specification. In other words, the fieldPath variable
   *                  can be a substring of the path specified, and that will be
   *                  sufficient to proceed with the application of the value to the field.
   * @param fieldName The name of the field/variable that the value will be applied to.
   *                  It may not be necessary that this name is in the data repository, as the user
   *                  supplied lambda may apply whatever it condition it would like. However, this
   *                  name is used for comparing against the value given in the specification.
   * @param lambda A lambda function that defines the application of the field.
   *
   * This function loops through all available fields, checks to see if they
   * should be applied, and applies them. More specifically, this function simply checks
   * values of fieldPath,fieldName, against each FieldSpecificationBase object contained in the
   * FieldSpecificationManager and decides on whether or not to call the user defined lambda.
   */
  template< typename LAMBDA >
  void Apply( real64 const time,
              dataRepository::ManagedGroup * domain,
              string const & fieldPath,
              string const & fieldName,
              LAMBDA && lambda ) const
  {
    GEOSX_MARK_FUNCTION;
    for( auto & subGroup : this->GetSubGroups() )
    {
      FieldSpecificationBase const * fs = subGroup.second->group_cast<FieldSpecificationBase const *>();
      int const isInitialCondition = fs->initialCondition();

      if( ( isInitialCondition && fieldPath=="" ) ||
          ( !isInitialCondition && fs->GetObjectPath().find(fieldPath) != string::npos ) )
      {
        string_array const targetPath = stringutilities::Tokenize( fs->GetObjectPath(), "/" );
        localIndex const targetPathLength = integer_conversion<localIndex>(targetPath.size());
        string const targetName = fs->GetFieldName();

        if( ( isInitialCondition && fieldName=="" ) ||
            ( !isInitialCondition && time >= fs->GetStartTime() && time < fs->GetEndTime() && targetName==fieldName ) )
        {
          MeshLevel * const meshLevel = domain->group_cast<DomainPartition*>()->
                                        getMeshBody( 0 )->getMeshLevel( 0 );

          dataRepository::ManagedGroup * targetGroup = meshLevel;

          string processedPath;
          for( localIndex pathLevel=0 ; pathLevel<targetPathLength ; ++pathLevel )
          {
            targetGroup = targetGroup->GetGroup( targetPath[pathLevel] );
            processedPath += "/" + targetPath[pathLevel];

            GEOS_ERROR_IF( targetGroup == nullptr,
                         "ApplyBoundaryCondition(): Last entry in objectPath ("<<processedPath<<") is not found" );
          }

          dataRepository::ManagedGroup const * setGroup = targetGroup->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
          string_array setNames = fs->GetSetNames();
          for( auto & setName : setNames )
          {
            dataRepository::ViewWrapper<set<localIndex> > const * const setWrapper = setGroup->getWrapper<set<localIndex> >( setName );
            if( setWrapper != nullptr )
            {
              set<localIndex> const & targetSet = setWrapper->reference();
              lambda( fs, setName, targetSet, targetGroup, targetName );
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
  FieldSpecificationManager( string const & name, dataRepository::ManagedGroup * const parent );
  virtual ~FieldSpecificationManager() override;

};

template< typename LAMBDA >
void
FieldSpecificationManager::
ApplyFieldValue( real64 const time,
                               dataRepository::ManagedGroup * domain,
                               string const & fieldPath,
                               string const & fieldName,
                               LAMBDA && lambda ) const
{
  FieldSpecificationBase const * fsBase = nullptr;
  set<localIndex> const * targetSetCopy = nullptr;
  Apply( time, domain, fieldPath, fieldName,
        [&]( FieldSpecificationBase const * const fs,
        string const &,
        set<localIndex> const & targetSet,
        ManagedGroup * const targetGroup,
        string const & targetField )
    {
      fs->ApplyFieldValue<FieldSpecificationEqual>( targetSet, time, targetGroup, targetField );
      lambda( fs, targetSet );
    } );


}


} /* namespace geosx */

#endif /*
          SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
        */
