/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "FieldSpecificationManager.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{

FieldSpecificationManager * FieldSpecificationManager::m_instance = nullptr;

using namespace dataRepository;
using namespace constitutive;
FieldSpecificationManager::FieldSpecificationManager( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );

  GEOSX_ERROR_IF( m_instance != nullptr, "Only one FieldSpecificationManager can exist at a time." );
  m_instance = this;

}

FieldSpecificationManager::~FieldSpecificationManager()
{
  GEOSX_ERROR_IF( m_instance != this, "m_instance != this should not be possible." );
  m_instance = nullptr;
}


FieldSpecificationManager & FieldSpecificationManager::getInstance()
{
  GEOSX_ERROR_IF( m_instance == nullptr,
                  "FieldSpecificationManager has not been constructed, or is already been destructed." );
  return *m_instance;
}

Group * FieldSpecificationManager::createChild( string const & childKey, string const & childName )
{
  std::unique_ptr< FieldSpecificationBase > bc = FieldSpecificationBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup( childName, std::move( bc ) );
}


void FieldSpecificationManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BoundaryConditionBase here
  for( auto & catalogIter: FieldSpecificationBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

void FieldSpecificationManager::validateBoundaryConditions( MeshLevel & mesh ) const
{
  // loop over all the FieldSpecification of the XML file
  this->forSubGroups< FieldSpecificationBase >( [&] ( FieldSpecificationBase const & fs )
  {
    localIndex isFieldNameFound = 0;
    // map from set name to a flag (1 if targetSet not empty, 0 otherwise)
    map< string, localIndex > isTargetSetNonEmpty;
    // map from set name to a flag (1 if targetSet has been created, 0 otherwise)
    map< string, localIndex > isTargetSetCreated;

    // Step 1: collect all the set names in a map (this is made necessary by the "apply" loop pattern

    array1d< string > const & setNames = fs.getSetNames();
    for( localIndex i = 0; i < setNames.size(); ++i )
    {
      isTargetSetNonEmpty[setNames[i]] = 0;
      isTargetSetCreated[setNames[i]] = 0;
    }

    // Step 2: apply the boundary condition

    fs.apply< dataRepository::Group >( mesh,
                                       [&]( FieldSpecificationBase const &,
                                            string const & setName,
                                            SortedArrayView< localIndex const > const & targetSet,
                                            Group & targetGroup,
                                            string const fieldName )
    {
      dataRepository::InputFlags const flag = fs.getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).getInputFlag();

      // 2.a) If we enter this loop, we know that the set has been created
      //      Fracture/fault sets are created later and the "apply" call silently ignores them
      isTargetSetCreated.at( setName ) = 1;

      // 2.b) If the fieldName is registered on this target, we record it
      //      Unfortunately, we need two exceptions:
      //       - FieldSpecification that do not target a field, like Aquifer, Traction, Equilibrium, etc. For these, the check is not
      // necessary (the user cannot mess up)
      //       - Face boundary conditions that target cell-based quantities, like the face BC of the flow solvers
      if( targetGroup.hasWrapper( fieldName ) ||
          flag == InputFlags::FALSE || // no need to check input if the input flag is false (Aquifer, Traction, Equilibrium do not target a
                                       // field)
          targetGroup.getName() == MeshLevel::groupStructKeys::faceManagerString() ) // the field names of the face BCs are not always
                                                                                     // registered on
                                                                                     // the faceManager...
      {
        isFieldNameFound = 1;
      }

      // 2.c) If the target set is not empty, we record it
      //      Fracture/fault sets are sometimes at initialization and the "apply" call silently ignores them
      if( targetSet.size() > 0 )
      {
        isTargetSetNonEmpty.at( setName ) = 1;
      }
    } );

    // Step 3: MPI synchronization

    isFieldNameFound = MpiWrapper::max( isFieldNameFound );

    bool areAllSetsEmpty = true;
    for( std::pair< string const, localIndex > & mapEntry : isTargetSetNonEmpty )
    {
      mapEntry.second = MpiWrapper::max( mapEntry.second );
      if( mapEntry.second == 1 )
      {
        areAllSetsEmpty = false;
      }
    }
    bool areAllSetsMissing = true;
    for( std::pair< string const, localIndex > & mapEntry : isTargetSetCreated )
    {
      mapEntry.second = MpiWrapper::max( mapEntry.second );
      if( mapEntry.second == 1 )
      {
        areAllSetsMissing = false;
      }
    }

    // Step 4: issue an error or a warning if the field was not found

    // if all sets are missing, we stop the simulation.
    if( areAllSetsMissing )
    {
      // still loop over the map to get the set name
      for( auto const & mapEntry : isTargetSetCreated )
      {
        GEOSX_THROW( GEOSX_FMT( "\n{}: there is no set named `{}` under the {} `{}`, check the XML input\n",
                                fs.getName(), mapEntry.first, FieldSpecificationBase::viewKeyStruct::objectPathString(), fs.getObjectPath() ),
                     InputError );
      }
    }

    // if a target set is empty, we issue a warning
    // ideally we would just stop the simulation, but the SurfaceGenerator relies on this behavior
    for( auto const & mapEntry : isTargetSetNonEmpty )
    {
      GEOSX_LOG_RANK_0_IF( mapEntry.second == 0,
                           GEOSX_FMT( "\nWarning!"
                                      "\n{}: this FieldSpecification targets (an) empty set(s)"
                                      "\nIf the simulation does not involve the SurfaceGenerator, check the content of the set `{}` in `{}`."
                                      "\n  -If `{}` is in ElementRegions, the set(s) must contain at least an element"
                                      "\n  -If `{}` is faceManager, the set(s) must contain at least a face"
                                      "\n  -If `{}` is nodeManager, the set(s) must contain at least a node"
                                      "\nIf the set of elements/faces is created using a Box in <Geometry>, make sure the nodes of the elements/faces are fully inside the Box\n",
                                      fs.getName(), mapEntry.first,
                                      fs.getObjectPath(), fs.getObjectPath(), fs.getObjectPath(), fs.getObjectPath() ) );
    }

    char const fieldNameNotFoundMessage[] =
      "\n{}: there is no {} named `{}` under the {} `{}`, check the XML input\n";

    // if the field name was not found and the sets are empty, we issue a warning (may be on a fracture region)
    GEOSX_LOG_RANK_0_IF( isFieldNameFound == 0 && areAllSetsEmpty,
                         GEOSX_FMT( fieldNameNotFoundMessage,
                                    fs.getName(), FieldSpecificationBase::viewKeyStruct::fieldNameString(),
                                    fs.getFieldName(), FieldSpecificationBase::viewKeyStruct::objectPathString(), fs.getObjectPath(), fs.getFieldName() ) );
    // if the field name was not found and some sets are not empty, the user misspelled the field name
    GEOSX_THROW_IF( isFieldNameFound == 0 && !areAllSetsEmpty,
                    GEOSX_FMT( fieldNameNotFoundMessage,
                               fs.getName(), FieldSpecificationBase::viewKeyStruct::fieldNameString(),
                               fs.getFieldName(), FieldSpecificationBase::viewKeyStruct::objectPathString(), fs.getObjectPath(), fs.getFieldName() ),
                    InputError );
  } );
}

void FieldSpecificationManager::applyInitialConditions( MeshLevel & mesh ) const
{
  this->forSubGroups< FieldSpecificationBase >( [&] ( FieldSpecificationBase const & fs )
  {
    if( fs.initialCondition() )
    {
      fs.apply< dataRepository::Group >( mesh,
                                         [&]( FieldSpecificationBase const & bc,
                                              string const &,
                                              SortedArrayView< localIndex const > const & targetSet,
                                              Group & targetGroup,
                                              string const fieldName )
      {
        bc.applyFieldValue< FieldSpecificationEqual >( targetSet, 0.0, targetGroup, fieldName );
      } );
    }
  } );
}

} /* namespace geosx */
