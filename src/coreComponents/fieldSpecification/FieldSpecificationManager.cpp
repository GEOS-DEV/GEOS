/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "FieldSpecificationManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshObjectPath.hpp"

namespace geos
{

FieldSpecificationManager * FieldSpecificationManager::m_instance = nullptr;

using namespace dataRepository;
using namespace constitutive;

FieldSpecificationManager::FieldSpecificationManager( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );

  GEOS_ERROR_IF( m_instance != nullptr, "Only one FieldSpecificationManager can exist at a time." );
  m_instance = this;

}

FieldSpecificationManager::~FieldSpecificationManager()
{
  GEOS_ERROR_IF( m_instance != this, "m_instance != this should not be possible." );
  m_instance = nullptr;
}


FieldSpecificationManager & FieldSpecificationManager::getInstance()
{
  GEOS_ERROR_IF( m_instance == nullptr,
                 "FieldSpecificationManager has not been constructed, or is already been destructed." );
  return *m_instance;
}

Group * FieldSpecificationManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< FieldSpecificationBase > bc =
    FieldSpecificationBase::CatalogInterface::factory( childKey, getDataContext(), childName, this );
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

/// @brief alias for the map containing the relationship between the fieldSpefication "setNames" and the/their associated manager(s)
/// Represented as : { { "setName", { ObjectPathType (= A Manager) , localIndex (0,1) }, ... }, ... }
using SetNameToTypesMap = std::map< std::string, std::vector< std::pair< MeshObjectPath::ObjectTypes, localIndex > > >;

/**
 * @brief Iterate through all the manager of meshLevel
 * @tparam TYPE The current type to be tested
 * @tparam NEXT_TYPES The remaining Types
 * @param mesh Holds all the managers.
 * @param targetType The target fieldSpecification object path type
 * @param setName The current setName to be evaluated in all the managers
 * @param setTypesMap The map containing the relation between the setNames and their managers
 */
template< typename TYPE, typename ... NEXT_TYPES >
void forManagerInMeshLevel( MeshLevel const & mesh,
                            MeshObjectPath::ObjectTypes const targetType,
                            string const & setName, SetNameToTypesMap & setTypesMap )
{
  mesh.forSubGroups< TYPE >( [&]( Group const & targetManager )
  {
    TYPE const * manager = dynamic_cast< TYPE const * >( &targetManager );

    if( manager != nullptr )
    {
      if( manager->sets().hasWrapper( setName ))
      {
        auto const & targetSet = manager->getSet( setName );

        if( std::is_same_v< TYPE, NodeManager > &&
            targetType !=  MeshObjectPath::ObjectTypes::nodes &&
            targetSet.size() > 0 )
        {
          setTypesMap[setName].push_back( { MeshObjectPath::ObjectTypes::nodes, 0 } );
        }
        else if( std::is_same_v< TYPE, EdgeManager >  &&
                 targetType !=  MeshObjectPath::ObjectTypes::edges &&
                 targetSet.size() > 0 )
        {
          setTypesMap[setName].push_back( { MeshObjectPath::ObjectTypes::edges, 0 } );
        }
        else if( std::is_same_v< TYPE, FaceManager > &&
                 targetType !=  MeshObjectPath::ObjectTypes::faces &&
                 targetSet.size() > 0 )
        {
          setTypesMap[setName].push_back( { MeshObjectPath::ObjectTypes::faces, 0 } );
        }
      }
    }
  } );

  if constexpr ( sizeof...(NEXT_TYPES) > 0 )
    forManagerInMeshLevel< NEXT_TYPES... >( mesh, targetType, setName, setTypesMap );
}

void FieldSpecificationManager::validateBoundaryConditions( MeshLevel & mesh ) const
{
  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  Group const & meshBodies = domain.getMeshBodies();
  // loop over all the FieldSpecification of the XML file
  this->forSubGroups< FieldSpecificationBase >( [&] ( FieldSpecificationBase const & fs )
  {
    localIndex isFieldNameFound = 0;
    // map from set name to a flag (1 if targetSet has been created, 0 otherwise)
    map< string, localIndex > isTargetSetCreated;
    // map from all the targeted setNames where we store the associated pair [ObjectType, Bool]
    // 0 if the target objectType exists, 1 otherwise
    SetNameToTypesMap setTypesMap;
    // The fs target objectPath type
    MeshObjectPath::ObjectTypes const expectedSetType = fs.getMeshObjectPaths().getObjectType();


    // Step 1: collect all the set names in a map (this is made necessary by the "apply" loop pattern

    string_array const & setNames = fs.getSetNames();
    for( size_t i = 0; i < setNames.size(); ++i )
    {
      isTargetSetCreated[setNames[i]] = 0;
      setTypesMap[setNames[i]].push_back( {expectedSetType, 1 } );
    }

    // We have to make sure that the meshLevel is in the target of the boundary conditions
    // This is important for multi-level simulations, such as high-order wave propagation
    MeshObjectPath const & objectPath = fs.getMeshObjectPaths();
    if( !objectPath.containsMeshLevel( mesh ) )
    {
      return;
    }

    // Step 2: apply the boundary condition

    fs.apply< Group >( mesh,
                       [&]( FieldSpecificationBase const &,
                            string const & setName,
                            SortedArrayView< localIndex const > const & targetSet,
                            Group & targetGroup,
                            string const fieldName )
    {
      InputFlags const flag = fs.getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).getInputFlag();
      isTargetSetCreated.at( setName ) = 1;
      if( targetGroup.hasWrapper( fieldName ) ||flag == InputFlags::FALSE ||
          targetGroup.getName() == MeshLevel::groupStructKeys::faceManagerString() )     // the field names of the face BCs are not always
                                                                                         // registered on
                                                                                         // the faceManager...
      {
        isFieldNameFound = 1;
      }

      if( targetSet.size() > 0 )
        std::get< 1 >( setTypesMap[setName].front()) = 0;
    } );

    // Step 3: MPI synchronization
    isFieldNameFound = MpiWrapper::max( isFieldNameFound );

    string_array setNamesToErase;
    for( auto & [setName, typesValue] : setTypesMap )
    {
      auto & typeValue =  std::get< 1 >( typesValue.front() );
      typeValue =  MpiWrapper::min( typeValue );
      if( typeValue == 0 )
      {
        setNamesToErase.push_back( setName );
      }
    }

    for( auto const & setName : setNamesToErase )
    {
      setTypesMap.erase( setName );
    }

    for( std::pair< string const, localIndex > & mapEntry : isTargetSetCreated )
    {
      mapEntry.second = MpiWrapper::max( mapEntry.second );
    }
    bool areAllSetsMissing = true;
    for( std::pair< string const, localIndex > & mapEntry : isTargetSetCreated )
    {
      if( mapEntry.second == 1 ) // target set has been created
      {
        areAllSetsMissing = false;
      }
    }

    // Step 4: issue an error or a warning if the field was not found

    // if all sets are missing, we stop the simulation.
    if( areAllSetsMissing )
    {
      // loop again over the map to collect the set names
      string_array missingSetNames;
      for( auto const & mapEntry : isTargetSetCreated )
      {
        missingSetNames.emplace_back( mapEntry.first );
      }

      std::set< string > registeredSets;
      Group const & meshBody = domain.getMeshBodies();
      objectPath.forObjectsInPath( meshBody,
                                   [&]( Group const & targetGroup )
      {
        string_array availableRegions;
        ObjectManagerBase const & targetObject = dynamic_cast< ObjectManagerBase const & >(targetGroup);
        targetObject.sets().forWrappers( [&] ( dataRepository::WrapperBase const & wrapper )
        {
          registeredSets.insert( wrapper.getName());
        } );
      } );

      std::ostringstream errorMessageBuilder;
      errorMessageBuilder << GEOS_FMT( "\n{}: there are no set(s) named `{}` under the {} `{}`.\n",
                                       fs.getWrapperDataContext( FieldSpecificationBase::viewKeyStruct::objectPathString() ),
                                       fmt::join( missingSetNames, ", " ),
                                       FieldSpecificationBase::viewKeyStruct::objectPathString(), fs.getObjectPath() );
      errorMessageBuilder << ( !registeredSets.empty() ?
                               GEOS_FMT( "Available set(s) are: {}",
                                         stringutilities::join( registeredSets, ", " ) ) :
                               GEOS_FMT( "No set are available for the targeted `{}`",
                                         FieldSpecificationBase::viewKeyStruct::objectPathString() ) );

      GEOS_THROW( errorMessageBuilder.str(), InputError );
    }

    if( !setTypesMap.empty() && MpiWrapper::commRank() == 0 )
    {
      std::ostringstream message;
      message << GEOS_FMT( "{}: this FieldSpecification targets (an) empty set(s).\n",
                           fs.getDataContext() );

      for( auto const & [setName, objectTypesPair] : setTypesMap )
      {
        forManagerInMeshLevel< FaceManager,
                               EdgeManager,
                               NodeManager >( mesh, expectedSetType, setName, setTypesMap );

        string_array capturedTypes;
        for( auto const & pair : objectTypesPair )
        {
          if( pair.first != expectedSetType )
          {
            capturedTypes.push_back( EnumStrings< MeshObjectPath::ObjectTypes >::toString( pair.first ) );
          }
        }

        if( !capturedTypes.empty())
        {
          message << GEOS_FMT( "Set '{}':\n"
                               "  - Does not capture: {}\n", setName, fs.getObjectPath());
          message << GEOS_FMT( "  - Instead, captures: {}\n", stringutilities::join( capturedTypes, ", " ));
        }
        else
        {
          message << GEOS_FMT( "Set '{}' does not capture anything in the mesh ", setName );
        }
      }

      Wrapper< integer > const & wrapper = fs.getWrapper< integer >(
        FieldSpecificationBase::viewKeyStruct::emptySetErrorModeString());
      switch( wrapper.getDefaultValue() )
      {
        case  FieldSpecificationBase::SetErrorMode::Silent:
          break;
        case  FieldSpecificationBase::SetErrorMode::Error:
          GEOS_THROW( message.str(), InputError );
          break;
        case  FieldSpecificationBase::SetErrorMode::SurfaceGeneratorWarning:
          message << GEOS_FMT( "As the simulation includes a SurfaceGenerator, the set may be modified later",
                               fs.getDataContext() );
          GEOS_WARNING( message.str() );
          break;
      }
    }

    MpiWrapper::barrier( MPI_COMM_GEOS );

    if( isFieldNameFound == 0 )
    {
      std::ostringstream errorMessageBuilder;
      errorMessageBuilder << GEOS_FMT( "\n{}: there are no field named `{}` under the region `{}`.\n",
                                       fs.getWrapperDataContext( FieldSpecificationBase::viewKeyStruct::fieldNameString() ),
                                       fs.getFieldName(), fs.getObjectPath() );

      std::set< string > registeredFields;
      objectPath.forObjectsInPath( meshBodies,
                                   [&]( Group const & targetGroup )
      {
        targetGroup.forWrappers( [&]( dataRepository::WrapperBase const & wrapper ) {
          if( wrapper.getPlotLevel() != dataRepository::PlotLevel::NOPLOT )
            registeredFields.insert( wrapper.getName());
        } );
      } );

      if( !registeredFields.empty())
        errorMessageBuilder << ( !registeredFields.empty() ?
                                 GEOS_FMT( "Available fields in {} are:\n{{ {} }}", fs.getObjectPath(),
                                           stringutilities::join( registeredFields, ", " )) :
                                 GEOS_FMT( "No available field in {}.",
                                           fs.getObjectPath() ) );

      GEOS_THROW( errorMessageBuilder.str(), InputError );
    }
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
                                              SortedArrayView< localIndex const > const & targetObject,
                                              Group & targetGroup,
                                              string const fieldName )
      {
        bc.applyFieldValue< FieldSpecificationEqual >( targetObject, 0.0, targetGroup, fieldName );
      } );
    }
  } );
}

}   /* namespace geos */
