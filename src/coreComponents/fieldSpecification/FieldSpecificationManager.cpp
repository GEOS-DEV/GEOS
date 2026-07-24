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
#include "FieldSpecificationABC.hpp"
#include "FieldSpecificationFactory.hpp"
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
  std::unique_ptr< FieldSpecificationABC > bc =
    FieldSpecificationABC::CatalogInterface::factory( childKey, getDataContext(), childName, this );
  return &this->registerGroup( childName, std::move( bc ) );
}


void FieldSpecificationManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BoundaryConditionBase here
  for( auto & catalogIter: FieldSpecificationABC::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

void FieldSpecificationManager::postInputInitialization()
{
  using ProcessorRegistry = FieldSpecificationProcessorRegistry;

  // as the subgroup list can change during expansion, we need an immutable list
  stdVector< FieldSpecificationABC const * > fieldSpecifications;
  this->forSubGroups< FieldSpecificationABC >( [&] ( FieldSpecificationABC const & fs )
  {
    fieldSpecifications.push_back(&fs);
  } );

  for (FieldSpecificationABC const * fs : fieldSpecifications)
  {
    GEOS_LOG( GEOS_FMT( "----- Processing {}", fs->getName()));
    auto const & processors = ProcessorRegistry::getProcessors();
    auto it = processors.find( fs->getCatalogName());
    if( it != processors.end())
    {
      GEOS_LOG( GEOS_FMT( "      - Found Processor {}", it->first ));
      ProcessorRegistry::ProcessorBase const & processor = *it->second;
      processor.expandFieldSpecification( *fs, *this );
    } else {
      GEOS_LOG( GEOS_FMT( "      - NO PROCESSOR !! {}", it->first ));
    }
  }
}

void FieldSpecificationManager::validateBoundaryConditions( MeshLevel & mesh )
{
  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  Group const & meshBodies = domain.getMeshBodies();
  // loop over all the FieldSpecification of the XML file
  this->forSubGroups< FieldSpecification >( [&] ( FieldSpecification & fs )
  {
    localIndex isFieldNameFound = 0;
    // map from set name to a flag (1 if targetSet has been created, 0 otherwise)
    stdMap< string, localIndex > isTargetSetCreated;
    MeshObjectPath::SetNameToTypesMap setTypesMap;
    // The fieldSpecification target objectType
    MeshObjectPath::ObjectTypes const targetElementType = fs.getMeshObjectPaths().getObjectType();


    // Step 1: collect all the set names in a map (this is made necessary by the "apply" loop pattern

    string_array const & setNames = fs.getSetNames();
    for( size_t i = 0; i < setNames.size(); ++i )
    {
      isTargetSetCreated.get_inserted( setNames[i] ) = 0;
      setTypesMap.get_inserted( setNames[i] ).get_inserted( targetElementType ) = 1;
    }

    // We have to make sure that the meshLevel is in the target of the boundary conditions
    // This is important for multi-level simulations, such as high-order wave propagation
    MeshObjectPath const & objectPath = fs.getMeshObjectPaths();
    if( !objectPath.containsMeshLevel( mesh ) )
    {
      return;
    }

    // Step 2: apply the boundary condition

    FieldSpecificationImpl::apply< Group >( fs,
                                            mesh,
                                            [&]( FieldSpecification const &,
                                                 string const & setName,
                                                 SortedArrayView< localIndex const > const & targetSet,
                                                 Group & targetGroup,
                                                 string const fieldName )
    {
      InputFlags const flag = fs.getWrapper< string >( FieldSpecification::viewKeyStruct::fieldNameString() ).getInputFlag();

      // 2.a) If we enter this loop, we know that the set has been created
      //      Fracture/fault sets are created later and the "apply" call silently ignores them
      isTargetSetCreated.at( setName ) = 1;

      // 2.b) If the fieldName is registered on this target, we record it
      //      Unfortunately, we need two exceptions:
      //       - FieldSpecification that do not target a field, like Aquifer, Traction, Equilibrium, etc. For these, the check is not
      // necessary (the user cannot mess up)
      //       - Face boundary conditions that target cell-based quantities, like the face BC of the flow solvers
      if( targetGroup.hasWrapper( fieldName ) || flag == InputFlags::FALSE ||
          targetGroup.getName() == MeshLevel::groupStructKeys::faceManagerString() ) // the field names of the face BCs are not always
                                                                                     // registered on
                                                                                     // the faceManager...
      {
        isFieldNameFound = 1;
      }

      if( targetGroup.hasWrapper( fieldName ) )
      {
        localIndex numElem = targetGroup.getWrapperBase( fieldName ).numArrayComp();
        fs.validateNumArrayComp( numElem );
      }

      if( targetSet.size() > 0 )
      {
        // Use the single set name provided by the callback (setName), not the outer array (setNames)
        setTypesMap.get_inserted( setName ).get_inserted( targetElementType ) = 0;
      }
    } );

    // Step 3: MPI synchronization
    isFieldNameFound = MpiWrapper::max( isFieldNameFound );

    string_array setNamesToErase;
    for( auto & [setName, typesAssociated] : setTypesMap )
    {
      localIndex & elementTypeValue = typesAssociated[targetElementType];
      elementTypeValue =  MpiWrapper::min( elementTypeValue );
      if( elementTypeValue == 0 )
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

      objectPath.forObjectsInPath( meshBodies,
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
                                       fs.getWrapperDataContext( FieldSpecification::viewKeyStruct::objectPathString() ),
                                       fmt::join( missingSetNames, ", " ),
                                       FieldSpecification::viewKeyStruct::objectPathString(), fs.getObjectPath() );
      errorMessageBuilder << ( !registeredSets.empty() ?
                               GEOS_FMT( "Available set(s) are: {}",
                                         stringutilities::join( registeredSets, ", " ) ) :
                               GEOS_FMT( "No set are available for the targeted `{}`",
                                         FieldSpecification::viewKeyStruct::objectPathString() ) );

      GEOS_THROW( errorMessageBuilder.str(), InputError, getDataContext() );
    }

    if( !setTypesMap.empty() && MpiWrapper::commRank() == 0 )
    {
      std::ostringstream message;
      message << GEOS_FMT( "{}: this FieldSpecification targets (an) empty set(s).\n",
                           fs.getDataContext() );

      for( auto const & [setName, elementsType] : setTypesMap )
      {
        objectPath.forManagersForSetName< FaceManager,
                                          EdgeManager,
                                          NodeManager >( mesh, setName,
                                                         [&setTypesMap]( MeshObjectPath::ObjectTypes managerType,
                                                                         string const & capturedSetName ){
          setTypesMap.get_inserted( capturedSetName ).get_inserted( managerType ) = 0;
        } );
        string_array capturedTypes;
        for( auto const & [element, isExisting] : elementsType )
        {
          if( element != targetElementType )
          {
            capturedTypes.push_back( EnumStrings< MeshObjectPath::ObjectTypes >::toString( element ) );
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
      auto errMode = fs.getReference< FieldSpecification::SetErrorMode >(
        FieldSpecification::viewKeyStruct::errorSetModeString());

      if( errMode == FieldSpecification::SetErrorMode::error && m_isSurfaceGenerationCase )
      {
        errMode = FieldSpecification::SetErrorMode::warning;
      }

      switch( errMode )
      {
        case  FieldSpecification::SetErrorMode::silent:
          break;
        case  FieldSpecification::SetErrorMode::error:
          GEOS_THROW( message.str(), InputError, getDataContext() );
          break;
        case  FieldSpecification::SetErrorMode::warning:
          if( m_isSurfaceGenerationCase )
            message << "As the simulation includes a SurfaceGenerator, the set may be modified later";
          GEOS_WARNING( message.str() );
          break;
      }
    }

    if( isFieldNameFound == 0 )
    {
      std::ostringstream errorMessageBuilder;
      errorMessageBuilder << GEOS_FMT( "\n{}: there are no field named `{}` under the region `{}`.\n",
                                       fs.getWrapperDataContext( FieldSpecification::viewKeyStruct::fieldNameString() ),
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

      GEOS_THROW( errorMessageBuilder.str(), InputError, getDataContext() );
    }
  } );
}

void FieldSpecificationManager::applyInitialConditions( MeshLevel & mesh ) const
{
  this->forSubGroups< FieldSpecification >( [&] ( FieldSpecification const & fs )
  {
    if( fs.initialCondition() )
    {
      FieldSpecificationImpl::apply< dataRepository::Group >( fs,
                                                              mesh,
                                                              [&]( FieldSpecification const & bc,
                                                                   string const &,
                                                                   SortedArrayView< localIndex const > const & targetObject,
                                                                   Group & targetGroup,
                                                                   string const fieldName )
      {
        FieldSpecificationImpl::applyFieldValue< FieldSpecificationEqual >( bc, targetObject, 0.0, targetGroup, fieldName );
      } );
    }
  } );
}

}   /* namespace geos */
