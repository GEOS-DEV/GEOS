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

#include "PermeabilitySpecification.hpp"
#include "FieldSpecification.hpp"
#include "FieldSpecificationFactory.hpp"
#include "common/logger/Logger.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{
using namespace dataRepository;

PermeabilitySpecification::PermeabilitySpecification( string const & name, Group * parent ):
  FieldSpecificationABC( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::setNamesString(), &m_setNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of sets that the boundary condition is applied to.\n"
                    "A set can contain heterogeneous elements in the mesh (volumes, nodes, faces, edges).\n"
                    "A set can be be defined by a 'Geometry' component, or correspond to imported sets in case of an external mesh" );

  registerWrapper( viewKeyStruct::regionNamesString(), &m_regionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Names of the cell regions that boundary condition is applied to." );

  registerWrapper( viewKeyStruct::fieldNameString(), &m_fieldName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of field that boundary condition is applied to.\n"
                    "A field can represent a physical variable. (pressure, temperature, global composition fraction of the fluid, ...)" );

  registerWrapper( viewKeyStruct::componentString(), &m_component ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component of field (if tensor) to apply boundary condition to.\n"
                    "The component must use the order in which the phaseNames have been defined in the Constitutive Element." );

  getWrapper< int >( viewKeyStruct::initialConditionString() ).
    setApplyDefaultValue( 1 );
}


PermeabilitySpecification::~PermeabilitySpecification()
{}


void PermeabilitySpecification::postInputInitialization()
{
  FieldSpecificationABC::postInputInitialization();

  for( real64 scale : getScales() )
  {
    GEOS_THROW_IF( scale < 0,
                   GEOS_FMT( "Scale values for a permeability must be non-negative\nA value of {} was given in {} '{}'.",
                             scale, catalogName(), getName() ),
                   InputError,
                   getDataContext() );
  }

  GEOS_THROW_IF( m_component != -1 && m_scales.size() > 1,
                 GEOS_FMT ( "'{}' must not be set when '{}' has more than one value.",
                            viewKeyStruct::componentString(),
                            viewKeyStruct::scalesString() ),
                 InputError,
                 getDataContext() );
}

namespace
{

/**
 * @brief @return The element region named @p regionName in the domain if it exists, or a nullptr
 */
ElementRegionBase const * findElementRegion( DomainPartition const & domain,
                                             string const & regionName )
{
  ElementRegionBase const * region = nullptr;
  domain.forMeshBodies( [&]( MeshBody const & meshBody )
  {
    meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
    {
      if( region == nullptr && meshLevel.getElemManager().hasRegion( regionName ) )
      {
        region = &meshLevel.getElemManager().getRegion( regionName );
      }
    } );
  } );
  return region;
}

/**
 * @brief Validate that a region with the name @p regionName exists and is a valid cell region
 * @param domain The domain
 * @param regionName The name of the region to validated
 * @param permSpec Reference to the current object to print its DataContext
 * @note Throws if the region doesn't exists or isn't a cell region
 */
void expectValidCellRegion( DomainPartition const & domain,
                            string const & fullRegionName,
                            PermeabilitySpecification const & permSpec )
{
  // get only region from "region/subregion"
  string const regionName = fullRegionName.substr( 0, fullRegionName.find( '/' ) );

  ElementRegionBase const * const region = findElementRegion( domain, regionName );

  GEOS_THROW_IF( region == nullptr,
                 GEOS_FMT( "Region '{}' does not exist.", regionName ),
                 InputError,
                 permSpec.getDataContext() );

  GEOS_THROW_IF( dynamic_cast< CellElementRegion const * >( region ) == nullptr,
                 GEOS_FMT( "Region '{}' is a '{}', but must be a '{}'",
                           regionName,
                           region->getCatalogName(),
                           CellElementRegion::catalogName() ),
                 InputError,
                 permSpec.getDataContext() );
}

}


template<>
void expandFieldSpecification< PermeabilitySpecification >( PermeabilitySpecification const & permSpec,
                                                            dataRepository::Group & manager )
{
  Group const & problem = manager.getParent();
  DomainPartition const & domain = problem.getGroup< DomainPartition >( "domain" );

  for( string const & regionName : permSpec.getRegionNames() )
  {
    expectValidCellRegion( domain, regionName, permSpec );

    string const objectPath = "ElementRegions/" + regionName;
    string const childName = permSpec.getName() + "_" + regionName;

    FieldSpecification & fs = manager.registerGroup< FieldSpecification >( childName );
    fs.setDataContextReference( permSpec );
    fs.setFieldName( permSpec.getFieldName() );
    fs.setComponent( permSpec.getComponent() );
    fs.setObjectPath( objectPath );
    fs.initialCondition( permSpec.initialCondition() );
    fs.setScales( permSpec.getScales() );
    fs.setFunctionNames( permSpec.getFunctionNames() );
    fs.setStartTime( permSpec.getStartTime() );
    fs.setEndTime( permSpec.getEndTime() );
    fs.setErrorSetMode( permSpec.getErrorSetMode() );

    for( auto const & setName : permSpec.getSetNames() )
    {
      fs.addSetName( setName );
    }
  }
}

REGISTER_CATALOG_ENTRY( FieldSpecificationABC, PermeabilitySpecification, string const &, Group * const )
REGISTER_FIELD_SPECIFICATION_PROCESSOR( PermeabilitySpecification )

}
