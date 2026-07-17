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
    setDescription( "Names of the regions that boundary condition is applied to." );

  registerWrapper( viewKeyStruct::fieldNameString(), &m_fieldName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of field that boundary condition is applied to.\n"
                    "A field can represent a physical variable. (pressure, temperature, global composition fraction of the fluid, ...)" );

  registerWrapper( viewKeyStruct::initialConditionString(), &m_initialCondition ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Boundary condition is applied as an initial condition." );

  registerWrapper( viewKeyStruct::functionNamesString(), &m_functionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name(s) of function(s) that specifies variation of the boundary condition." );

  registerWrapper( viewKeyStruct::scalesString(), &m_scales ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Apply scaling factor(s) for the value(s) of the boundary condition." );

  registerWrapper( viewKeyStruct::beginTimeString(), &m_beginTime ).
    setApplyDefaultValue( -1.0e99 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Time at which the boundary condition will start being applied." );

  registerWrapper( viewKeyStruct::endTimeString(), &m_endTime ).
    setApplyDefaultValue( 1.0e99 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Time at which the boundary condition will stop being applied." );
}


PermeabilitySpecification::~PermeabilitySpecification()
{}


void PermeabilitySpecification::postInputInitialization()
{
  for( real64 scale : getScales() )
  {
    GEOS_THROW_IF( scale < 0,
                   GEOS_FMT( "Scale values for a permeability must be non-negative\nA value of {} was given in {} '{}'.",
                             scale, catalogName(), getName() ),
                   InputError,
                   getDataContext() );
  }

  GEOS_THROW_IF( !m_functionNames.empty() &&
                  m_functionNames.size() != 1 &&
                  m_functionNames.size() != static_cast< string_array::size_type >( m_scales.size() ),
                  GEOS_FMT ( "Size mismatch: '{}' has {} entries but '{}' has {}. "
                            "'{}' either must be empty, have a single entry, or be sized exactly like '{}'",
                            viewKeyStruct::functionNamesString(), m_functionNames.size(),
                            viewKeyStruct::scalesString(), m_scales.size(),
                            viewKeyStruct::functionNamesString(), viewKeyStruct::scalesString() ),
                  InputError,
                  getDataContext() );

  GEOS_THROW_IF( m_beginTime > m_endTime,
                 GEOS_FMT( "{} ({}) must be less than {} ({}) in {} '{}'",
                           viewKeyStruct::beginTimeString(), m_beginTime,
                           viewKeyStruct::endTimeString(), m_endTime,
                           catalogName(), getName() ),
                 InputError,
                 getDataContext() );
}


REGISTER_CATALOG_ENTRY( FieldSpecificationABC, PermeabilitySpecification, string const &, Group * const )

template<>
void generateFieldSpecifications< PermeabilitySpecification >( PermeabilitySpecification const & ps,
                                                               dataRepository::Group & manager )
{
  for( string const & regionName : ps.getRegionNames() )
  {
    string const objectPath = "ElementRegions/" + regionName;

    string const childName = ps.getName() + "_" + regionName;

    FieldSpecification & fs = manager.registerGroup< FieldSpecification >( childName );
    fs.setFieldName( ps.getFieldName() );
    fs.setObjectPath( objectPath );
    fs.initialCondition( ps.initialCondition() );
    fs.setScales( ps.getScales() );
    fs.setFunctionNames( ps.getFunctionNames() );
    fs.setStartTime( ps.getStartTime() );
    fs.setEndTime( ps.getEndTime() );

    for( auto const & setName : ps.getSetNames() )
    {
      fs.addSetName( setName );
    }
  }
}

}
