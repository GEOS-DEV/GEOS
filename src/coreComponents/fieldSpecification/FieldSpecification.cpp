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

#include "FieldSpecification.hpp"

namespace geos
{
using namespace dataRepository;

FieldSpecification::FieldSpecification( string const & name, Group * parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::setNamesString(), &m_setNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of sets that the boundary condition is applied to.\n"
                    "A set can contain heterogeneous elements in the mesh (volumes, nodes, faces, edges).\n"
                    "A set can be be defined by a 'Geometry' component, or correspond to imported sets in case of an external mesh" );

  registerWrapper( viewKeyStruct::objectPathString(), &m_objectPath ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Path to the target field" );

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

  registerWrapper( viewKeyStruct::directionString(), &m_direction ).
    setApplyDefaultValue( {0, 0, 0} ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Direction to apply boundary condition to." );

  registerWrapper( viewKeyStruct::functionNamesString(), &m_functionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name(s) of function(s) that specifies variation of the boundary condition." );

  registerWrapper( viewKeyStruct::scaleString(), &m_scale ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Apply scaling factor(s) for the value(s) of the boundary condition." );

  registerWrapper( viewKeyStruct::initialConditionString(), &m_initialCondition ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Boundary condition is applied as an initial condition." );

  registerWrapper( viewKeyStruct::beginTimeString(), &m_beginTime ).
    setApplyDefaultValue( -1.0e99 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Time at which the boundary condition will start being applied." );

  registerWrapper( viewKeyStruct::endTimeString(), &m_endTime ).
    setApplyDefaultValue( 1.0e99 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Time at which the boundary condition will stop being applied." );

  registerWrapper( viewKeyStruct::errorSetModeString(), &m_emptySetErrorMode ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( SetErrorMode::error ).
    setDescription( GEOS_FMT( "Set the log state when a “set” does not target any region\n"
                              "When set to \"{}\", no output.\n"
                              "When set to \"{}\", output a warning.\n"
                              "When set to \"{}\", output a throw.\n",
                              EnumStrings< SetErrorMode >::toString( SetErrorMode::silent ),
                              EnumStrings< SetErrorMode >::toString( SetErrorMode::warning ),
                              EnumStrings< SetErrorMode >::toString( SetErrorMode::error )  ));
}


FieldSpecification::~FieldSpecification()
{}

FieldSpecification::CatalogInterface::CatalogType &
FieldSpecification::getCatalog()
{
  static FieldSpecification::CatalogInterface::CatalogType catalog;
  return catalog;
}


void FieldSpecification::postInputInitialization()
{
  { // both conditions work together
    GEOS_THROW_IF( !m_functionNames.empty() &&
                   m_functionNames.size() != 1 &&
                   m_functionNames.size() != static_cast< string_array::size_type >( m_scale.size() ),
                   GEOS_FMT ( "Size mismatch: '{}' has {} entries but '{}' has {}. "
                              "'{}' either must be empty, have a single entry, or be sized exactly like '{}'",
                              viewKeyStruct::functionNamesString(), m_functionNames.size(),
                              viewKeyStruct::scaleString(), m_scale.size(),
                              viewKeyStruct::functionNamesString(), viewKeyStruct::scaleString() ),
                   InputError, getDataContext() );

                              viewKeyStruct::scalesString(), m_scales.size(),
                              viewKeyStruct::functionNamesString(), viewKeyStruct::scalesString() ),
                   InputError,
                   getDataContext() );

    GEOS_THROW_IF( m_component != -1 && m_scale.size() > 1,
                   GEOS_FMT ( "'{}' must not be set when '{}' has more than one value.",
                              viewKeyStruct::componentString(),
                              viewKeyStruct::scaleString() ),
                   InputError, getDataContext() );
  }
}

void FieldSpecification::validateNumArrayComp( localIndex numComp )
{
  if( m_component != -1 )
  {
    return;
  }

  auto expand = [&]( auto & values, string_view attributeName )
  {
    GEOS_THROW_IF( values.size() > 1 && static_cast< localIndex >( values.size() ) != numComp,
                   GEOS_FMT ( "'{}' has {} entries but the target field '{}' has {} components. "
                              "'{}' must either have a single entry (applied to every components) "
                              "or exactly have {} entries.",
                              attributeName, values.size(), m_fieldName, numComp,
                              attributeName, numComp ),
                   InputError,
                   getDataContext() );

    // "Broadcast"/duplicate values
    if( values.size() == 1 && numComp > 1 )
    {
      auto const value = values[ 0 ];
      values.resize( numComp );
      for( localIndex i = 0; i < numComp; ++i )
      {
        values[ i ] = value;
      }
    }
  };

  expand( m_scale, viewKeyStruct::scaleString() );
  if( !m_functionNames.empty() )
  {
    expand( m_functionNames, viewKeyStruct::functionNamesString() );
  }
}

void FieldSpecification::setMeshObjectPath( Group const & meshBodies )
{
  try
  {
    m_meshObjectPaths = std::make_unique< MeshObjectPath >( m_objectPath, meshBodies );
  }
  catch( std::exception const & e )
  {
    ErrorLogger::global().modifyCurrentExceptionMessage()
      .addToMsg( getWrapperDataContext( viewKeyStruct::objectPathString() ).toString() +
                 " is a wrong objectPath: " + m_objectPath + "\n" )
      .addContextInfo( getWrapperDataContext( viewKeyStruct::objectPathString() ).getContextInfo()
                         .setPriority( 2 ) );
    throw InputError( e, getWrapperDataContext( viewKeyStruct::objectPathString() ).toString() +
                      " is a wrong objectPath: " + m_objectPath + "\n" );
  }
}



REGISTER_CATALOG_ENTRY( FieldSpecification, FieldSpecification, string const &, Group * const )

}
