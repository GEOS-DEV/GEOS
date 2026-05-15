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

#include "FieldSpecificationBase.hpp"

#include "common/format/StringUtilities.hpp"
#include "common/logger/Logger.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"

namespace geos
{
using namespace dataRepository;

FieldSpecificationBase::FieldSpecificationBase( string const & name, Group * parent ):
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

  registerWrapper( viewKeyStruct::regionNamesString(), &m_regionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the regions where the field specification is applied." );

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

  registerWrapper( viewKeyStruct::functionNameString(), &m_functionName ).
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


FieldSpecificationBase::~FieldSpecificationBase()
{}

FieldSpecificationBase::CatalogInterface::CatalogType &
FieldSpecificationBase::getCatalog()
{
  static FieldSpecificationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


void FieldSpecificationBase::postInputInitialization()
{
  GEOS_THROW_IF( !m_functionName.empty() &&
                 m_functionName.size() != 1 &&
                 m_functionName.size() != static_cast< string_array::size_type >( m_scale.size() ),
                 GEOS_FMT ( "Size mismatch: '{}' has {} entries but '{}' has {}. "
                            "'{}' either must be empty, have a single entry, or be sized exactly like '{}'",
                            viewKeyStruct::functionNameString(), m_functionName.size(),
                            viewKeyStruct::scaleString(), m_scale.size(),
                            viewKeyStruct::functionNameString(), viewKeyStruct::scaleString() ),
                 InputError,
                 getDataContext() );

  if( usesNonScalarValues() )
  {
    GEOS_THROW_IF( m_component != -1,
                   GEOS_FMT ( "'{}' must not be set when '{}' has more than one value.",
                              viewKeyStruct::componentString(),
                              viewKeyStruct::scaleString() ),
                   InputError,
                   getDataContext() );
  }

  if( !m_regionNames.empty() )
  {
    GEOS_THROW_IF( !m_objectPath.empty(),
                   GEOS_FMT ( "'{}' must not be set when '{}' is set.",
                              viewKeyStruct::objectPathString(),
                              viewKeyStruct::regionNamesString() ),
                   InputError,
                   getDataContext() );
  }
}

void FieldSpecificationBase::setMeshObjectPath( Group const & meshBodies )
{
  string const path = m_regionNames.empty()
                      ? m_objectPath
                      : "ElementRegions/{" + stringutilities::join( m_regionNames, ' ' ) + "}";
  try
  {
    m_meshObjectPaths = std::make_unique< MeshObjectPath >( path, meshBodies );
  }
  catch( std::exception const & e )
  {
    ErrorLogger::global().modifyCurrentExceptionMessage()
      .addToMsg( getWrapperDataContext( viewKeyStruct::objectPathString() ).toString() +
                 " is a wrong objectPath: " + path + "\n" )
      .addContextInfo( getWrapperDataContext( viewKeyStruct::objectPathString() ).getContextInfo()
                         .setPriority( 2 ) );
    throw InputError( e, getWrapperDataContext( viewKeyStruct::objectPathString() ).toString() +
                      " is a wrong objectPath: " + path + "\n" );
  }
}



REGISTER_CATALOG_ENTRY( FieldSpecificationBase, FieldSpecificationBase, string const &, Group * const )

}
