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

#include "FieldSpecificationABC.hpp"

namespace geos
{
using namespace dataRepository;

FieldSpecificationABC::FieldSpecificationABC( string const & name, Group * parent ):
  Group( name, parent )
{
  registerWrapper( viewKeyStruct::functionNamesString(), &m_functionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name(s) of function(s) that specifies variation of the boundary condition." );

  registerWrapper( viewKeyStruct::scalesString(), &m_scales ).
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

FieldSpecificationABC::~FieldSpecificationABC()
{}

FieldSpecificationABC::CatalogInterface::CatalogType &
FieldSpecificationABC::getCatalog()
{
  static FieldSpecificationABC::CatalogInterface::CatalogType catalog;
  return catalog;
}


void FieldSpecificationABC::postInputInitialization()
{
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
                           getCatalogName(), getName() ),
                 InputError,
                 getDataContext() );
}

}
